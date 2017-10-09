from jinja2 import Environment, FileSystemLoader
import os
import numpy
from eboss_qso.fits.preparer import QSOFitPreparer
from pyRSD.rsd import QuasarSpectrum
import tempfile

class QSOFitDriver(object):
    """
    Class to drive a QSO power spectrum fit.

    Parameters
    ----------
    rsdfit_args : list
        list of arguments to pass to ``rsdfit`` command, i.e., "mcmc", etc
    spectra_file : str
        the data power spectrum file
    vary : list of str
        list of parameter names we want to vary in the fit, i.e., ``f``
    stats : list of str
        list of the names of the statistics we want to fit
    kmin : float, optional
        the minimum ``k`` value to include in the fit
    kmax : float, optional
        the maximum  ``k`` value to include in the fit
    overwrite : bool, optional
        if ``True`` overwrite the input fit files
    """
    def __init__(self, rsdfit_args, spectra_file, vary, stats,
                    kmin=0.0001, kmax=0.4, overwrite=False, output_only=False):

        quiet = False
        if output_only:
            overwrite = False
            quiet = True

        self.rsdfit_args = rsdfit_args
        self.vary = vary
        self.stats = stats

        # prepare the fit
        assert kmin >= 0.0001
        assert kmax <= 0.4
        self.preparer = QSOFitPreparer(spectra_file, stats, kmin=0.0001, kmax=0.4, overwrite=overwrite, quiet=quiet)
        self.config = self.preparer.config

        if output_only:
            print(self.output_dir)
            return

        # keywords we are going to add to parameter file template
        kws = {}
        kws['kmin'] = 1e-6
        kws['kmax'] = 1.0
        kws['covariance_file'] = self.preparer.covariance_file
        kws['data_file'] = self.preparer.data_file
        kws['window_file'] = self.preparer.window_file
        kws['ells'] = self.preparer.ells
        kws['stats'] = self.preparer.stat_names
        kws['z_eff'] = self.preparer.z_eff
        kws['fitting_range'] = [(kmin, kmax) for stat in stats]
        kws['max_ellprime'] = 4
        kws['theory_decorator'] = {}
        for i, stat in enumerate(stats):
            if stat == 'P0_sysfree':
                kws['theory_decorator'][kws['stats'][i]] = "systematic_free_P0"
        if 'mcmc' in ' '.join(self.rsdfit_args):
            kws['init_from'] = 'nlopt'
        else:
            kws['init_from'] = 'fiducial'

        # render the parameter file
        params = self._render_params(**kws)

        # make the parameter file
        with tempfile.NamedTemporaryFile(mode='wb') as ff:

            # write out the rendered template
            ff.write((params+"\n\n").encode())

            # write out the theory too
            model = QuasarSpectrum(z=self.preparer.z_eff)
            theorypars = model.default_params()

            # the pars which we are varying
            for par in theorypars:
                theorypars[par].vary = par in self.vary

            # update redshift-dependent quantities
            for par in ['f', 'sigma8_z']:
                value = getattr(model, par)
                theorypars[par].update(value=value, fiducial=value, lower=0.5*value, upper=1.5*value)

            for par in ['alpha_par', 'alpha_perp']:
                theorypars[par].update(lower=0.6, upper=1.4)

            # write to file
            theorypars.to_file(ff, mode='a')

            # run
            ff.seek(0)
            self._run(ff.name)


    @classmethod
    def initialize(cls):
        """
        Initialize the QSOFitDriver from command-line arguments.
        """
        import argparse

        descr = 'run a QSO power spectrum fit'
        parser = argparse.ArgumentParser(description=descr)

        h = 'the power spectrum file we wish to fit'
        parser.add_argument('-f', '--spectra_file', type=str, help=h, required=True)

        h = 'the minimum k value to include'
        parser.add_argument('--kmin', type=float, default=0.0001, help=h)

        h = 'the maximum k value to include'
        parser.add_argument('--kmax', type=float, default=0.4, help=h)

        h = 'the parameters to vary'
        choices = ['alpha_par', 'alpha_perp', 'f', 'sigma8_z', 'b1', 'sigma_fog', 'f_nl']
        parser.add_argument('--vary', type=str, nargs='+', choices=choices, help=h, required=True)

        h = 'the statistics to include'
        stats = ['P0', 'P2', 'P0_sysfree']
        parser.add_argument('--stats', nargs='*', type=str, choices=stats, default=['P0', 'P2'], help=h)

        h = 'whether to overwrite existing files'
        parser.add_argument('--overwrite', action='store_true', help=h)

        h = 'whether to only print the output directory'
        parser.add_argument('--output', dest='output_only', action='store_true', help=h)

        ns, unknown = parser.parse_known_args()
        return cls(rsdfit_args=unknown, **vars(ns))

    @property
    def hashinfo(self):
        """
        The meta-data that will be used to make the hash string.
        """
        # make the hash string
        meta = {}
        meta['vary'] = self.vary
        meta['spectra_file'] = self.preparer.spectra_file
        meta['stats'] = self.stats
        return meta

    @property
    def output_dir(self):
        """
        The output directory name.
        """
        try:
            return self._output_dir
        except AttributeError:
            from eboss_qso.measurements.utils import make_hash
            hashstr = make_hash(self.hashinfo)

            # the output directory name
            sample = self.config.sample
            stats = '+'.join(self.stats)
            tag = f'QSO-{sample}-{stats}-{hashstr}'
            self._output_dir = os.path.join(self.config.fits_results_dir, tag)
            return self._output_dir

    def _run(self, param_file):
        """
        Internal function that will call the ``rsdfit`` command.
        """
        # the arguments to pass to RSDFit
        args = self.rsdfit_args
        args += ['-p', param_file, '-o', self.output_dir]

        # run RSDFit
        cmd = 'rsdfit' + ' ' + ' '.join(args)
        print(f"calling '{cmd}'...")
        ret = os.system(cmd)

        # and save the hashkey info
        if os.path.isdir(self.output_dir):
            with open(os.path.join(self.output_dir, 'hashinfo.json'), 'w') as ff:
                import json
                json.dump(self.hashinfo, ff)


    def _render_params(self, **kwargs):
        """
        Return the rendered parameter file.
        """
        # jinja environ
        jinja_env = Environment(loader=FileSystemLoader(self.config.fits_params_dir))
        tpl = jinja_env.get_template('template.params')
        return tpl.render(**kwargs)


def __main__():
    driver = QSOFitDriver.initialize()
