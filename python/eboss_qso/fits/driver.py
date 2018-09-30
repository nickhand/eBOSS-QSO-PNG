from jinja2 import Environment, FileSystemLoader
import os
import numpy
from eboss_qso.fits.preparer import QSOFitPreparer
from pyRSD.rsd import QuasarSpectrum
from pyRSD.rsdfit.data import PoleCovarianceMatrix

import tempfile
import errno
import shutil
import numpy as np
from mpi4py import MPI


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class QSOFitDriver(object):
    """
    Class to drive a QSO power spectrum fit.
    """

    def __init__(self, rsdfit_args,
                 kind,
                 version,
                 sample,
                 ells,
                 z_eff,
                 vary,
                 theory_p,
                 data_p=None,
                 kmin=0.0001,
                 kmax=0.4,
                 cov_type='analytic',
                 error_rescale=1.0,
                 overwrite=False,
                 output_only=False,
                 comm=None,
                 use_temp_files=False,
                 tag=None,
                 box=None):

        quiet = False
        if output_only:
            overwrite = False
            quiet = True

        # determine the communicator
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # save the input args
        self.tag = tag
        self.rsdfit_args = rsdfit_args
        self.vary = vary

        # prepare the fit
        assert kmin >= 0.0001
        assert kmax <= 0.4
        self.prep = QSOFitPreparer(kind,
                                   version,
                                   sample,
                                   ells,
                                   p=data_p,
                                   z_eff=z_eff,
                                   cov_type=cov_type,
                                   error_rescale=error_rescale,
                                   kmin=0.0001,
                                   kmax=0.4,
                                   overwrite=overwrite,
                                   quiet=quiet,
                                   use_temp_files=use_temp_files,
                                   box=box)
        self.config = self.prep.config

        # check the p
        self.p = theory_p
        if self.prep.p is not None:
            if self.prep.p != self.p:
                args = (self.prep.p, self.p)
                raise ValueError("data p = %.1f, theory p = %.1f" % args)

        # log the p
        args = (str(self.prep.p), str(self.p))
        print("data p = %s, theory p = %s" % args)

        # store zmin/zmax
        self.zmin = self.prep.attrs['zmin']
        self.zmax = self.prep.attrs['zmax']

        if output_only:
            print(self.output_dir)
            return

        # make the directory output
        mkdir_p(self.output_dir)

        # keywords we are going to add to parameter file template
        kws = {}
        kws['kmin'] = 1e-6
        kws['kmax'] = 1.0
        kws['covariance_file'] = self.prep.covariance_file
        kws['data_file'] = self.prep.data_file
        kws['ells'] = self.prep.ells
        kws['stats'] = self.prep.stat_names
        kws['fitting_range'] = [(kmin, kmax) for _ in ells]
        kws['max_ellprime'] = 4
        kws['z_eff'] = self.prep.z_eff[0]

        # do the window function dict
        kws['window_file'] = {}
        for stat, window_file in zip(self.prep.stat_names, self.prep.window_file):
            kws['window_file'][stat] = window_file

        # add stat specific params
        d = {}
        if not np.allclose(kws['z_eff'], self.prep.z_eff):
            for stat_name, z_eff in zip(self.prep.stat_names, self.prep.z_eff):
                d[stat_name] = {'z': z_eff}
        kws['stat_specific_params'] = d

        # how to initialize?
        if 'mcmc' in ' '.join(self.rsdfit_args):
            kws['init_from'] = 'nlopt'
            kws['lbfgs_numerical_from_lnlike'] = False
        else:
            kws['init_from'] = 'fiducial'
            kws['lbfgs_numerical_from_lnlike'] = True

        # re-scale mock covariance?
        if self.prep.kind == 'data':
            kws['covariance_Nmocks'] = 0
        else:
            C = PoleCovarianceMatrix.from_plaintext(kws['covariance_file'])
            kws['covariance_Nmocks'] = C.attrs.get('Nmock', 0)

        # render the parameter file
        params = self._render_params(**kws)

        # make the parameter file
        with tempfile.NamedTemporaryFile(mode='wb') as ff:

            # write out the rendered template
            ff.write((params+"\n\n").encode())

            # write out the theory too
            model = QuasarSpectrum(z=self.prep.z_eff[0])
            theorypars = model.default_params()

            # the pars which we are varying
            for par in theorypars:
                theorypars[par].vary = par in self.vary

            # update b1 bounds
            theorypars['b1'].update(lower=0.1, upper=6.0)

            # update p
            theorypars['p'].update(vary=False, value=self.p, fiducial=self.p)

            # update redshift-dependent quantities
            for par in ['f', 'sigma8_z']:
                value = getattr(model, par)
                theorypars[par].update(
                    value=value, fiducial=value, lower=0, upper=2.0)

            for par in ['alpha_par', 'alpha_perp']:
                theorypars[par].update(lower=0.3, upper=1.8)

            theorypars['N'].update(lower=-5000, upper=5000)

            # write to file
            theorypars.to_file(ff, mode='a')

            # run
            ff.seek(0)
            self._run(ff.name)

    @classmethod
    def run_from_args(cls, args=None, comm=None):
        """
        Initialize the QSOFitDriver from command-line arguments.
        """
        import argparse

        descr = 'run a QSO power spectrum fit'
        parser = argparse.ArgumentParser(description=descr)

        h = 'the kind of measurement we are fitting'
        parser.add_argument('--kind', choices=['data', 'ezmock'],
                            type=str, required=True, help=h)

        h = 'the version of the data we are fitting'
        parser.add_argument('--version', type=str,
                            required=True, help=h)

        h = 'fit either the north or the south sample'
        parser.add_argument('--sample', choices=['N', 'S'],
                            type=str, required=True, help=h)

        h = 'the multipoles to fit'
        parser.add_argument('--ells', nargs="+", type=int,
                            help=h, required=True)

        h = 'the effective redshift'
        parser.add_argument('--z_eff', type=float, nargs="+",
                            required=True, help=h)

        h = 'the p value used for the theory'
        parser.add_argument('--theory_p', choices=[1.0, 1.6],
                            required=True, type=float, help=h)

        h = 'the p value used for the data'
        parser.add_argument('--data_p', choices=[1.0, 1.6],
                            default=None, type=float, help=h)

        h = 'the minimum k value to include'
        parser.add_argument('--kmin', type=float, default=0.0001, help=h)

        h = 'the maximum k value to include'
        parser.add_argument('--kmax', type=float, default=0.4, help=h)

        h = 'the type of covariance to use'
        parser.add_argument('--cov', dest='cov_type',
                            choices=['analytic', 'mock'], default='analytic', help=h)

        h = 'tag the output directory'
        parser.add_argument('--tag', type=str, default=None, help=h)

        h = 'rescale the errors by this amount'
        parser.add_argument('--error-rescale', type=float, default=1.0, help=h)

        h = 'the parameters to vary'
        choices = ['alpha_par', 'alpha_perp', 'f',
                   'sigma8_z', 'b1', 'sigma_fog', 'f_nl', 'N']
        parser.add_argument('--vary', type=str, nargs='+',
                            choices=choices, help=h, required=True)

        h = 'whether to overwrite existing files'
        parser.add_argument('--overwrite', action='store_true', help=h)

        h = 'whether to only print the output directory'
        parser.add_argument('--output', dest='output_only',
                            action='store_true', help=h)

        h = 'whether to use temporary files'
        parser.add_argument('--use-temp-files', action='store_true', help=h)

        ns, unknown = parser.parse_known_args(args=args)
        return cls(rsdfit_args=unknown, comm=comm, **vars(ns))

    @property
    def hashinfo(self):
        """
        The meta-data that will be used to make the hash string.
        """
        # make the hash string
        meta = {}
        meta['vary'] = self.vary
        meta['stats'] = self.prep.stats
        meta['hashstr'] = self.prep.hashstr
        if 'f_nl' in self.vary:
            meta['p'] = self.p
        else:
            meta['p'] = None

        return meta

    @property
    def output_dir(self):
        """
        The output directory name.

        Format is krange/params/zbounds
        """
        try:
            return self._output_dir
        except AttributeError:
            from eboss_qso.measurements.utils import make_hash

            usekeys = list(self.hashinfo.keys())
            if getattr(self.prep, 'box', None) is not None:
                usekeys.pop(usekeys.index('spectra_file'))
            hashstr = make_hash(self.hashinfo, usekeys=usekeys)

            # k-range
            path = "%s-%s" % (self.prep.kmin, self.prep.kmax)

            # params path
            tmp = "basemodel"
            if 'f' in self.vary and 'sigma8_z' in self.vary:
                tmp += "-fs8"
            if 'alpha_par' in self.vary or 'alpha_perp' in self.vary:
                tmp += '-alphas'
            if 'N' in self.vary:
                tmp += '-N'
            if 'f_nl' in self.vary:
                tmp += '-fnl'
            path = os.path.join(path, tmp)

            # z bounds
            path = os.path.join(path, "%.1f-%.1f" % (self.zmin, self.zmax))

            # the output directory name
            sample = self.config.sample
            stats = '+'.join(self.prep.stats)
            box = getattr(self.prep, 'box', None)

            cov_type = self.prep.cov_type+'-cov'
            tag = f'QSO-{sample}-'
            if box is not None:
                tag += box + '-'
            tag += f'{stats}-{cov_type}-{hashstr}'

            # full path
            self._output_dir = os.path.join(
                self.config.fits_results_dir, path, tag)
            if self.tag is not None:
                self._output_dir += '_' + self.tag
            return self._output_dir

    def _run(self, param_file):
        """
        Internal function that will call the ``rsdfit`` command.
        """
        from pyRSD.rsdfit.util import rsdfit_parser
        from pyRSD.rsdfit import rsdfit
        from pyRSD.rsdfit.parameters import ParameterSet

        # the arguments to pass to RSDFit
        args = self.rsdfit_args
        args += ['-p', param_file, '-o', self.output_dir, '--no-save-model']

        # parse the args passed to rsdfit
        args = rsdfit_parser().parse_args(args)
        args = vars(args)

        # initialize and run the RSDFit driver
        mode = args.pop('subparser_name')
        driver = rsdfit.RSDFitDriver(self.comm, mode, **args)
        driver.run()

        # and save the hashkey info
        if os.path.isdir(self.output_dir):
            with open(os.path.join(self.output_dir, 'hashinfo.json'), 'w') as ff:
                import json
                json.dump(self.hashinfo, ff)

        # remove temporary files
        if self.prep.use_temp_files:

            # read the original parameter file
            params_file = os.path.join(self.output_dir, 'params.dat')
            lines = open(params_file, 'r').readlines()

            # rename the input files
            for name in ['data_file', 'covariance_file', 'window_file']:
                f = getattr(self.prep, name)
                newf = getattr(self.prep, '_' + name)

                if not isinstance(f, list):
                    f = [f]
                if not isinstance(newf, list):
                    newf = [newf]

                # rename the file
                for a, b in zip(f, newf):
                    if os.path.exists(a):
                        shutil.move(a, b)

            # search and fix NERSC-specific paths
            for i, line in enumerate(lines):
                tag = None
                if line.startswith('data.covariance ='):
                    tag = "covariance"
                    newf = self.prep._covariance_file
                elif line.startswith('data.window_file ='):
                    tag = 'window_file'
                    newf = self.prep._window_file
                elif line.startswith('data.data_file ='):
                    tag = 'data_file'
                    newf = self.prep._data_file

                # replace line
                if tag is not None:
                    if not isinstance(newf, list):
                        newf = [newf]

                    for a in newf:
                        newf = os.path.join('$(EBOSS_DIR)', os.path.relpath(
                            a, os.environ['EBOSS_DIR']))
                        lines[i] = lines[i].replace(a, newf)

            # write out new parameter file
            with open(params_file, 'w') as ff:
                ff.write("".join(lines))

    def _render_params(self, **kwargs):
        """
        Return the rendered parameter file.
        """
        # jinja environ
        jinja_env = Environment(
            loader=FileSystemLoader(self.config.fits_params_dir))
        tpl = jinja_env.get_template('template.params')
        return tpl.render(**kwargs)


def __main__():
    driver = QSOFitDriver.run_from_args()
