import os
import numpy
from nbodykit.lab import ConvolvedFFTPower
import tempfile

from ..measurements.results import info_from_filename
from ..measurements.utils import find_window_measurement
from ..measurements import get_hashkeys
from ..measurements.zweights import bias_model

from . import eBOSSConfig
from .data import write_data
from .window import compute_window
from .covariance import compute_analytic_covariance, compute_ezmock_covariance

def get_spectra_type(filename):
    """
    Return the spectrum type, i.e., 'data', 'ezmock', etc
    """
    if 'data' in filename:
        return 'data'
    elif 'mocks' in filename:
        home_dir = os.environ['EBOSS_DIR']
        spectra_dir = os.path.join(home_dir, 'measurements', 'spectra', 'mocks')
        return os.path.relpath(filename, spectra_dir).split(os.path.sep)[0]
    else:
        raise ValueError("do not understand file name '%s'" %filename)


class QSOFitPreparer(object):
    """
    Class to prepare a QSO power spectrum fit.
    """
    def __init__(self, spectra_file, stats, kmin=0.0001, kmax=0.4, cov_type='analytic',
                    overwrite=False, error_rescale=1.0, quiet=False, use_temp_files=False):

        self.spectra_file = os.path.abspath(spectra_file)
        self.stats = stats
        self.kmin = kmin
        self.kmax = kmax
        self.overwrite = overwrite
        self.quiet = quiet
        self.error_rescale = error_rescale
        self.use_temp_files = use_temp_files
        self.cov_type = cov_type

        # dictionary of info from the filename
        self.kind = get_spectra_type(self.spectra_file)
        info = info_from_filename(self.kind, os.path.basename(self.spectra_file))

        # make sure we have all of the keys
        assert 'version' in info
        assert 'hashstr' in info
        assert 'sample' in info and info['sample'] in 'NS'

        # copy over file info
        for k in info:
            setattr(self, k, info[k])

        # check statistics
        self.ells = []
        for stat in self.stats:
            if stat == 'P0_sysfree':
                self.ells.append([0,2])
            else:
                self.ells.append(int(stat[-1]))

        # the input for the hash string
        self.hashinput = get_hashkeys(self.spectra_file, 'ConvolvedFFTPower')
        assert all(k in self.hashinput for k in ['zmin', 'zmax'])

        # the configuration
        self.config = eBOSSConfig(self.sample, self.version, self.kind)

        # determine the statistic names
        tag = 'ngc' if self.config.sample == 'N' else 'sgc'
        self.stat_names = ['pole_' + tag + '_%s' % stat.split('_')[0][1] for stat in stats]

        # load z_eff and nbar_eff
        res = ConvolvedFFTPower.load(self.spectra_file)
        self.z_eff = res.attrs.get('z_eff', None)
        self.nbar_eff = res.attrs.get('nbar_eff', None)

        if self.z_eff is None:
            raise ValueError("need 'z_eff' in 'attrs' of spectra result")
        if self.nbar_eff is None:
            raise ValueError("need 'nbar_eff' in 'attrs' of spectra result")

        # and run
        self.write_data()
        self.write_window()
        self.write_covariance()

    @classmethod
    def initialize(cls):
        import argparse

        descr = 'prepare a QSO power spectrum fit'
        parser = argparse.ArgumentParser(description=descr)

        h = 'the power spectrum file we wish to fit'
        parser.add_argument('-f', '--spectra_file', type=str, help=h, required=True)

        h = 'the minimum k value to include'
        parser.add_argument('--kmin', type=float, default=0.0001, help=h)

        h = 'the maximum k value to include'
        parser.add_argument('--kmax', type=float, default=0.4, help=h)

        h = 'the statistics to include'
        stats = ['P0', 'P2', 'P0_sysfree']
        parser.add_argument('--stats', nargs='*', type=str, choices=stats, default=['P0', 'P2'], help=h)

        h = 'whether to overwrite existing files'
        parser.add_argument('--overwrite', action='store_true', help=h)

        h = 'whether to use temporary files'
        parser.add_argument('--use-temp-files', action='store_true', help=h)

        h = 'rescale the errors by this amount'
        parser.add_argument('--error-rescale', type=float, default=1.0, help=h)

        h = 'the type of covariance to use'
        parser.add_argument('--cov', dest='cov_type', choices=['analytic', 'mock'],
                                default='analytic', help=h)

        ns = parser.parse_args()
        return cls(**vars(ns))

    def write_data(self):
        """
        Write the necessary data file.
        """
        # the output data file
        stats = '+'.join(self.stats)
        filename = f"poles_{self.version}-QSO-{self.sample}"
        box = getattr(self, 'box', None)
        if box is not None:
            filename += '-' + box
        filename += f"_{stats}_{self.hashstr}.dat"
        output = os.path.join(self.config.fits_data_dir, filename)

        if self.use_temp_files:
            self._data_file = output
            output = tempfile.mktemp()

        # make the data file
        if not os.path.exists(output) or self.overwrite:
            write_data(self.spectra_file, self.stats, self.stat_names, output,
                        kmin=self.kmin, kmax=self.kmax, quiet=self.quiet)
        else:
            if not self.quiet:
                print('skipping data preparation...')
        self.data_file = output

    def write_window(self):
        """
        Write the necessary window file.
        """
        from ..measurements.utils import make_hash

        # the window file
        meta = {'zmin':self.hashinput['zmin'], 'zmax':self.hashinput['zmax']}
        hashstr = make_hash(meta)
        filename = f"poles_{self.version}-QSO-{self.sample}_{hashstr}.dat"
        output = os.path.join(self.config.fits_window_dir, filename)

        if self.use_temp_files:
            self._window_file = output
            output = tempfile.mktemp()

        # make the window
        if not os.path.exists(output) or self.overwrite:

            # the name of the (unformatted) window file
            if self.kind == 'ezmock':
                version = self.version[:4]
            elif self.kind == 'data':
                version = self.version
            else:
                raise ValueError("do not understand 'kind' = '%s'" %self.kind)
            window_file = find_window_measurement(version, self.sample,
                                                    self.hashinput['zmin'],
                                                    self.hashinput['zmax'],
                                                    self.hashinput['p'])

            print('using window file %s...' % window_file)
            ells = [0,2,4,6,8,10]
            compute_window(window_file, ells, output, smin=1e-2, smax=1e4, quiet=self.quiet)
        else:
            if not self.quiet:
                print('skipping window preparation...')
        self.window_file = output

    def write_covariance(self):
        """
        Write the necessary covariance file.
        """
        # the output data file
        stats = '+'.join(self.stats)
        filename = f"poles_{self.version}-QSO-{self.sample}-{self.cov_type}"
        if self.cov_type == 'analytic':
            box = getattr(self, 'box', None)
            if box is not None and box == 'mean':
                filename += '-mean'
        filename += f"_{stats}_{self.hashstr}.dat"
        output = os.path.join(self.config.fits_covariance_dir, filename)

        if self.use_temp_files:
            self._covariance_file = output
            output = tempfile.mktemp()

        # make the covariance
        if not os.path.exists(output) or self.overwrite:

            if self.cov_type == 'analytic':
                from pyRSD.rsd import QuasarSpectrum

                P0_FKP = self.hashinput.get('P0_FKP', None)
                if P0_FKP is None: P0_FKP = 0.

                # load the model
                model = QuasarSpectrum(z=self.z_eff, params=config.cosmo)

                # compute
                C = compute_analytic_covariance(self.config, self.stats, model,
                                                self.hashinput['zmin'], self.hashinput['zmax'],
                                                P0_FKP, kmin=self.kmin, kmax=self.kmax, dk=0.005)

                C *= self.error_rescale

                # and save
                if not self.quiet:
                    print("saving %s..." %output)
                C.attrs.clear()
                C.to_plaintext(output)

            else:
                # compute the covariance from the mocks
                C = compute_ezmock_covariance('v1.8e-fph', self.sample, self.stats,
                                                kmin=self.kmin, kmax=self.kmax, p=self.hashinput['p'])

                # rescale
                C *= self.error_rescale

                # and save
                C.to_plaintext(output)

        else:
            if not self.quiet:
                print('skipping covariance preparation...')
        self.covariance_file = output
