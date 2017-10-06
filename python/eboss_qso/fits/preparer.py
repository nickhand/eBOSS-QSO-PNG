import os
import numpy
from nbodykit.lab import ConvolvedFFTPower

from ..measurements.results import info_from_filename
from ..measurements import get_hashkeys
from ..measurements.zweights import bias_model

from . import eBOSSConfig
from .data import write_data
from .window import compute_window
from .covariance import compute_covariance

def get_spectra_type(filename):
    """
    Return the spectrum type, i.e., 'data', 'ezmock', etc
    """
    if 'data' in filename:
        return 'data'
    elif 'mocks' in filename:
        home_dir = os.path.join(os.environ['THESIS_DIR'], 'eBOSS-QSO-PNG')
        spectra_dir = os.path.join(home_dir, 'measurements', 'spectra', 'mocks')
        return os.path.relpath(filename, spectra_dir).split(os.path.sep)[0]
    else:
        raise ValueError("do not understand file name '%s'" %filename)


class QSOFitPreparer(object):
    """
    Class to prepare a QSO power spectrum fit.
    """
    def __init__(self, spectra_file, stats, kmin=0.0001, kmax=0.4, overwrite=False, quiet=False):

        self.spectra_file = os.path.abspath(spectra_file)
        self.stats = stats
        self.kmin = kmin
        self.kmax = kmax
        self.overwrite = overwrite
        self.quiet = quiet

        # dictionary of info from the filename
        kind = get_spectra_type(self.spectra_file)
        info = info_from_filename(kind, os.path.basename(self.spectra_file))

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
        self.config = eBOSSConfig(self.sample, self.version)

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

        ns = parser.parse_args()
        return cls(**vars(ns))

    def write_data(self):
        """
        Write the necessary data file.
        """
        # the output data file
        stats = '+'.join(self.stats)
        filename = f"poles_{self.version}-QSO-{self.sample}_{stats}_{self.hashstr}.dat"
        output = os.path.join(self.config.fits_data_dir, filename)

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

        if not os.path.exists(output) or self.overwrite:

            # the name of the (unformatted) window file
            window_file = self._find_window_measurement()

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
        filename = f"poles_{self.version}-QSO-{self.sample}_{stats}_{self.hashstr}.dat"
        output = os.path.join(self.config.fits_covariance_dir, filename)

        if not os.path.exists(output) or self.overwrite:

            P0_FKP = self.hashinput.get('P0_FKP', None)
            if P0_FKP is None: P0_FKP = 0.
            b1 = bias_model(self.z_eff)
            sigma_fog = 4.0
            compute_covariance(self.config, self.stats, self.z_eff, b1, sigma_fog,
                                self.hashinput['zmin'], self.hashinput['zmax'],
                                P0_FKP, output, kmin=self.kmin, kmax=self.kmax, dk=0.005,
                                quiet=self.quiet)
        else:
            if not self.quiet:
                print('skipping covariance preparation...')
        self.covariance_file = output

    def _find_window_measurement(self):
        """
        Try to find and return a matching window function file
        """
        from glob import glob

        # the directory holding any window results
        home_dir = os.path.join(os.environ['THESIS_DIR'], 'eBOSS-QSO-PNG')
        dirname = os.path.join(home_dir, 'measurements', 'window', self.version)

        filename = f"RR_eboss_{self.version}-QSO-{self.sample}-*.json"
        pattern = os.path.join(dirname, filename)

        # search all file matches
        for f in glob(pattern):
            hashinput = get_hashkeys(f, 'SurveyDataPairCount')

            # compare zmin and zmax
            x = [self.hashinput[k] for k in ['zmin', 'zmax']]
            y = [hashinput[k] for k in ['zmin', 'zmax']]
            if numpy.allclose(x, y):
                return f

        raise ValueError(f"no window file match found for pattern '{pattern}'")
