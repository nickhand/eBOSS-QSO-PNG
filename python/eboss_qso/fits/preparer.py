import os
import numpy as np
from nbodykit.lab import ConvolvedFFTPower
import tempfile

from ..measurements.utils import find_window_measurement, make_hash
from ..measurements import get_hashkeys, load_data_spectra, load_ezmock_spectra
from ..measurements.weights import bias_model

from . import eBOSSConfig
from .window import compute_window
from .covariance import compute_analytic_covariance, compute_ezmock_covariance


class QSOFitPreparer(object):
    """
    Class to prepare a QSO power spectrum fit.
    """

    def __init__(self, kind,
                 version,
                 sample,
                 ells,
                 z_eff,
                 p=None,
                 kmin=0.0001,
                 kmax=0.4,
                 cov_type='analytic',
                 overwrite=False,
                 error_rescale=1.0,
                 quiet=False,
                 use_temp_files=False,
                 box=None):

        assert kind in ['data', 'ezmock']
        assert sample in ['N', 'S']
        assert len(z_eff) == len(ells)
        if kind == 'ezmock':
            assert box is not None

        self.kind = kind
        self.version = version
        self.sample = sample
        self.ells = ells
        self.z_eff = z_eff
        self.box = box
        self.p = p

        self.kmin = kmin
        self.kmax = kmax
        self.overwrite = overwrite
        self.quiet = quiet
        self.error_rescale = error_rescale
        self.use_temp_files = use_temp_files
        self.cov_type = cov_type

        # make the hash string
        cols = ['kind', 'version', 'sample', 'ells', 'z_eff', 'p', 'kmin', 'kmax']
        attrs = {k: getattr(self, k) for k in cols}
        self.hashstr = make_hash(attrs)

        # which function do we need to load measurements
        if self.kind == 'data':
            self.loader = load_data_spectra
        else:
            self.loader = load_ezmock_spectra

        # get the loader kwargs
        self.loader_kws = {k: getattr(self, k)
                           for k in ['version', 'sample', 'p']}
        if self.box is not None:
            self.loader_kws['box'] = self.box

        # data attrs from the spectra file
        self.attrs = self.loader(ell=self.ells[0], **self.loader_kws).attrs

        # the configuration
        self.config = eBOSSConfig(self.sample, self.version, self.kind)

        # determine the statistic names
        tag = 'ngc' if self.config.sample == 'N' else 'sgc'
        self.stat_names = ['pole_' + tag + '_%d' % ell for ell in ells]
        self.stats = ['P%d' % ell for ell in ells]

        # and run
        self.write_data()
        self.write_window()
        self.write_covariance()

    def write_data(self):
        """
        Write the necessary data file.
        """
        # the output data file
        stats = '+'.join(self.stats)
        filename = f"poles_{self.version}-QSO-{self.sample}"
        if self.box is not None:
            filename += '-%04d' % self.box
        filename += f"_{stats}_{self.hashstr}.dat"
        output = os.path.join(self.config.fits_data_dir, filename)

        if self.use_temp_files:
            self._data_file = output
            output = tempfile.mktemp()

        # make the data file
        if not os.path.exists(output) or self.overwrite:
            self._write_data(output)
        else:
            if not self.quiet:
                print('skipping data preparation...')
        self.data_file = output

    def _write_data(self, output):
        """
        Internal function to write out the multipole data to 
        a plaintext file.

        Parameters
        ----------
        output : str
            the output file name
        """
        from pyRSD.rsdfit.data import PowerMeasurements

        power = []
        for ell in self.ells:
            # get the ConvolvedFFTPower object
            r = self.loader(ell=ell, subtract_shot_noise=True,
                            **self.loader_kws)

            # trim to the correct k range
            valid = (r.poles['k'] >= self.kmin) & (r.poles['k'] <= self.kmax)
            poles = r.poles.data[valid]

            # add the power array
            power.append(poles['power_%s' % ell].real)

        # convert to array
        power = np.vstack(power).T

        # make the data array
        dtype = [('k', 'f8'), ('power', 'f8')]
        data = np.zeros_like(power, dtype=dtype)

        # store k and P(k)
        data['power'] = power[:]
        data['k'] = np.repeat(poles['k'][:, None], data.shape[1], axis=1)

        # initialize the power measurements object
        measurements = PowerMeasurements.from_array(self.stat_names, data)

        # write to file
        if not self.quiet:
            print(f'saving {output}...')
        measurements.to_plaintext(output)

    def write_window(self):
        """
        Write the necessary window file.
        """
        # get zmin/zmax
        zmin = self.attrs['zmin']
        zmax = self.attrs['zmax']

        self.window_file = []
        self._window_file = []

        # the window file
        for ell in self.ells:

            meta = {'zmin': zmin, 'zmax': zmax, 'p': self.p, 'ell': ell}
            hashstr = make_hash(meta)
            filename = f"poles_{self.version}-QSO-{self.sample}_{hashstr}.dat"
            output = os.path.join(self.config.fits_window_dir, filename)

            if self.use_temp_files:
                self._window_file.append(output)
                output = tempfile.mktemp()

            # make the window
            if not os.path.exists(output) or self.overwrite:

                # the name of the (unformatted) window file
                if self.kind == 'ezmock':
                    version = self.version[:4]
                elif self.kind == 'data':
                    version = self.version
                else:
                    raise ValueError(
                        "do not understand 'kind' = '%s'" % self.kind)
                window_file = find_window_measurement(version, self.sample,
                                                      zmin, zmax, self.p, ell)

                print('using window file %s...' % window_file)
                ells = [0, 2, 4, 6, 8, 10]
                compute_window(window_file, ells, output, smin=1e-2,
                               smax=1e4, quiet=self.quiet)
            else:
                if not self.quiet:
                    print('skipping window preparation...')
            self.window_file.append(output)

    def write_covariance(self):
        """
        Write the necessary covariance file.
        """
        # the output data file
        stats = '+'.join(self.stats)
        filename = f"poles_{self.version}-QSO-{self.sample}-{self.cov_type}"
        if self.cov_type == 'analytic':
            if self.box is not None and self.box == 'mean':
                filename += '-mean'
        filename += f"_{stats}_{self.hashstr}.dat"
        output = os.path.join(self.config.fits_covariance_dir, filename)

        if self.use_temp_files:
            self._covariance_file = output
            output = tempfile.mktemp()

        # make the covariance
        if not os.path.exists(output) or self.overwrite:

            if self.cov_type == 'analytic':
                raise ValueError("analytic covariance not currently supported")
                # from pyRSD.rsd import QuasarSpectrum

                # P0_FKP = self.attrs.get('P0_FKP', None)
                # if P0_FKP is None:
                #     P0_FKP = 0.

                # # load the model
                # model = QuasarSpectrum(z=self.z_eff, params=config.cosmo)

                # # compute
                # C = compute_analytic_covariance(self.config, self.stats, model,
                #                                 self.attrs['zmin'], self.attrs['zmax'],
                #                                 P0_FKP, kmin=self.kmin, kmax=self.kmax, dk=0.005)

                # C *= self.error_rescale

                # # and save
                # if not self.quiet:
                #     print("saving %s..." % output)
                # C.attrs.clear()
                # C.to_plaintext(output)

            else:
                # compute the covariance from the mocks
                C = compute_ezmock_covariance('v1.8e-fph',
                                              self.sample,
                                              self.ells,
                                              kmin=self.kmin,
                                              kmax=self.kmax,
                                              p=self.p)

                # rescale
                C *= self.error_rescale

                # and save
                C.to_plaintext(output)

        else:
            if not self.quiet:
                print('skipping covariance preparation...')
        self.covariance_file = output
