import os
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import itertools

def parametrize(params):
    """
    Execute a function for each of the input parameters.
    """
    keys = list(params.keys())
    params = list(itertools.product(*[params[k] for k in params]))
    def wrapped(func):
        def func_wrapper(*args, **kwargs):
            for p in params:
                kwargs.update(dict(zip(keys, p)))
                func(*args, **kwargs)

        return func_wrapper
    return wrapped

def get_mock_version(kind, version0):
    if kind == 'ezmock':
        version = version0[:4]
    elif kind == 'data':
        version = version0
    else:
        raise ValueError("do not understand 'kind' = '%s'" %kind)

    return version

class eBOSSConfig(object):
    """
    Configuration class for the eBOSS results + fits

    Parameters
    ----------
    sample : 'N', 'S'
        the eBOSS sample we are fitting
    version : str
        the version we are fitting
    """
    def __init__(self, sample, version, kind):

        assert sample in ['N', 'S']
        self.version = version
        self.sample = sample
        self.kind = kind
        self._kind_tag = kind if kind == 'data' else 'mocks/' + kind

    def _get_nbar_file(self):
        version = get_mock_version(self.kind, self.version)
        filename = f'nbar-eboss_{version}-QSO-{self.sample}-eboss_{version}.dat'
        return os.path.join(self.data_dir, 'meta', filename)

    def _get_effective_area(self):
        """
        Return the effective area in squared degrees from the relevant nbar file.
        """
        f = self._get_nbar_file()
        lines = open(f, 'r').readlines()
        return float(lines[1].split()[0])

    def _get_nbar(self):
        """
        Return a spline fit to n(z).
        """
        f = self._get_nbar_file()
        nbar = numpy.loadtxt(f, skiprows=3)
        return spline(nbar[:,0], nbar[:,3])

    @property
    def home_dir(self):
        return os.path.join(os.environ['THESIS_DIR'], 'eBOSS-QSO-PNG')

    @property
    def fits_input_dir(self):
        return os.path.join(os.environ['THESIS_DIR'], 'eBOSS-QSO-PNG', 'fits', 'input', self._kind_tag)

    @property
    def fits_covariance_dir(self):
        return os.path.join(self.fits_input_dir, 'covariance', self.version)

    @property
    def fits_window_dir(self):
        return os.path.join(self.fits_input_dir, 'window', self.version)

    @property
    def fits_data_dir(self):
        return os.path.join(self.fits_input_dir, 'spectra', self.version)

    @property
    def fits_params_dir(self):
        return os.path.join(self.fits_input_dir, 'params')

    @property
    def fits_results_dir(self):
        return os.path.join(os.environ['THESIS_DIR'], 'eBOSS-QSO-PNG', 'fits', 'results', self._kind_tag, self.version)

    @property
    def spectra_dir(self):
        return os.path.join(self.home_dir, 'measurements', 'spectra')

    @property
    def data_dir(self):
        version = get_mock_version(self.kind, self.version)
        return os.path.join(self.home_dir, 'data', version)

    @property
    def fsky(self):
        try:
            return self._fsky
        except AttributeError:
            eff_area = self._get_effective_area()
            self._fsky = eff_area * (numpy.pi/180.)**2 / (4*numpy.pi)
            return self._fsky

    @property
    def nbar(self):
        try:
            return self._nbar
        except AttributeError:
            self._nbar = self._get_nbar()
            return self._nbar

    @property
    def cosmo(self):
        from pyRSD.rsd.cosmology import Cosmology
        if self.kind == 'data':
            params = {'H0': 67.6, 'Neff': 3.046, 'Ob0': 0.04814257203879415,'Om0': 0.31,'Tcmb0': 2.7255,'m_nu': 0.0,'n_s': 0.97,'sigma8': 0.80}
        elif self.kind == 'ezmock':
            params = {'H0': 67.77, 'Neff': 3.046, 'Ob0': 0.048206, 'Om0': 0.307115, 'Tcmb0': 2.7255, 'm_nu': 0.0, 'n_s': 0.9611,'sigma8': 0.8225}
        else:
            raise ValueError("cannot compute cosmology in eBOSSConfig")

        return Cosmology(flat=True, **params)


from .runner import RSDFitRunner
