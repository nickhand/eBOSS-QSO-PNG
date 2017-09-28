import os
from nbodykit.cosmology import Cosmology

# data and result directories
data_dir = '/global/cscratch1/sd/nhand/eBOSS'
results_dir = '/global/cscratch1/sd/nhand/Research/eBOSS/measurements'

# allowed data versions
DATA_VERSIONS = ['v1.6', 'v1.8', 'v1.9f']

# the fiducial DR12 cosmology
_h = 0.676
_Ob0 = 0.022/_h**2
_Ocdm0 = 0.31 - _Ob0
fidcosmo = Cosmology(h=_h, Omega_b=_Ob0, Omega_cdm=_Ocdm0, m_ncdm=None, T_cmb=2.7255)

from .data import read_data, read_randoms, save_data_spectra
from .utils import finalize_source, trim_redshift_range, redshift_range_type
from .zweights import fnl_weight, bias_weight
