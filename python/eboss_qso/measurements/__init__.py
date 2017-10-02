import os
from nbodykit.cosmology import Cosmology

# data and result directories
data_dir = '/global/cscratch1/sd/nhand/eBOSS'
results_dir = '/global/cscratch1/sd/nhand/Research/eBOSS/measurements'

# allowed data versions
DATA_VERSIONS = ['v1.6', 'v1.8', 'v1.9f']

# EZmock subversions
EZMOCK_VERSIONS = ['v1.8e']
EZMOCK_SUBVERSIONS = ['reg', 'no', 'fph']

# the fiducial DR12 cosmology
_h = 0.676
_Ob0 = 0.022/_h**2
_Ocdm0 = 0.31 - _Ob0
fidcosmo = Cosmology(h=_h, Omega_b=_Ob0, Omega_cdm=_Ocdm0, m_ncdm=None, T_cmb=2.7255)

# data functions
from .data import read_data, read_randoms, finalize_data

# ezmock function
from .ezmock import read_ezmock_data, read_ezmock_randoms, finalize_ezmock

# results
from .results import save_data_spectra, save_ezmock_spectra, save_RR_paircount

# utilities and weights
from .utils import trim_redshift_range, redshift_range_type, get_hashkeys,
                    compute_effective_redshift, compute_effective_nbar
from .zweights import fnl_weight, bias_weight
