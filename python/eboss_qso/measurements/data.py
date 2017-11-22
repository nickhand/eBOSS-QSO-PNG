
from . import data_dir, DATA_VERSIONS
import os
from nbodykit.lab import FITSCatalog

def read_data(sample, version, focal_weights=False):
    """
    Read a eBOSS QSO data file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    focal_weights : bool, optional
        whether we are using focal plane redshift error corrections
    """
    assert version in DATA_VERSIONS

    # get the file path
    filename = f'eboss_{version}-QSO-{sample}-eboss_{version}'
    if focal_weights: filename += '-focal'
    filename += '.dat.fits'
    path = os.path.join(data_dir, 'data', version, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    if focal_weights: usecols += ['WEIGHT_FOCAL']
    return FITSCatalog(path)[usecols]


def read_randoms(sample, version):
    """
    Read a eBOSS QSO randoms file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    """
    assert version in DATA_VERSIONS

    # get the file path
    filename = f'eboss_{version}-QSO-{sample}-eboss_{version}.ran.fits'
    path = os.path.join(data_dir, 'data', version, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    return FITSCatalog(path)[usecols]

def finalize_data(s, cosmo, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', 'Weight', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    cosmo : Cosmology
        the cosmology parameters
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """
    from nbodykit.transform import SkyToCartesian

    # add the Position column
    s['Position'] = SkyToCartesian(s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)

    # add systematic weights
    if 'WEIGHT_NOZ' in s and 'WEIGHT_FOCAL' not in s:
        s['Weight'] = s['WEIGHT_SYSTOT'] * (s['WEIGHT_NOZ'] + s['WEIGHT_CP'] - 1.)
    elif 'WEIGHT_FOCAL' in s:
        s['Weight'] = s['WEIGHT_SYSTOT'] * s['WEIGHT_CP'] * s['WEIGHT_FOCAL']

    # FKP WEIGHT
    if P0_FKP is not None:
        s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
