
from . import data_dir, DATA_VERSIONS, fidcosmo
import os
from nbodykit.lab import FITSCatalog

def read_data(version, sample, focal_weights=False):
    """
    Read a eBOSS QSO data file.

    Parameters
    ----------
    version : str
        the string specifying which version to load
    sample : 'N' or 'S'
        the sample to load
    P0_FKP : float, optional
        the FKP P0 value to use
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
    return FITSCatalog(path, use_cache=True)[usecols]


def read_randoms(version, sample):
    """
    Read a eBOSS QSO randoms file.

    Parameters
    ----------
    version : str
        the string specifying which version to load
    sample : 'N' or 'S'
        the sample to load
    P0_FKP : float, optional
        the FKP P0 value to use
    focal_weights : bool, optional
        whether we are using focal plane redshift error corrections
    """
    assert version in DATA_VERSIONS

    # get the file path
    filename = f'eboss_{version}-QSO-{sample}-eboss_{version}.ran.fits'
    path = os.path.join(data_dir, 'data', version, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    return FITSCatalog(path, use_cache=True)[usecols]

def save_data_spectra(r, sample, version, **kwargs):
    """
    Save the power spectrum result with a unique hash string based on the
    meta-data of the result.
    """
    from .utils import make_hash
    from . import results_dir

    # output path
    output_dir =  os.path.join(results_dir, 'spectra', 'data', version)
    filename = f"poles_eboss_{version}-QSO-{sample}"

    # save the extra meta-data
    for k in kwargs:
        r.attrs[k] = kwargs[k]
    r.attrs['z-weighted'] = r.attrs.get('p', None) is not None

    # make the hash
    usekeys = ['Nmesh', 'P0_FKP', 'dk', 'kmin', 'mesh.window', 'mesh.interlaced',
                'poles', 'p', 'BoxPad', 'zmin', 'zmax', 'z-weighted']

    id_str = make_hash(r.attrs, usekeys=usekeys)
    r.save(filename + '-' + id_str + '.json')
