
from . import data_dir, EZMOCK_VERSIONS, EZMOCK_SUBVERSIONS
import os
from nbodykit.lab import CSVCatalog

def read_ezmock_data(box, sample, version, subversion):
    """
    Read an eBOSS QSO EZ mock data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    subversion : str
        the EZmock sub version to load, i.e., one of 'reg', 'no', 'fph'
    """
    assert version in EZMOCK_VERSIONS
    assert subversion in EZMOCK_SUBVERSIONS

    sample = 'ngc' if sample == 'N' else 'sgc'

    # get the file path
    filename = f'zevoEZmock_QSO_{version}_veto_{sample}_{box:04d}_{subversion}.dat'
    path = os.path.join(data_dir, 'mocks', '-'.join([version,subversion]), filename)

    # load the source
    names = ['RA', 'DEC', 'Z', 'WEIGHT_FKP', 'COMP', 'NZ', 'WEIGHT_COMP']
    return CSVCatalog(path, names=names)


def read_ezmock_randoms(sample, version):
    """
    Read a eBOSS QSO EZ mock randoms file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    """
    assert version in EZMOCK_VERSIONS

    sample = 'ngc' if sample == 'N' else 'sgc'

    # get the file path
    filename = f'randx100_QSO_{version}_veto_{sample}.dat'
    path = os.path.join(data_dir, 'mocks', 'randoms', filename)

    # load the source
    names = ['RA', 'DEC', 'Z', 'WEIGHT_FKP', 'COMP', 'NZ', 'VETO']
    s = CSVCatalog(path, names=names)

    return s


def finalize_ezmock(s, cosmo, P0_FKP=None):
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
    if 'WEIGHT_COMP' in s:
        s['Weight'] = s['WEIGHT_COMP']
    else:
        s['Weight'] = 1.0

    # FKP WEIGHT
    if P0_FKP is not None:
        s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
