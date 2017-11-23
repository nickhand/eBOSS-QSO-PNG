
from . import data_dir, MOCK_VERSIONS, MOCK_SUBVERSIONS
import os
from nbodykit.lab import CSVCatalog

def read_qpm_data(box, sample, version, subversion):
    """
    Read an eBOSS QSO QPM mock data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    subversion : str
        the mock sub version to load, i.e., one of 'reg', 'no', 'fph'
    """
    assert version in MOCK_VERSIONS
    assert subversion in MOCK_SUBVERSIONS

    sample = 'ngc' if sample == 'N' else 'sgc'

    # get the file path
    filename = f"qpmmock_QSO_{version}_veto_{sample}_{box:04d}_{subversion}.dat"
    path = os.path.join(data_dir, 'mocks', 'qpm', '-'.join([version,subversion]), filename)

    # load the source
    names = ['RA', 'DEC', 'Z', 'WEIGHT_FKP', 'COMP', 'NZ', 'WEIGHT_COMP']
    return CSVCatalog(path, names=names)


def read_qpm_randoms(sample, version):
    """
    Read a eBOSS QSO QPM mock randoms file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    """
    assert version in MOCK_VERSIONS

    sample = 'ngc' if sample == 'N' else 'sgc'

    # get the file path
    filename = f"qpm_random_v1.8.1_{sample}_vetoed_smeared.dat"
    path = os.path.join(data_dir, 'mocks', 'qpm', 'randoms', filename)

    # load the source
    names = ['RA', 'DEC', 'Z', 'WEIGHT_FKP', 'COMP', 'NZ', 'VETO']
    s = CSVCatalog(path, names=names)

    return s


def finalize_qpm(s, cosmo, P0_FKP=None):
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
