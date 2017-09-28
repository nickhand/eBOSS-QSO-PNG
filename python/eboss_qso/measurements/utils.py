import hashlib
import json
from nbodykit.utils import JSONEncoder

def redshift_range_type(s):
    try:
        return tuple(map(float, s.split(',')))
    except:
        raise TypeError("redshift range must be zmin,zmax")

def trim_redshift_range(s, zmin=None, zmax=None):
    """
    Trim the redshift range of a CatalogSource.

    Parameters
    ----------
    zmin : float, optional
        minimum redshift to include (exclusive)
    zmax : float, optional
        maximum redshift to include (exclusive)
    """
    # trim the redshift range
    if zmin is None: zmin = 0.
    if zmax is None: zmax = 10.0

    return s[(s['Z'] > zmin)&(s['Z'] < zmax)]

def finalize_source(s, cosmo, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource by adding 'Position', 'Weight',
    and 'FKPWeight'.
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

def make_hash(attrs, usekeys=None, N=10):
    """
    Return a unique hash string for the subset of ``attrs`` specified
    by ``usekeys``.

    Parameters
    ----------
    attrs : dict
        the dictionary of meta-data to use to make the hashlib
    usekeys : list, optional
        only include these keys in the calculation
    N : int, optional
        return the first ``N`` characters from the hash string
    """
    if usekeys is None:
        d = attrs
    else:
        d = {k:attrs[k] for k in usekeys}

    s = json.dumps(d, sort_keys=True, cls=JSONEncoder).encode()
    return hashlib.sha1(s).hexdigest()[:N]
