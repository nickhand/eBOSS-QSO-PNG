import hashlib
import json
from nbodykit.utils import JSONEncoder

def redshift_range_type(s):
    """
    Allow input redshift ranges via the command line
    """
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

def get_hashkeys(filename, cls):
    """
    Return a dict of key/values that generated a filename with a unique hash ID

    Parameters
    ----------
    filename : str
        the name of file to load
    cls : str
        the result class
    """
    from nbodykit import lab

    # get the result class
    cls = getattr(lab, cls)
    r = cls.load(filename)

    # need hashkeys
    assert 'hashkeys' in r.attrs, "result filename does not have 'hashkeys' attribute"

    # the dict
    d = {k:r.attrs[k] for k in r.attrs['hashkeys']}
    return d

def echo_hash():
    """
    Echo the key/values that generated the hash in the input filename
    """
    import argparse
    desc = 'echo the key/values that generated the hash in the input filename'
    parser = argparse.ArgumentParser(description=desc)

    h = 'the input file name'
    parser.add_argument('filename', type=str, help=h)

    h = 'the result class'
    parser.add_argument('--cls', type=str, default='ConvolvedFFTPower', help=h)

    ns = parser.parse_args()

    d = get_hashkeys(ns.filename, ns.cls)
    for k in sorted(d.keys()):
        print("%-10s = %s" %(k, str(d[k])))
