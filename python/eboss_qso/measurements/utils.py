import hashlib
import json
from nbodykit.utils import JSONEncoder
import os

def compute_effective_redshift(cat):
    """
    Compute the effective redshift of a CatalogSource.
    """
    # the total weight
    total_weight =  cat['Weight']*cat['FKPWeight']

    # effective redshift
    zeff = (total_weight*cat['Z']).sum() / total_weight.sum()

    return cat.compute(zeff)

def compute_effective_nbar(cat):
    """
    Compute the effective number density of a CatalogSource.
    """
    # the total weight
    total_weight =  cat['Weight']*cat['FKPWeight']

    # effective nbar
    nbar = (total_weight*cat['NZ']).sum() / total_weight.sum()

    return cat.compute(nbar)


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

    # filename is a directory --> FIT result
    if os.path.isdir(filename):
        filename = os.path.join(os.path.abspath(filename), 'hashinfo.json')
        if not os.path.exists(filename):
            return
        import json

        # use json to load
        d = {}
        with open(filename, 'r') as ff:
            d.update(json.load(ff))

        # echo hash info for the spectra file too
        spectra_file = d.pop('spectra_file')
        d.update(get_hashkeys(spectra_file, 'ConvolvedFFTPower'))
    else:
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
    parser.add_argument('filenames', type=str, nargs='+', help=h)

    h = 'the result class'
    parser.add_argument('--cls', type=str, default='ConvolvedFFTPower', help=h)

    h = 'only show entries with this min z value'
    parser.add_argument('--zmin', type=float, help=h)

    h = 'only show entries with this max z value'
    parser.add_argument('--zmax', type=float, help=h)

    h = 'only show entries with this z-weighted values'
    parser.add_argument('--z-weighted', choices=[0,1], type=int, help=h)

    ns, unknown = parser.parse_known_args()


    for filename in ns.filenames:
        d = get_hashkeys(filename, ns.cls)

        # filter
        if ns.zmin is not None and d.get('zmin', None) != ns.zmin:
            continue
        if ns.zmax is not None and d.get('zmax', None) != ns.zmax:
            continue
        if ns.z_weighted is not None and d.get('z-weighted', None) != bool(ns.z_weighted):
            continue

        # print
        print(f"{filename}" + '\n' + '-'*40)
        for k in sorted(d.keys()):
            print("%-10s = %s" %(k, str(d[k])))
