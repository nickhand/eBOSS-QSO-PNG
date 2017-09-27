import hashlib
import json
from nbodykit.utils import JSONEncoder

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
    return hashlib.sha1(json.dumps(d, sort_keys=True, cls=JSONEncoder).encode()).hexdigest()[:N]
