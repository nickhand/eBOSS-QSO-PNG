import os
from glob import glob
import numpy
from ..measurements.utils import get_hashkeys


def load_ezmock_results(version, sample, krange, params, p=None):
    """
    Load a set of ezmock fit results.

    Returns a structued numpy array holding best-fit values for all free
    parameters all mocks.
    """
    from pyRSD.rsdfit.results import LBFGSResults
    from collections import defaultdict

    assert sample in ['N', 'S']
    assert version in ['v1.8e-no', 'v1.8e-fph', 'v1.8e-reg']

    d = os.path.join(os.environ['EBOSS_DIR'], 'fits', 'results', 'mocks', 'ezmock', version)
    d = os.path.join(d, krange, params, '0.8-2.2')
    assert os.path.isdir(d), "'%s' directory not found" % d

    matches = glob(os.path.join(d, f'QSO-{sample}-0001-*'))
    match = None
    for f in matches:
        hashkeys = get_hashkeys(f, None)
        if hashkeys['p'] == p:
            match = f

    assert match is not None, "no matches found!"
    pattern = match.replace('0001', '*')

    data = defaultdict(list)
    matches = glob(pattern)
    for f in matches:
        r = sorted(glob(os.path.join(f, '*.npz')), key=os.path.getmtime, reverse=True)
        assert len(r) > 0, "no npz results found in directory '%s'" %os.path.normpath(f)
        r = LBFGSResults.from_npz(r)
        for param in r.free_names:
            data[param].append(r[param])

    params = list(data.keys())
    out = numpy.empty(len(matches), dtype=list(zip(params, ['f8']*len(params))))
    for param in params:
        out[param] = numpy.array(data[param])

    return out
