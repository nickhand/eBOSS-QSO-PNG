import os
from glob import glob
import numpy
from ..measurements.utils import get_hashkeys

def fix_result_paths():
    """
    Fix paths local to NERSC in fitting results.
    """
    from argparse import ArgumentParser
    import os

    parser = ArgumentParser(description='fix EBOSS_DIR path')
    parser.add_argument('dirname', type=str)
    ns = parser.parse_args()
    assert os.path.isdir(ns.dirname)

    for dirpath, dirnames, filenames in os.walk(ns.dirname):
        if 'params.dat' in filenames:
            f = os.path.join(dirpath, 'params.dat')
            content = open(f, 'r').read()
            new_content = content.replace("/global/cscratch1/sd/nhand/eBOSS", "$EBOSS_DIR")
            with open(f, 'w') as ff:
                ff.write(new_content)


def load_ezmock_results(version, sample, krange, params, p=None):
    """
    Load a set of ezmock fit results.

    Returns a structued numpy array holding best-fit values for all free
    parameters all mocks.
    """
    from pyRSD.rsdfit.results import LBFGSResults
    from pyRSD.rsdfit import FittingDriver
    from pyRSD.rsdfit.parameters import ParameterSet
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

    # load the driver
    driver = FittingDriver.from_directory(match)

    assert match is not None, "no matches found!"
    pattern = match.replace('0001', '*')

    data = defaultdict(list)
    matches = glob(pattern)
    for f in matches:
        r = sorted(glob(os.path.join(f, '*.npz')), key=os.path.getmtime, reverse=True)
        assert len(r) > 0, "no npz results found in directory '%s'" %os.path.normpath(f)

        th = ParameterSet.from_file(os.path.join(f, 'params.dat'), tags='theory')
        r = LBFGSResults.from_npz(r[0])
        for param in r.free_names:
            data[param].append(r[param])
            th[param].value = r[param]

        # add fsigma8
        if 'f' in r.free_names and 'sigma8_z' in r.free_names:
            data['fsigma8'].append(r['f'] * r['sigma8_z'])

        # the prior to add back
        lnprior = sum(par.lnprior for par in th.free)

        # add the reduced chi2
        red_chi2 = (2*(r.min_chi2 + lnprior)) / driver.dof
        data['red_chi2'].append(red_chi2)

    params = list(data.keys())
    dtype = list(zip(params, ['f8']*len(params)))
    out = numpy.empty(len(matches), dtype=dtype)
    for param in out.dtype.names:
        out[param] = numpy.array(data[param])

    return out
