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
            new_content = content.replace(
                "/global/cscratch1/sd/nhand/eBOSS", "$EBOSS_DIR")
            with open(f, 'w') as ff:
                ff.write(new_content)


def load_ezmock_driver(box_num, version, sample, krange, params, z_weighted, p=None):
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

    d = os.path.join(os.environ['EBOSS_DIR'], 'fits',
                     'results', 'mocks', 'ezmock', version)
    d = os.path.join(d, krange, params, '0.8-2.2')
    assert os.path.isdir(d), "'%s' directory not found" % d

    if p is None or p == 1.6:
        p = [None, 1.6]
    else:
        p = [p]

    matches = glob(os.path.join(d, f'QSO-{sample}-{box_num:04d}-*'))
    match = None
    for f in matches:
        hashkeys = get_hashkeys(f, None)
        if hashkeys['p'] in p and hashkeys['z-weighted'] == z_weighted:
            match = f

    if match is None:
        raise ValueError((f"no matches found: version={version}, sample={sample}, "
                          f"box={box_num}, krange={krange}, params={params}, z_weighted={z_weighted}, p={p}"))

    # load the driver
    driver = FittingDriver.from_directory(match)
    r = sorted(glob(os.path.join(match, '*.npz')),
               key=os.path.getmtime, reverse=True)
    if len(r) == 0:
        raise ValueError(
            "warning: no npz results found in directory '%s'" % os.path.normpath(match))
    driver.results = r[0]

    return driver


def load_ezmock_results(version, sample, krange, params, z_weighted, p=None):
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

    d = os.path.join(os.environ['EBOSS_DIR'], 'fits',
                     'results', 'mocks', 'ezmock', version)
    d = os.path.join(d, krange, params, '0.8-2.2')
    assert os.path.isdir(d), "'%s' directory not found" % d

    if p is None or p == 1.6:
        p = [None, 1.6]
    else:
        p = [p]

    matches = glob(os.path.join(d, f'QSO-{sample}-0001-*'))
    match = None
    for f in matches:
        hashkeys = get_hashkeys(f, None)
        if hashkeys['p'] in p and hashkeys['z-weighted'] == z_weighted:
            match = f

    if match is None:
        raise ValueError((f"no matches found: version={version}, sample={sample}, "
                          f"krange={krange}, params={params}, z_weighted={z_weighted}, p={p}"))

    # load the driver
    driver = FittingDriver.from_directory(match)
    dirname, basename = os.path.split(match)
    pattern = os.path.join(dirname, basename.replace('0001', '*'))

    data = defaultdict(list)
    matches = glob(pattern)

    for f in matches:
        r = sorted(glob(os.path.join(f, '*.npz')),
                   key=os.path.getmtime, reverse=True)
        if len(r) == 0:
            raise ValueError(
                "warning: no npz results found in directory '%s'" % os.path.normpath(f))

        th = ParameterSet.from_file(
            os.path.join(f, 'params.dat'), tags='theory')
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


def load_joint_data_results(kmin, z_weighted, p, ells=[0, 2]):
    """
    Load a set of data joint NGC + SGC fit results.
    """
    from eboss_qso.measurements.utils import make_hash
    from pyRSD.rsdfit.results import EmceeResults

    # the data to load
    kws = {}
    kws['version'] = 'v1.9f'
    kws['krange'] = '%s-0.3' % kmin
    kws['params'] = 'basemodel-N-fnl'
    kws['zrange'] = '0.8-2.2'
    kws['z_weighted'] = z_weighted
    kws['p'] = p

    if ells == [0]:
        kws['ells'] = ells

    hashstr = make_hash(kws)

    d = os.path.join(os.environ['EBOSS_DIR'],
                     'fits', 'results', 'data', 'v1.9f')
    d = os.path.join(d, kws['krange'], kws['params'], kws['zrange'])
    assert os.path.isdir(d), "'%s' directory not found" % d

    matches = glob(os.path.join(d, f'QSO-N+S-{hashstr}'))
    assert len(matches) == 1
    match = matches[0]

    if match is None:
        raise ValueError("no matches found for joint NGC + SGC data fits!")

    r = sorted(glob(os.path.join(match, '*.npz')),
               key=os.path.getmtime, reverse=True)
    assert len(
        r) > 0, "no npz results found in directory '%s'" % os.path.normpath(f)

    return EmceeResults.from_npz(r[0])


def load_data_results(version, sample, krange, params, zrange, z_weighted, ells=[0, 2], p=None):
    """
    Load a set of data fit results.

    Returns a structued numpy array holding best-fit values for all free
    parameters all mocks.
    """
    from pyRSD.rsdfit import FittingDriver

    assert sample in ['N', 'S']
    assert version in ['v1.8', 'v1.9f']

    d = os.path.join(os.environ['EBOSS_DIR'],
                     'fits', 'results', 'data', version)
    d = os.path.join(d, krange, params, zrange)
    assert os.path.isdir(d), "'%s' directory not found" % d

    if p is None or p == 1.6:
        p = [None, 1.6]
    else:
        p = [p]

    matches = glob(os.path.join(d, f'QSO-{sample}-*'))
    match = None
    stats = ['P%d' % ell for ell in ells]
    for f in matches:
        hashkeys = get_hashkeys(f, None)
        if hashkeys['p'] in p and hashkeys['z-weighted'] == z_weighted and hashkeys['stats'] == stats:
            match = f

    if match is None:
        raise ValueError((f"no matches found: version={version}, sample={sample}, "
                          f"krange={krange}, params={params}, zrange={zrange}, "
                          f"z_weighted={z_weighted}, p={p}"))

    # load the driver
    driver = FittingDriver.from_directory(match)

    r = sorted(glob(os.path.join(match, '*.npz')),
               key=os.path.getmtime, reverse=True)
    assert len(
        r) > 0, "no npz results found in directory '%s'" % os.path.normpath(f)
    driver.results = r[0]
    return driver
