import os
from .utils import get_hashkeys
from glob import glob
import numpy

def load_data_spectra(version, sample, p=None, focal_weights=True):
    """
    Load a data measurement result.
    """
    from nbodykit.lab import ConvolvedFFTPower

    assert sample in ['N', 'S']
    assert version in ['v1.8', 'v1.9f']

    # the directory
    eboss_dir = os.environ['EBOSS_DIR']
    d = os.path.join(eboss_dir, 'measurements', 'spectra', 'data', version)

    filename = f'poles_eboss_{version}'
    if focal_weights: filename += '-focal-*' + sample + '-*'
    matches = glob(os.path.join(d, filename))
    assert len(matches) > 0

    for f in matches:
        hashkeys = get_hashkeys(f, 'ConvolvedFFTPower')
        if hashkeys['p'] == p:
            r = ConvolvedFFTPower.load(f)
            r.poles['power_0'].real -= r.attrs['shotnoise']
            return r

    raise ValueError("no matches found!")

def load_ezmock_spectra(version, sample, p=None, box=None, average=True):
    """
    Load a ezmock measurement result.
    """
    from nbodykit.lab import ConvolvedFFTPower

    assert sample in ['N', 'S']
    assert version in ['v1.8e-fph', 'v1.8e-no', 'v1.8e-reg']

    # the directory
    eboss_dir = os.environ['EBOSS_DIR']
    d = os.path.join(eboss_dir, 'measurements', 'spectra', 'mocks', 'ezmock', version)

    if box is not None:
        filename = f'poles_zevoEZmock_{version}_QSO-{sample}_{box:04d}-*.json'
    else:
        filename = f'poles_zevoEZmock_{version}_QSO-{sample}_0001-*.json'

    matches = glob(os.path.join(d, filename))
    assert len(matches) > 0

    hashstr = None
    for f in matches:
        hashkeys = get_hashkeys(f, 'ConvolvedFFTPower')
        if hashkeys['p'] == p:
            hashstr = os.path.splitext(f)[0][-10:]
            break

    if hashstr is None:
        raise ValueError("no matches found!")

    if box is not None:
        f = os.path.join(d, f'poles_zevoEZmock_{version}_QSO-{sample}_{box:04d}-{hashstr}.json')
        r = ConvolvedFFTPower.load(f)
        r.poles['power_0'].real -= r.attrs['shotnoise']
        return r
    else:
        files = glob(os.path.join(d, f'poles_zevoEZmock_{version}_QSO-{sample}_*-{hashstr}.json'))
        results = [ConvolvedFFTPower.load(f) for f in files]

        for r in results:
            r.poles['power_0'].real -= r.attrs['shotnoise']

        if not average:
            data = numpy.asarray([r.poles.data for r in results])
            return data
        else:
            data = results[0].poles.copy()
            for k in results[0].poles.variables:
                data[k] = numpy.asarray([r.poles[k] for r in results]).mean(axis=0)
            return data

    raise ValueError("no matches found!")

def info_from_filename(kind, filename):
    """
    Return a dictionary of information from a file name.
    """
    import re

    assert kind in ['data', 'ezmock']

    if kind == 'data':
        pattern = r'poles\_eboss\_(?P<version>[\w\.]+)\-(?P<focal>[a-z]*)\-*QSO\-(?P<sample>[SN])\-(?P<hashstr>[\w+]{10})\.json'
    elif kind == 'ezmock':
        pattern = r'poles\_zevoEZmock\_(?P<version>[\w\.\-]+)_QSO\-(?P<sample>[SN])\_(?P<box>[0-9]{4})\-(?P<hashstr>[\w+]{10})\.json'
    else:
        raise ValueError("invalid kind '%s'" %kind)
    return re.match(pattern, filename).groupdict()

def save_RR_paircount(r, sample, version, **kwargs):
    """
    Save the randoms paircount result with a unique hash string based on the
    meta-data of the result.

    Parameters
    ----------
    r :
        the nbodykit result object; should have a save() method
    sample : str
        the eBOSS QSO sample
    version : str
        the version number
    **kwargs :
        additional meta-data to save
    """
    from .utils import make_hash
    from . import results_dir

    # output path
    output_dir =  os.path.join(results_dir, 'window', version)
    filename = f"RR_eboss_{version}-QSO-{sample}"

    # save the extra meta-data
    for k in kwargs:
        r.attrs[k] = kwargs[k]

    # make the hash
    usekeys = ['redges_str', 'zmin', 'zmax', 'N', 'subsample']

    r.attrs['hashkeys'] = usekeys
    id_str = make_hash(r.attrs, usekeys=usekeys)
    r.save(os.path.join(output_dir, filename + '-' + id_str + '.json'))


def save_data_spectra(r, sample, version, focal_weights, **kwargs):
    """
    Save the power spectrum result with a unique hash string based on the
    meta-data of the result.

    Parameters
    ----------
    r :
        the nbodykit result object; should have a save() method
    sample : str
        the eBOSS QSO sample
    version : str
        the version number
    focal_weights : bool
        whether the focal plane weights were used
    **kwargs :
        additional meta-data to save
    """
    from .utils import make_hash
    from . import results_dir

    # output path
    output_dir =  os.path.join(results_dir, 'spectra', 'data', version)
    if focal_weights: version += '-focal'
    filename = f"poles_eboss_{version}-QSO-{sample}"

    # save the extra meta-data
    for k in kwargs:
        r.attrs[k] = kwargs[k]
    r.attrs['z-weighted'] = r.attrs.get('p', None) is not None

    # make the hash
    usekeys = ['Nmesh', 'P0_FKP', 'dk', 'kmin', 'mesh.window', 'mesh.interlaced',
                'poles', 'p', 'BoxPad', 'zmin', 'zmax', 'z-weighted']

    r.attrs['hashkeys'] = usekeys
    id_str = make_hash(r.attrs, usekeys=usekeys)
    r.save(os.path.join(output_dir, filename + '-' + id_str + '.json'))

def save_ezmock_spectra(r, box, sample, version, subversion, **kwargs):
    """
    Save the EZ mock power spectrum result with a unique hash string based on
    the meta-data of the result.

    Parameters
    ----------
    r :
        the nbodykit result object; should have a save() method
    box : int
        the box number
    sample : str
        the eBOSS QSO sample
    version : str
        the version number
    subversion : str
        the sub version
    **kwargs :
        additional meta-data to save
    """
    from .utils import make_hash
    from . import results_dir

    version = '-'.join([version, subversion])

    # output path
    output_dir =  os.path.join(results_dir, 'spectra', 'mocks', 'ezmock', version)
    filename = f"poles_zevoEZmock_{version}_QSO-{sample}_{box:04d}"

    # save the extra meta-data
    for k in kwargs:
        r.attrs[k] = kwargs[k]
    r.attrs['z-weighted'] = r.attrs.get('p', None) is not None

    # make the hash
    usekeys = ['Nmesh', 'P0_FKP', 'dk', 'kmin', 'mesh.window', 'mesh.interlaced',
                'poles', 'p', 'BoxPad', 'zmin', 'zmax', 'z-weighted']

    r.attrs['hashkeys'] = usekeys
    id_str = make_hash(r.attrs, usekeys=usekeys)
    r.save(os.path.join(output_dir, filename + '-' + id_str + '.json'))
