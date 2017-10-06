from ebossqso import config

import os
import numpy
from pyRSD.rsdfit.data import PowerMeasurements
from nbodykit.algorithms import ConvolvedFFTPower

def write_data(filename, stats, names, output, kmin=0., kmax=0.7, quiet=False):
    """
    Write out the multipole data to a plaintext file

    Parameters
    ----------
    filename : str
        the input spectra to load
    zmin : float
        the minimum redshift of the sample
    zmax : float
        the maximum redshift of the sample
    lmax : int
        the maximum multipole number to include
    output : str
        the output file name
    """
    # load the data
    r = ConvolvedFFTPower.load(filename)

    # trim to the correct k range
    valid = (r.poles['k'] >= kmin)&(r.poles['k'] <= kmax)
    poles = r.poles.data[valid]

    # subtract shot noise
    poles['power_0'] -= r.attrs['shotnoise']

    # compute the power
    power = []
    for stat in stats:
        if stat in ['P0', 'P2']:
            power.append(poles['power_%s' %stat[-1]].real)
        elif stat in ['P0_sysfree']:
            power.append(poles['power_0'].real + 2./5. * poles['power_2'].real)
    power = numpy.vstack(power).T

    # make the data array
    dtype = [('k', 'f8'), ('power', 'f8')]
    data = numpy.zeros_like(power, dtype=dtype)

    # store k and P(k)
    data['power'] = power[:]
    data['k'] = numpy.repeat(poles['k'][:,None], data.shape[1], axis=1)

    # initialize the power measurements object
    measurements = PowerMeasurements.from_array(names, data)

    # write to file
    if not quiet:
        print(f'saving {output}...')
    measurements.to_plaintext(output)
