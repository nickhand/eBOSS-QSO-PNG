
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss
import os
import argparse

def main(ns):

    # load the randoms
    randoms = eboss.read_randoms(ns.version, ns.sample)

    # subsample?
    if ns.subsample > 1:
        randoms = randoms[::ns.subsample]

    # trim redshift range
    r = eboss.trim_redshift_range(randoms, zmin=ns.zmin, zmax=ns.zmax)

    # do the pair count
    redges = numpy.logspace(0, 4, 500)
    result = SurveyDataPairCount('2d', r, redges, eboss.fidcosmo, Nmu=100, ra='RA', dec='DEC', redshift='Z')

    # and save!
    kws = {'subsample':ns.subsample, 'zmin':zmin, 'zmax':zmax, 'N':r.csize, 'redges':'logspace(0,4,500)'}
    eboss.save_RR_paircount(result, ns.sample, ns.version, **kws)

if __name__ == '__main__':

    desc = 'compute the pair count of an eBOSS QSO randoms file'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str, choices=['N', 'S'], help=h, required=True)

    h = 'the version to load'
    group.add_argument('--version', type=str, choices=eboss.DATA_VERSIONS, help=h, required=True)

    h = 'the minimum redshift to include'
    group.add_argument('--zmin', type=float, help=h, required=True)

    h = 'the maximum redshift to include'
    group.add_argument('--zmax', type=float, help=h, required=True)

    h = 'the subsample factor'
    parser.add_argument('--subsample', type=int, default=1, help=h)

    # and go!
    main(parser.parse_args())
