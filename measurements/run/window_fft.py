
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss
import os
import argparse

setup_logging()

def select_subsample(insize, outsize):

    index = numpy.zeros(insize, dtype=bool)
    valid = numpy.random.choice(range(insize), replace=False, size=outsize)
    index[valid] = True
    return index


def main(ns):

    # load the randoms
    randoms = eboss.read_randoms(ns.sample, ns.version)

    # set the seed
    numpy.random.seed(42*(randoms.comm.rank+1000))

    # trim redshift range
    r = eboss.trim_redshift_range(randoms, zmin=ns.zmin, zmax=ns.zmax)

    # finalize
    eboss.finalize_ezmock(r, eboss.ezmock_cosmo)

    # FKP
    fkp = FKPCatalog(data=None, randoms=r, BoxPad=0.1)

    # to mesh
    mesh_kwargs = {'Nmesh':1024, 'interlaced':True, 'window':'tsc', 'dtype':'f8'}
    mesh = fkp.to_mesh(nbar='NZ', fkp_weight='FKPWeight', comp_weight='Weight', **mesh_kwargs)

    # power
    ells = [0, 2, 4, 6, 8, 10]
    result = ConvolvedFFTPower(first=mesh, poles=ells, dk=0.005, kmin=0.)

    # and save!
    kws = {'zmin':ns.zmin, 'zmax':ns.zmax}
    eboss.save_RR_poles(result, ns.sample, ns.version, **kws)

if __name__ == '__main__':

    desc = 'compute multipoles of an eBOSS QSO randoms file'
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

    # and go!
    main(parser.parse_args())
