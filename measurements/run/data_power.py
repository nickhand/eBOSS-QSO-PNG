
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss

import os
import argparse

setup_logging()

def main(ns):

    # load the data first
    data = eboss.read_data(ns.version, ns.sample, focal_weights=ns.focal_weights)

    # load the randoms
    randoms = eboss.read_randoms(ns.version, ns.sample)

    # compute for every redshift bin
    for (zmin, zmax) in ns.zbins:

        # trim redshift range
        d = eboss.trim_redshift_range(data, zmin=zmin, zmax=zmax)
        r = eboss.trim_redshift_range(randoms, zmin=zmin, zmax=zmax)

        # finalize columns
        eboss.finalize_source(d, eboss.fidcosmo, P0_FKP=ns.P0_FKP)
        eboss.finalize_source(r, eboss.fidcosmo, P0_FKP=ns.P0_FKP)

    # combine data and randoms into the FKP source
    fkp = FKPCatalog(data=d, randoms=r, BoxPad=0.1, use_cache=True)

    # the mesh kwargs to use
    mesh_kwargs = {'Nmesh':1024, 'interlaced':True, 'window':'tsc', 'dtype':'f8'}

    # compute unweighted results
    if ns.do_unweighted:
        unweighted_mesh = fkp.to_mesh(nbar='NZ', fkp_weight='FKPWeight', comp_weight='Weight', **mesh_kwargs)
        result = ConvolvedFFTPower(first=unweighted_mesh, poles=[0,2], dk=0.005, kmin=0.)
        eboss.save_data_spectra(result, ns.sample, ns.version, ns.focal_weights,
                                    p=None, zmin=zmin, zmax=zmax, P0_FKP=ns.P0_FKP)

    # the bias weight for the first field
    fkp['data/BiasWeight'] = d['FKPWeight'] * eboss.bias_weight(d['Z'], eboss.fidcosmo)
    fkp['randoms/BiasWeight'] = r['FKPWeight'] * eboss.bias_weight(r['Z'], eboss.fidcosmo)

    # the fnl weight for the second field
    fkp['data/FnlWeight'] = d['FKPWeight'] * eboss.fnl_weight(d['Z'], p=ns.p)
    fkp['randoms/FnlWeight'] = r['FKPWeight'] * eboss.fnl_weight(r['Z'], p=ns.p)

    # convert to mesh
    mesh1 = fkp.to_mesh(nbar='NZ', fkp_weight='BiasWeight', comp_weight='Weight', **mesh_kwargs)
    mesh2 = fkp.to_mesh(nbar='NZ', fkp_weight='FnlWeight', comp_weight='Weight', **mesh_kwargs)

    # get power and save
    result = ConvolvedFFTPower(first=mesh1, second=mesh2, poles=[0,2], dk=0.005, kmin=0.)
    eboss.save_data_spectra(result, ns.sample, ns.version, ns.focal_weights,
                                p=ns.p, zmin=zmin, zmax=zmax, P0_FKP=ns.P0_FKP)


if __name__ == '__main__':

    desc = 'compute the redshift weighted power spectra of the eBOSS QSO data sample'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str, choices=['N', 'S'], help=h, required=True)

    h = 'the version to load'
    group.add_argument('--version', type=str, choices=eboss.DATA_VERSIONS, help=h, required=True)

    h = 'the redshift bins to compute in the form of ZMIN, ZMAX'
    group.add_argument('--zbins', type=eboss.redshift_range_type, nargs='+', help=h)

    h = 'the value of p to use'
    group.add_argument('--p', type=float, help=h, required=True)

    h = 'the P0 FKP version to use'
    parser.add_argument('--P0_FKP', type=float, default=3e4, help=h)

    h = 'whether to compute unweighted results'
    parser.add_argument('--do-unweighted', action='store_true', help=h)

    h = 'whether to use focal plane weights'
    parser.add_argument('--focal-weights', action='store_true', help=h)

    # and go!
    main(parser.parse_args())
