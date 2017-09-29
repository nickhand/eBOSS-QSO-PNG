
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss

import os
import argparse

setup_logging()

def main(ns):

    # initialize the task manager
    with TaskManager(ns.cpus_per_task, use_all_cpus=True) as tm:

        # load the randoms
        randoms = eboss.read_ezmock_randoms(ns.sample, ns.version)
        eboss.finalize_ezmock(randoms, eboss.fidcosmo, P0_FKP=ns.P0_FKP)

        for box_num in tm.iterate(range(ns.start, ns.stop, ns.step)):

            # load the data
            data = eboss.read_ezmock_data(box_num, ns.sample, ns.version, ns.subversion)
            eboss.finalize_ezmock(data, eboss.fidcosmo, P0_FKP=ns.P0_FKP)

            # combine data and randoms into the FKP source
            fkp = FKPCatalog(data=data, randoms=randoms, BoxPad=0.1, use_cache=True)

            # mesh kwargs
            mesh_kwargs = {'Nmesh':1024, 'interlaced':True, 'window':'tsc', 'dtype':'f8'}

            # compute unweighted results
            if ns.do_unweighted:
                unweighted_mesh = fkp.to_mesh(nbar='NZ', fkp_weight='FKPWeight', comp_weight='Weight', **mesh_kwargs)

                # run
                result = ConvolvedFFTPower(first=unweighted_mesh, poles=[0,2], dk=0.005, kmin=0.)

                # save
                meta = {'p':None, 'zmin':0.9, 'zmax':2.2, 'P0_FKP':ns.P0_FKP}
                eboss.save_ezmock_spectra(result, box_num, ns.sample, ns.version, ns.subversion, **meta)

            # the bias weight for the first field
            fkp['data/BiasWeight'] = data['FKPWeight'] * eboss.bias_weight(data['Z'], eboss.fidcosmo)
            fkp['randoms/BiasWeight'] = randoms['FKPWeight'] * eboss.bias_weight(randoms['Z'], eboss.fidcosmo)

            # the fnl weight for the second field
            fkp['data/FnlWeight'] = data['FKPWeight'] * eboss.fnl_weight(data['Z'], p=ns.p)
            fkp['randoms/FnlWeight'] = randoms['FKPWeight'] * eboss.fnl_weight(randoms['Z'], p=ns.p)

            # convert to mesh
            mesh1 = fkp.to_mesh(nbar='NZ', fkp_weight='BiasWeight', comp_weight='Weight', **mesh_kwargs)
            mesh2 = fkp.to_mesh(nbar='NZ', fkp_weight='FnlWeight', comp_weight='Weight', **mesh_kwargs)

            # run
            result = ConvolvedFFTPower(first=mesh1, second=mesh2, poles=[0,2], dk=0.005, kmin=0.)

            # save
            meta = {'p':ns.p, 'zmin':0.9, 'zmax':zmax, 'P0_FKP':ns.P0_FKP}
            eboss.save_ezmock_spectra(result, box_num, ns.sample, ns.version, ns.subversion, **meta)


if __name__ == '__main__':

    desc = 'compute the redshift weighted power spectra of the eBOSS QSO EZ mocks'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'number of cpus per task'
    group.add_argument('cpus_per_task', type=int, help=h)

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str, choices=['N', 'S'], help=h, required=True)

    h = 'the version to load'
    group.add_argument('--version', type=str, choices=eboss.EZMOCK_VERSIONS, help=h, required=True)

    h = 'the version to load'
    group.add_argument('--subversion', type=str, choices=eboss.EZMOCK_SUBVERSIONS, help=h, required=True)

    h = 'the value of p to use'
    group.add_argument('--p', type=float, help=h, required=True)

    h = 'the P0 FKP version to use'
    parser.add_argument('--P0_FKP', type=float, default=3e4, help=h)

    h = 'whether to compute unweighted results'
    parser.add_argument('--do-unweighted', action='store_true', help=h)

    h = 'the start box number'
    parser.add_argument('--start', default=1, type=int, help=h)

    h = 'the stop box number'
    parser.add_argument('--stop', default=1000, type=int, help=h)

    h = 'the step box number'
    parser.add_argument('--step', default=1, type=int, help=h)

    # and go!
    main(parser.parse_args())
