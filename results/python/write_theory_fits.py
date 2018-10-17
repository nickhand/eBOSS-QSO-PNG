
from eboss_qso.fits import load_data_results
from eboss_qso.fits import load_ezmock_driver
import os
import numpy as np
import shutil
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank


def ensure_dir(*dirnames, remove=True):
    for dirname in dirnames:
        if rank == 0 and remove and os.path.exists(dirname):
            shutil.rmtree(dirname)

        try:
            os.makedirs(dirname)
        except:
            pass

        comm.barrier()


# make output directory
output_dir = os.path.join("..", "products", "theory_fits")
ensure_dir(output_dir)


cols = ['k', 'power_0', 'power_2']
for sample in ['N', 'S']:
    tag = sample.lower()

    # output dirs
    mock_dir = os.path.join(output_dir, 'ezmock')
    data_dir = os.path.join(output_dir, 'data')
    ensure_dir(mock_dir, data_dir, remove=False)

    for kmin in ['0.0001', '0.005']:
        if kmin == '0.0001':
            k_tag = 'all_bins'
        else:
            k_tag = 'first_bin_excluded'

        mock_dir_ = os.path.join(mock_dir, k_tag)
        data_dir_ = os.path.join(data_dir, k_tag)
        ensure_dir(mock_dir_, data_dir_, remove=False)

        results = [(True, 1.0), (True, 1.6), (False, 1.0), (False, 1.6)]
        for j in range(rank, comm.size, comm.size):
            (z_weighted, p) = results[j]

            if not z_weighted:
                mock_dir__ = os.path.join(mock_dir_, 'fkp_only_%.1f' % p)
                data_dir__ = os.path.join(data_dir_, 'fkp_only_%.1f' % p)
            else:
                mock_dir__ = os.path.join(mock_dir_, 'weighted_%.1f' % p)
                data_dir__ = os.path.join(data_dir_, 'weighted_%.1f' % p)
            ensure_dir(mock_dir__, data_dir__, remove=False)

            # mocks
            for i in range(1, 1001):
                print(i)
                d = load_ezmock_driver(i, 'v1.8e-fph', sample,
                                        f'{kmin}-0.3', 'basemodel-N-fnl',
                                        z_weighted, p=p)
                d.results.burnin = d.results.iterations//2
                d.set_fit_results()

                slices = d.data.flat_slices
                theory = d.combined_model

                filename = os.path.join(
                    mock_dir__, "mock_%04d_%sgc.dat" % (i, tag))
                with open(filename, 'wb') as ff:
                    ff.write(("# k P0 P2\n").encode())
                    ff.write(("# reduced chi2 = %.4f\n" %
                                d.reduced_chi2()).encode())
                    out = np.vstack([d.data[0].k] + [theory[sl]
                                                        for sl in slices]).T
                    np.savetxt(ff, out, fmt='%-15.4e')

            # data
            d = load_data_results('v1.9f', sample, f'{kmin}-0.3',
                                    'basemodel-N-fnl', '0.8-2.2', z_weighted, p=p)
            d.results.burnin = d.results.iterations//2
            d.set_fit_results()

            slices = d.data.flat_slices
            theory = d.combined_model

            filename = os.path.join(data_dir__, "data_%sgc.dat" % tag)
            with open(filename, 'wb') as ff:
                ff.write(("# k P0 P2\n").encode())
                ff.write(("# reduced chi2 = %.4f\n" %
                            d.reduced_chi2()).encode())
                out = np.vstack([d.data[0].k] + [theory[sl]
                                                    for sl in slices]).T
                np.savetxt(ff, out, fmt='%-15.4e')
