from eboss_qso.measurements import load_ezmock_spectra, load_data_spectra
import os
import numpy as np
import shutil


def ensure_dir(*dirnames, remove=True):
    for dirname in dirnames:
        if remove and os.path.exists(dirname):
            shutil.rmtree(dirname)

        if not os.path.exists(dirname):
            os.makedirs(dirname)


# make output directory
output_dir = os.path.join("..", "products", "power_spectra")
ensure_dir(output_dir)

cols = ['k', 'power_0', 'power_2']
for sample in ['N', 'S']:
    tag = sample.lower()

    # output dirs
    mock_dir = os.path.join(output_dir, 'ezmock')
    data_dir = os.path.join(output_dir, 'data')
    ensure_dir(mock_dir, data_dir, remove=False)

    for p in [None, 1.6, 1.0]:

        if p is None:
            mock_dir_ = os.path.join(mock_dir, 'fkp_only')
            data_dir_ = os.path.join(data_dir, 'fkp_only')
        else:
            mock_dir_ = os.path.join(mock_dir, 'weighted_%.1f' % p)
            data_dir_ = os.path.join(data_dir, 'weighted_%.1f' % p)
        ensure_dir(mock_dir_, data_dir_, remove=False)

        # load the mocks
        P0 = load_ezmock_spectra(
            'v1.8e-fph', sample, p=p, average=False, ell=0)
        P2 = load_ezmock_spectra(
            'v1.8e-fph', sample, p=p, average=False, ell=2)

        for i in range(1, len(P0)+1):

            x = P0[i-1]
            y = P2[i-1]
            d = np.vstack([x['k'], x['power_0'].real,
                           y['power_2'].real, x['modes']]).T

            filename = os.path.join(mock_dir_, "mock_%04d_%sgc.dat" % (i, tag))
            with open(filename, 'wb') as ff:
                ff.write(("# k P0 P2 modes\n").encode())
                np.savetxt(ff, d, fmt='%-15.4e')

        # mean mock spectra
        x = load_ezmock_spectra('v1.8e-fph', sample, p=p, average=True, ell=0)
        y = load_ezmock_spectra('v1.8e-fph', sample, p=p, average=True, ell=2)
        d = np.vstack([x['k'], x['power_0'].real,
                       y['power_2'].real, x['modes']]).T

        filename = os.path.join(mock_dir_, "mean_%sgc.dat" % tag)
        with open(filename, 'wb') as ff:
            ff.write(("# k P0 P2 modes\n").encode())
            np.savetxt(ff, d, fmt='%-15.4e')

        # mean data
        x = load_data_spectra(
            'v1.9f', sample, p=p, zmin=0.8, zmax=2.2, focal_weights=True, ell=0).poles.data
        y = load_data_spectra(
            'v1.9f', sample, p=p, zmin=0.8, zmax=2.2, focal_weights=True, ell=2).poles.data
        d = np.vstack([x['k'], x['power_0'].real,
                       y['power_2'].real, x['modes']]).T

        filename = os.path.join(data_dir_, "data_%sgc.dat" % tag)
        with open(filename, 'wb') as ff:
            ff.write(("# k P0 P2 modes\n").encode())
            np.savetxt(ff, d, fmt='%-15.4e')
