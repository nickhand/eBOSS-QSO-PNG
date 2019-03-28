
from eboss_qso.fits import load_data_results
from eboss_qso.fits import load_ezmock_driver
from eboss_qso.fits.covariance import compute_ezmock_covariance
from eboss_qso.measurements import load_ezmock_spectra, load_data_spectra
from pyRSD.rsdfit.data import PowerMeasurements
import os
import numpy as np
import shutil


def check_covariance(d, sample, p, kmin, z_weighted):

    if d.data.covariance_matrix.N != len(d.data.combined_k):
        output = d.data.covariance
        C = compute_ezmock_covariance('v1.8e-fph',
                                      sample,
                                      [0, 2],
                                      kmin=0.0001,
                                      kmax=0.3,
                                      p=p)
        C.to_plaintext(output)

    d = load_ezmock_driver(i, 'v1.8e-fph', sample,
                           f'{kmin}-0.3', 'basemodel-N-fnl',
                           z_weighted, p=p)
    return d


def write_spectra(output, loader, **kwargs):

    power = []
    for ell in [0, 2]:
        r = loader(ell=ell, subtract_shot_noise=True,
                   **kwargs)

        # trim to the correct k range
        valid = (r.poles['k'] >= 0.0001) & (r.poles['k'] <= 0.3)
        poles = r.poles.data[valid]

        # add the power array
        power.append(poles['power_%s' % ell].real)

    # convert to array
    power = np.vstack(power).T

    # make the data array
    dtype = [('k', 'f8'), ('power', 'f8')]
    data = np.zeros_like(power, dtype=dtype)

    # store k and P(k)
    data['power'] = power[:]
    data['k'] = np.repeat(poles['k'][:, None], data.shape[1], axis=1)

    # initialize the power measurements object
    measurements = PowerMeasurements.from_array(['pole_0', 'pole_2'], data)
    measurements.to_plaintext(output)


def ensure_dir(*dirnames, remove=True):
    for dirname in dirnames:
        if remove and os.path.exists(dirname):
            shutil.rmtree(dirname)

        try:
            os.makedirs(dirname)
        except:
            pass


# make output directory
output_dir = os.path.join("..", "products", "theory_fits")
ensure_dir(output_dir)


cols = ['k', 'power_0', 'power_2']
for sample in ['N', 'S']:
    tag = sample.lower()

    # output dirs
    data_dir = os.path.join(output_dir, 'data')
    ensure_dir(data_dir, remove=False)

    for kmin in ['0.0001', '0.005']:
        print(kmin)
        if kmin == '0.0001':
            k_tag = 'all_bins'
        else:
            k_tag = 'first_bin_excluded'

        data_dir_ = os.path.join(data_dir, k_tag)
        ensure_dir(data_dir_, remove=False)

        results = [(True, 1.0), (True, 1.6), (False, 1.0), (False, 1.6)]
        for (z_weighted, p) in results:
            if not z_weighted:
                data_dir__ = os.path.join(data_dir_, 'fkp_only_%.1f' % p)
            else:
                data_dir__ = os.path.join(data_dir_, 'weighted_%.1f' % p)
            ensure_dir(data_dir__, remove=False)

            # data
            d = load_data_results('v1.9f', sample, f'{kmin}-0.3',
                                  'basemodel-N-fnl', '0.8-2.2', z_weighted, p=p, ells=[0])

            # fix kmin
            if kmin == '0.0001' and d.data[0].k[0] > 0.005:
                output = d.data.params['data_file'].value
                write_spectra(output, load_data_spectra,
                              version='v1.9f', sample=sample, p=p)
                d = load_data_results('v1.9f', sample, f'{kmin}-0.3',
                                      'basemodel-N-fnl', '0.8-2.2', z_weighted, p=p, ells=[0])

            if kmin == '0.0001':
                assert d.data[0].k[0] < 0.005

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
