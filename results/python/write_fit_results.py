from eboss_qso.fits import load_data_results, load_joint_data_results
from eboss_qso.fits import load_ezmock_results
import os
import numpy as np
import shutil


def load_fit_results(kind, sample, kmin, z_weighted, p):

    if kind == 'Mocks':
        d = load_ezmock_results('v1.8e-fph', str(sample),
                                kmin+'-0.3', 'basemodel-N-fnl', z_weighted, p=p)
    else:
        if sample == 'N+S':
            r = load_joint_data_results(kmin, z_weighted, p)
            r.burnin = int(r.iterations//2)
        else:

            d = load_data_results(
                'v1.9f', str(sample), kmin+'-0.3', 'basemodel-N-fnl', '0.8-2.2', z_weighted, p=p)
            d.set_fit_results()
            r = d.results
            r.burnin = int(r.iterations//2)

        d = np.vstack([r[name].flat_trace for name in r.free_names]).T
        d = d.flatten().view([(name, 'f8') for name in r.free_names])
    return d


def ensure_dir(*dirnames, remove=True):
    for dirname in dirnames:
        if remove and os.path.exists(dirname):
            shutil.rmtree(dirname)

        try:
            os.makedirs(dirname)
        except:
            pass


# make output directory
output_dir = os.path.join("..", "products", "bestfit_params")
ensure_dir(output_dir)

for kind in ['Mocks', 'Data']:

    samples = ['N', 'S']
    if kind == 'Data':
        samples += ['N+S']

    for sample in samples:

        if kind == 'Mocks':
            dirname = os.path.join(output_dir, 'ezmock')
        else:
            dirname = os.path.join(output_dir, 'data')
        ensure_dir(dirname, remove=False)

        for kmin in ['0.0001', '0.005']:
            if kmin == '0.0001':
                k_tag = 'all_bins'
            else:
                k_tag = 'first_bin_excluded'

            dirname_ = os.path.join(dirname, k_tag)
            ensure_dir(dirname_, remove=False)

            results = [(True, 1.0), (True, 1.6), (False, 1.0), (False, 1.6)]
            for (z_weighted, p) in results:

                if not z_weighted:
                    dirname__ = os.path.join(dirname_, 'fkp_only_%.1f' % p)
                else:
                    dirname__ = os.path.join(dirname_, 'weighted_%.1f' % p)
                ensure_dir(dirname__, remove=False)

                # get results
                r = load_fit_results(kind, sample, kmin, z_weighted, p)

                if sample != 'N+S':
                    tag = '%sgc' % sample.lower()
                else:
                    tag = "ngc+sgc"
                filename = os.path.join(dirname__, "params_%s.npy" % tag)
                np.save(filename, r)
