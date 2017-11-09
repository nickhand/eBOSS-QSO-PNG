#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

NWALKERS = 100
ITERATIONS = 1000

kmins = [1e-4, 0.01]
ZBINS = [(0.8, 2.2), (0.5, 3.0)]
fixed_params = ['b1', 'sigma_fog', 'f', 'sigma8_z', 'f_nl']
all_params = fixed_params + ['alpha_par', 'alpha_perp']
stats = [['P0', 'P2'], ['P0_sysfree']]

@parametrize({'kmin':kmins, 'sample':['N', 'S'], 'params': [all_params, fixed_params], 'stats':stats})
def add_commands(kmin, sample, stats, params):
    VERSION = 'v1.9f'
    HASHES = ['bba5aabfa6', '983c59a111']
    stats = " ".join(stats); params = " ".join(params)

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for i, hashstr in enumerate(HASHES):
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc --kmin {kmin} -f {filename} --vary {params} --stats {stats} -i {ITERATIONS} -w {NWALKERS} --overwrite"
        tag = {'kmin':kmin, 'sample':sample, 'params':params, 'stats':stats, 'zbin':ZBINS[i]}
        RSDFitRunner.register(command, tag=tag)

if __name__ == '__main__':
    add_commands()
    RSDFitRunner.execute()
