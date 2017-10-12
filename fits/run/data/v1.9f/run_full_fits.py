#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

ZBINS = [(0.9, 2.2), (0.5, 3.0)]
fixed_params = ['b1', 'sigma_fog', 'f', 'sigma8_z', 'f_nl']
all_params = fixed_params + ['alpha_par', 'alpha_perp']
stats = [['P0', 'P2'], ['P0_sysfree']]

@parametrize({'sample':['N', 'S'], 'params': [all_params, fixed_params], 'stats':stats})
def add_commands(sample, stats, params):
    VERSION = 'v1.9f'
    HASHES = ['8b9aab839d', '983c59a111']
    stats = " ".join(stats); params = " ".join(params)

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for i, hashstr in enumerate(HASHES):
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc -f {filename} --vary {params} --stats {stats} -i 1000 -w 20"
        tag = {'sample':sample, 'params':params, 'stats':stats, 'zbin':ZBINS[i]}
        RSDFitRunner.register(command, tag=tag)

if __name__ == '__main__':
    add_commands()
    RSDFitRunner.execute()
