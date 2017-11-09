#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

NWALKERS = 50
ITERATIONS = 500

kmins = [1e-4, 0.01]
ZBINS = [(0.8, 2.2), (0.5, 3.0)]
stats = [['P0', 'P2'], ['P0_sysfree']]

@parametrize({'kmin': kmins, 'sample':['N', 'S'], 'stats':stats})
def add_commands(kmin, sample, stats):
    VERSION = 'v1.9f'
    HASHES = ['bba5aabfa6', '983c59a111']
    PARAMS = " ".join(['b1', 'sigma_fog', 'f_nl'])
    stats = " ".join(stats)

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for i, hashstr in enumerate(HASHES):
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc --kmin {kmin} -f {filename} --vary {PARAMS} --stats {stats} -i {ITERATIONS} -w {NWALKERS} --overwrite"

        tag = {'kmin':kmin, 'sample':sample, 'stats':stats, 'zbin':ZBINS[i]}
        RSDFitRunner.register(command, tag=tag)


if __name__ == '__main__':
    add_commands()
    RSDFitRunner.execute()
