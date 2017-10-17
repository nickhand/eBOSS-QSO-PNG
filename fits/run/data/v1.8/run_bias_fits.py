#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

ZBINS = [(0.9, 1.2), (1.2, 1.5), (1.5, 1.8), (1.8, 2.2)]

@parametrize({'sample':['N', 'S']})
def add_commands(sample):
    VERSION = 'v1.8'
    HASHES = ['a34d7840d9', 'd96a267921', '89f0eee23b', 'f8a1e4e844']

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for i, hashstr in enumerate(HASHES):
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc -f {filename} --vary b1 sigma_fog --stats P0 P2 -i 1000 -w 10"
        RSDFitRunner.register(command, tag={'sample':sample, 'zbin':ZBINS[i]})


if __name__ == '__main__':
    add_commands()
    RSDFitRunner.execute()
