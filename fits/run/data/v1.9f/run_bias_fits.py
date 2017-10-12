#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

ZBINS = [(0.5, 0.9), (0.9, 1.2), (1.2, 1.5), (1.5, 1.8), (1.8, 2.2), (2.2, 3.0)]

@parametrize({'sample':['N', 'S']})
def add_commands(sample):
    VERSION = 'v1.9f'
    HASHES = ['09ccf0e4e0', '3e399248c3', '1f5e5a84d5', '6a332d6423', '1ec7564db2', '30a1823158']

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for i, hashstr in enumerate(HASHES):
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc -f {filename} --vary b1 sigma_fog --stats P0 P2 -i 1000 -w 10"
        RSDFitRunner.register(command, tag={'sample':sample, 'zbin':ZBINS[i]})


if __name__ == '__main__':

    add_commands()
    RSDFitRunner.execute()
