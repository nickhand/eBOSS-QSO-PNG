#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner
from eboss_qso import EBOSS_SPECTRA
import os

VERSION = 'v1.8'
SAMPLES = ['N', 'S']
HASHES = ['a34d7840d9', 'd96a267921', '89f0eee23b', 'f8a1e4e844']

dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)

for sample in SAMPLES:
    for hashstr in HASHES:
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc -f {filename} --vary b1 sigma_fog --stats P0 P2 -i 1000 -w 10"
        RSDFitRunner.register(command)


if __name__ == '__main__':

    RSDFitRunner.execute()
