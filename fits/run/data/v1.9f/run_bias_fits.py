#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner
from eboss_qso import EBOSS_SPECTRA
import os
from argparse import Action

class RegisterSample(Action):
    def __call__(self, parser, namespace, sample, option_string=None):
        add_commands(sample)
        setattr(namespace, self.dest, sample)

def add_commands(sample):
    VERSION = 'v1.9f'
    HASHES = ['09ccf0e4e0', '3e399248c3', '1f5e5a84d5', '6a332d6423', '1ec7564db2', '30a1823158']

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for hashstr in HASHES:
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc -f {filename} --vary b1 sigma_fog --stats P0 P2 -i 1000 -w 10"
        RSDFitRunner.register(command)


if __name__ == '__main__':

    # the sample
    RSDFitRunner.update_preparser('--sample', choices=['N', 'S'], required=True, type=str, action=RegisterSample)

    # run
    RSDFitRunner.execute()
