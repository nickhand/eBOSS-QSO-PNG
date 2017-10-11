#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner
from eboss_qso import EBOSS_SPECTRA
import os
from argparse import Action

class RegisterSample(Action):
    def __call__(self, parser, namespace, sample, option_string=None):
        stats = getattr(namespace, 'stats', None)
        if stats is None:
            raise ValueError("please specify 'stats' first on command line")

        # the parameters 
        PARAMS = ['b1', 'sigma_fog', 'f_nl']

        # register the commands
        add_commands(sample, " ".join(stats), " ".join(PARAMS))
        setattr(namespace, self.dest, sample)

def add_commands(sample, stats, params):
    VERSION = 'v1.9f'
    HASHES = ['8b9aab839d', '983c59a111']

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for hashstr in HASHES:
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')
        command = f"eboss-qso-fit mcmc -f {filename} --vary {params} --stats {stats} -i 1000 -w 20"
        RSDFitRunner.register(command)


if __name__ == '__main__':

    # the statistics
    choices = ['P0', 'P2', 'P0_sysfree']
    RSDFitRunner.update_preparser('--stats', choices=choices, required=True, nargs='*', type=str)

    # the sample
    RSDFitRunner.update_preparser('--sample', choices=['N', 'S'], required=True, type=str, action=RegisterSample)

    # run
    RSDFitRunner.execute()
