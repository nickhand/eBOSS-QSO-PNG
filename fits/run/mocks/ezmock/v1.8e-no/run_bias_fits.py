#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os
from argparse import Action

class BoxNumber(Action):
    def __call__(self, parser, namespace, box, option_string=None):
        add_commands(box=box)
        setattr(namespace, self.dest, box)

ZBINS = [(0.8, 2.2)]

@parametrize({'sample':['N', 'S']})
def add_commands(box, sample):
    VERSION = 'v1.8e-no'
    HASHES = ['bba5aabfa6'] # no z-weights

    dirname = os.path.join(EBOSS_SPECTRA, 'mocks', 'ezmock', VERSION)
    for i, hashstr in enumerate(HASHES):
        filename = os.path.join(dirname, f'poles_zevoEZmock_{VERSION}_QSO-{sample}_{box:04d}-{hashstr}.json')
        command = f"eboss-qso-fit nlopt -f {filename} --vary b1 sigma_fog --stats P0 P2 -i 500"
        RSDFitRunner.register(command, tag={'sample':sample, 'zbin':ZBINS[i]})


if __name__ == '__main__':

    # the box number
    RSDFitRunner.update_preparser('--box', required=True, type=int, action=BoxNumber)
    RSDFitRunner.execute()
