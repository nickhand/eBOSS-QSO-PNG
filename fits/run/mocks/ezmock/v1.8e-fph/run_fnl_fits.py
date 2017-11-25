#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os
from argparse import Action

ITERATIONS = 500
ZBINS = [(0.8, 2.2)]
all_params = ['b1', 'sigma_fog', 'f_nl']

class BoxNumber(Action):
    def __call__(self, parser, namespace, box, option_string=None):
        global all_params
        if namespace.vary_shot_noise:
            all_params += ['N']
        add_commands(box=box)
        setattr(namespace, self.dest, box)

ZBINS = [(0.8, 2.2)]

@parametrize({'sample':['N', 'S']})
def add_commands(box, sample):
    VERSION = 'v1.8e-fph'
    HASHES = ['bba5aabfa6', '8ba79d78df', '0ef5c14a14']
    params = " ".join(all_params)

    dirname = os.path.join(EBOSS_SPECTRA, 'mocks', 'ezmock', VERSION)
    for i, hashstr in enumerate(HASHES):

        filename = os.path.join(dirname, f'poles_zevoEZmock_{VERSION}_QSO-{sample}_{box:04d}-{hashstr}.json')
        command = f"eboss-qso-fit nlopt -f {filename} --vary {params} --stats P0 P2 -i {ITERATIONS} --kmax 0.3"
        RSDFitRunner.register(command, tag={'sample':sample, 'zbin':ZBINS[0]})


if __name__ == '__main__':

    # the box number
    RSDFitRunner.update_preparser('--vary-shot-noise', choices=[0,1], type=int, required=True)
    RSDFitRunner.update_preparser('--box', required=True, type=int, action=BoxNumber)
    RSDFitRunner.execute()
