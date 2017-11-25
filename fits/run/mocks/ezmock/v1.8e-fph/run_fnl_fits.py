#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

ITERATIONS = 500
ZBINS = [(0.8, 2.2)]
VERSION = 'v1.8e-fph'
HASHES = ['bba5aabfa6', '8ba79d78df', '0ef5c14a14']
p = {'0ef5c14a14': 1.0, 'bba5aabfa6':None, '8ba79d78df':1.6}

@parametrize({'sample':['N', 'S'], 'hashstr':HASHES})
def add_commands(sample, hashstr, box, vary_shot_noise=True):

    # determine the params we are fitting
    all_params = ['b1', 'sigma_fog', 'f_nl']
    if vary_shot_noise:
        all_params += ['N']
    params = " ".join(all_params)

    # filename of spectra we are fitting
    dirname = os.path.join(EBOSS_SPECTRA, 'mocks', 'ezmock', VERSION)
    filename = os.path.join(dirname, f'poles_zevoEZmock_{VERSION}_QSO-{sample}_{box:04d}-{hashstr}.json')

    # make the command
    command = f"eboss-qso-fit nlopt -f {filename} --vary {params} --stats P0 P2 -i {ITERATIONS} --kmax 0.3 --overwrite"
    RSDFitRunner.register(command, tag={'sample':sample, 'p':p[hashstr], 'zbin':ZBINS[0]})


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--vary-shot-noise', choices=[0,1], type=int, required=True)
    parser.add_argument('--box', required=True, type=int)

    ns, args = parser.parse_known_args()

    add_commands(vary_shot_noise=ns.vary_shot_noise, box=ns.box)
    RSDFitRunner.execute(args=args)
