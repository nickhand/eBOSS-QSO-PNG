#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

ITERATIONS = 500
ZBINS = [(0.8, 2.2)]
VERSION = 'v1.8e-fph'

# hash values for different p
HASHES = ['bba5aabfa6', '8ba79d78df', '0ef5c14a14']
p = {'0ef5c14a14': 1.0, 'bba5aabfa6':None, '8ba79d78df':1.6}

# params we want to vary
fixed_params = ['b1', 'sigma_fog', 'f', 'sigma8_z']
all_params = fixed_params + ['alpha_par', 'alpha_perp']

@parametrize({'sample':['N', 'S'], 'hashstr':HASHES, 'params':[all_params, fixed_params]})
def add_commands(sample, hashstr, params, box, vary_shot_noise=True):

    # determine the params we are fitting
    if vary_shot_noise:
        params = params + ['N']
    params = " ".join(params)

    # filename of spectra we are fitting
    dirname = os.path.join(EBOSS_SPECTRA, 'mocks', 'ezmock', VERSION)
    filename = os.path.join(dirname, f'poles_zevoEZmock_{VERSION}_QSO-{sample}_{box:04d}-{hashstr}.json')

    # make the command
    command = f"eboss-qso-fit nlopt -f {filename} --vary {params} --stats P0 P2 -i {ITERATIONS} --kmax 0.3 --overwrite"
    tag = {'sample':sample, 'p':p[hashstr], 'params':params, 'zbin':ZBINS[0]}
    RSDFitRunner.register(command, tag=tag)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--vary-shot-noise', choices=[0,1], type=int, required=True)
    parser.add_argument('--box', required=True, type=int)

    ns, args = parser.parse_known_args()

    add_commands(vary_shot_noise=ns.vary_shot_noise, box=ns.box)
    RSDFitRunner.execute(args=args)
