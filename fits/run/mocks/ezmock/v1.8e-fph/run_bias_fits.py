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
def add_commands(sample, hashstr, box, vary_shot_noise=True, cov='analytic',
                    use_temp_files=False, overwrite=False, kmin=1e-4, kmax=0.3):

    # determine the params we are fitting
    all_params = ['b1', 'sigma_fog']
    if vary_shot_noise:
        all_params += ['N']
    params = " ".join(all_params)

    # filename of spectra we are fitting
    dirname = os.path.join(EBOSS_SPECTRA, 'mocks', 'ezmock', VERSION)
    filename = os.path.join(dirname, f'poles_zevoEZmock_{VERSION}_QSO-{sample}_{box:04d}-{hashstr}.json')

    # make the command
    command = f"eboss-qso-fit nlopt -f {filename} --vary {params} --stats P0 P2"
    command += f" -i {ITERATIONS} --kmin {kmin} --kmax {kmax} --cov {cov}"
    if use_temp_files:
        command += " --use-temp-files"
    if overwrite:
        command += " --overwrite"

    # register it
    tag = {'sample':sample, 'p':p[hashstr], 'zbin':ZBINS[0]}
    RSDFitRunner.register(command, tag=tag)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--box', required=True, type=int)
    parser.add_argument('--vary-shot-noise', choices=[0,1], type=int, required=True)
    parser.add_argument('--cov', choices=['mock', 'analytic'], required=True)
    parser.add_argument('--overwrite', action='store_true', default=False)
    parser.add_argument('--kmin', type=float, default=1e-4)
    parser.add_argument('--kmax', type=float, default=0.3)

    ns, args = parser.parse_known_args()

    add_commands(**vars(ns))
    RSDFitRunner.execute(args=args)
