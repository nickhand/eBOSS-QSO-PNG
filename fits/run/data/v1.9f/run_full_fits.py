#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

NWALKERS = 200
ITERATIONS = 2000
ZBINS = [(0.8, 2.2)]
VERSION = 'v1.9f'
HASHES = ['bba5aabfa6', '8ba79d78df', '0ef5c14a14']
p = {'0ef5c14a14': 1.0, 'bba5aabfa6':None, '8ba79d78df':1.6}

# params we want to vary
fixed_params = ['b1', 'sigma_fog', 'f', 'sigma8_z', 'f_nl']
all_params = fixed_params + ['alpha_par', 'alpha_perp']


@parametrize({'sample':['N', 'S'], 'hashstr':HASHES, 'params':[all_params, fixed_params]})
def add_commands(sample, hashstr, params, vary_shot_noise=True, cov='analytic',
                    use_temp_files=False, overwrite=False, kmin=1e-4, kmax=0.3):

    # determine the params we are fitting
    if vary_shot_noise:
        params = params + ['N']
    params = " ".join(params)

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')

    # make the command
    command = f"eboss-qso-fit mcmc -f {filename} --vary {params} --stats P0 P2"
    command += f" -i {ITERATIONS} -w {NWALKERS} --kmin {kmin} --kmax {kmax} --cov {cov}"
    if use_temp_files:
        command += " --use-temp-files"
    if overwrite:
        command += " --overwrite"

    # and register it
    tag = {'sample':sample, 'p':p[hashstr], 'zbin':ZBINS[0], 'params':params}
    RSDFitRunner.register(command, tag=tag)

if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--vary-shot-noise', choices=[0,1], type=int, required=True)
    parser.add_argument('--cov', choices=['mock', 'analytic'], required=True)
    parser.add_argument('--overwrite', action='store_true', default=False)
    parser.add_argument('--kmin', type=float, default=1e-4)
    parser.add_argument('--kmax', type=float, default=0.3)

    ns, args = parser.parse_known_args()

    add_commands(**vars(ns))
    RSDFitRunner.execute(args=args)
