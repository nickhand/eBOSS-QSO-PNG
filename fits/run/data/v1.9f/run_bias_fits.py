#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

ITERATIONS = 1000
WALKERS = 10
ZBINS = [(0.9, 1.2), (1.2, 1.5), (1.5, 1.8), (1.8, 2.2), (0.9, 2.2)]

@parametrize({'sample':['N', 'S']})
def add_commands(sample, vary_shot_noise=True, cov='analytic',
                  overwrite=False, kmin=1e-4, kmax=0.3):

    VERSION = 'v1.9f'
    HASHES = ['3e399248c3', '1f5e5a84d5', '6a332d6423', '1ec7564db2', '8b9aab839d']

    # determine the params we are fitting
    all_params = ['b1', 'sigma_fog']
    if vary_shot_noise:
        all_params += ['N']
    params = " ".join(all_params)

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    for i, hashstr in enumerate(HASHES):

        # the data file
        filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')

        # make the command
        command = f"eboss-qso-fit mcmc -f {filename} --vary {params} --stats P0 P2"
        command += f" -i {ITERATIONS} -w {WALKERS} --kmin {kmin} --kmax {kmax} --cov {cov}"
        if overwrite:
            command += " --overwrite"

        # register it
        tag = {'sample':sample, 'zbin':ZBINS[i]}
        RSDFitRunner.register(command, tag=tag)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--vary-shot-noise', choices=[0,1], type=int, required=True)
    parser.add_argument('--overwrite', action='store_true', default=False)
    parser.add_argument('--kmin', type=float, default=1e-4)
    parser.add_argument('--kmax', type=float, default=0.3)

    ns, args = parser.parse_known_args()

    add_commands(**vars(ns))
    RSDFitRunner.execute(args=args)
