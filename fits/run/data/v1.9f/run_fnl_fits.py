#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

NWALKERS = 50
ITERATIONS = 500
ZBINS = [(0.8, 2.2)]
VERSION = 'v1.9f'

# data and theory values for p
p_values = [(None, 1.0),
            (None, 1.6),
            (1.0, 1.0),
            (1.6, 1.6)]

# effective redshifts based on data p value
effective_redshifts = {None: (1.557, 1.557),
                       1.0: (1.728, 1.671),
                       1.6: (1.829, 1.794)}


@parametrize({'sample': ['N', 'S'], 'p_values': p_values})
def add_commands(sample, p_values, vary_shot_noise=True, cov='analytic',
                 use_temp_files=False, overwrite=False, kmin=1e-4, kmax=0.3):

    # the p values to use
    data_p, theory_p = p_values

    # determine the params we are fitting
    all_params = ['b1', 'sigma_fog', 'f_nl']
    if vary_shot_noise:
        all_params += ['N']
    params = " ".join(all_params)

    # the effective redshift
    z_eff = effective_redshifts[data_p]

    # arguments for reading the data
    command = "eboss-qso-fit mcmc"
    command += f" --kind data --version {VERSION} --sample {sample} --ells 0 2"
    command += f" --z_eff {z_eff[0]} {z_eff[1]}"
    command += f" --theory_p {theory_p} --kmin {kmin} --kmax {kmax}"
    if data_p is not None:
        command += f" --data_p {data_p}"

    # theory arguments
    command += f" --vary {params} -i {ITERATIONS} -w {NWALKERS} --cov {cov}"
    if use_temp_files:
        command += " --use-temp-files"
    if overwrite:
        command += " --overwrite"

    # and register it
    weighted = data_p is not None
    tag = {'sample': sample, 'p': theory_p,
           'weighted': weighted, 'zbin': ZBINS[0]}
    RSDFitRunner.register(command, tag=tag)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--vary-shot-noise',
                        choices=[0, 1], type=int, required=True)
    parser.add_argument('--cov', choices=['mock', 'analytic'], required=True)
    parser.add_argument('--overwrite', action='store_true', default=True)
    parser.add_argument('--kmin', type=float, default=1e-4)
    parser.add_argument('--kmax', type=float, default=0.3)

    ns, args = parser.parse_known_args()

    add_commands(**vars(ns))
    RSDFitRunner.execute(args=args)
