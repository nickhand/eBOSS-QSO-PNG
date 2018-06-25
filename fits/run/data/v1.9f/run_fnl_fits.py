#! /usr/bin/env python
from eboss_qso.fits import RSDFitRunner, parametrize
from eboss_qso import EBOSS_SPECTRA
import os

NWALKERS = 50
ITERATIONS = 500
ZBINS = [(0.8, 2.2)]
VERSION = 'v1.9f'

# hash values for different p
HASHES = [('bba5aabfa6', 1.0),
          ('bba5aabfa6', 1.6),
          ('8ba79d78df', 1.6),
          ('0ef5c14a14', 1.0)
          ]
effective_redshifts = {'bba5aabfa6': 1.557,
                       '0ef5c14a14': 1.728,
                       '8ba79d78df': 1.829}

@parametrize({'sample':['N', 'S'], 'hashstr':HASHES})
def add_commands(sample, hashstr, vary_shot_noise=True, cov='analytic',
                    use_temp_files=False, overwrite=False, kmin=1e-4, kmax=0.3):

    # unpack the tuple of hashstring
    hashstr, p = hashstr

    # determine the params we are fitting
    all_params = ['b1', 'sigma_fog', 'f_nl']
    if vary_shot_noise:
        all_params += ['N']
    params = " ".join(all_params)

    dirname = os.path.join(EBOSS_SPECTRA, 'data', VERSION)
    filename = os.path.join(dirname, f'poles_eboss_{VERSION}-focal-QSO-{sample}-{hashstr}.json')

    # the effective redshift
    zeff = effective_redshifts[hashstr]

    # make the command
    command = f"eboss-qso-fit mcmc -f {filename} --vary {params} --stats P0 P2"
    command += f" -i {ITERATIONS} -w {NWALKERS} --kmin {kmin} --kmax {kmax} --cov {cov} --zeff {zeff} -p {p}"
    if use_temp_files:
        command += " --use-temp-files"
    if overwrite:
        command += " --overwrite"

    # and register it
    weighted = hashstr != "bba5aabfa6"
    tag = {'sample':sample, 'p':p, 'weighted':weighted, 'zbin':ZBINS[0]}
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
