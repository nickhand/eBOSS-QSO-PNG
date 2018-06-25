from eboss_qso.fits.joint import run_joint_mcmc_fit
from eboss_qso.measurements.utils import make_hash
import os.path as osp
import os
from glob import glob


ARGS = [(False, 1.0),
          (False, 1.6),
          (True, 1.6),
          (True, 1.0)
          ]
ITERATIONS = 500
WALKERS = 100

def main(argnum, kmin):

    z_weighted, p = ARGS[argnum]

    # the data to load
    kws = {}
    kws['version'] = 'v1.9f'
    kws['krange'] = '%s-0.3' % kmin
    kws['params'] = 'basemodel-N-fnl'
    kws['zrange'] = '0.8-2.2'
    kws['z_weighted'] = z_weighted
    kws['p'] = p

    hashstr = make_hash(kws)

    # output directory
    output = osp.join(os.environ['EBOSS_FITS'], 'data')
    output = osp.join(output, kws['version'], kws['krange'], kws['params'], kws['zrange'])
    output = osp.join(output, 'QSO-N+S-%s' %hashstr)

    if not osp.exists(output):
        os.makedirs(output)

    # output file name
    i = len(glob(osp.join(output, '*npz')))
    output = osp.join(output, 'chain_%dx%d_%d.npz' %(ITERATIONS,WALKERS,i))
    print(output)

    # run
    run_joint_mcmc_fit('data', ITERATIONS, WALKERS, output, kws, joint_params=['f_nl'])


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("argnum", type=int, choices=[0,1,2,3])
    parser.add_argument('kmin', type=str, choices=["0.0001", "0.005"])

    ns = parser.parse_args()
    main(ns.argnum, ns.kmin)

