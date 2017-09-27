from pylab import *

import argparse
import seaborn as sns
from pyRSD import pygcl
from pyRSD.rsd.window import WindowConvolution, convolve_multipoles


def main(ns):
    
    W = np.loadtxt(ns.filename)
    cosmo = pygcl.Cosmology('runPB.ini', pygcl.Cosmology.EH_NoWiggle)
    Plin = pygcl.LinearPS(cosmo, 0.)
    kaiser = pygcl.Kaiser(Plin, cosmo.f_z(2.5), 3.5)

    k = np.logspace(np.log10(ns.kmin), np.log10(ns.kmax), 4096)
    Pell = np.vstack([kaiser.P_ell(ell, k) for ell in ns.ells]).T

    convolver = WindowConvolution(W[:,0], W[:,1:], max_ellprime=4)
    kk, Pell_conv = convolve_multipoles(k, Pell, [0, 2, 4], convolver, dry_run=False, qbias=0., legacy=False)

    # save the transfer function
    out = [kk] + [Pell_conv[:,i]/Pell[:,i] for i in range(len(ns.ells))]
    out = np.vstack(out).T
    np.savetxt(ns.output, out)
    
    if ns.plot:
        c = sns.color_palette('Paired')
        for i in range(3):
            ell = 2*i
            semilogx(k, Pell[:,i], color=c[2*i], label=r'input $P_%d$' %ell)
            semilogx(kk, Pell_conv[:,i], color=c[2*i+1], label=r'conv. $P_%d$' %ell)
    
        ax = gca()
        ax.set_xlim(1e-3, 5.0)
        xlabel(r"$k$ $[h/\mathrm{Mpc}]$", fontsize=16)
        ylabel(r"$P_\ell$ $(\mathrm{Mpc}/h)^3$", fontsize=16)
        legend(loc=0, ncol=1, fontsize=12)
        show()
    
if __name__ == '__main__':
    
    desc = "compute the window transfer function, optionally plotting as well"
    parser = argparse.ArgumentParser(description=desc)
    
    h = 'the file holding the window multipoles'
    parser.add_argument('filename', type=str, help=h)
    
    h = 'the output transfer function file name'
    parser.add_argument('-o', '--output', type=str, required=True, help=h)
    
    h = 'whether to plot'
    parser.add_argument('--plot', action='store_true', help=h)
    
    h = 'which ells to output'
    parser.add_argument('--ells', type=int, nargs='*', default=[0,2,4], help=h)
    
    h = 'the minimum k value'
    parser.add_argument('--kmin', type=float, default=1e-4, help=h)
    
    h = 'the maximum k value'
    parser.add_argument('--kmax', type=float, default=1e4, help=h)
    
    ns = parser.parse_args()
    main(ns)