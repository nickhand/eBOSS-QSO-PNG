from pylab import *

import seaborn as sns
import argparse

def main(ns):
    
    
    W = np.loadtxt(ns.filename)
    for i in range(1, W.shape[1]):
        ell = 2*(i-1)
        semilogx(W[:,0], W[:,i], label=r'$Q_{%d}$' %ell)
        
    ax = gca()
    ax.set_xlim(left=0.1)
    xlabel(r"$s$ $[\mathrm{Mpc}/h]$", fontsize=16)
    ylabel(r"$Q_\ell$", fontsize=16)
    legend(loc=0, ncol=3, fontsize=16)
    if ns.output is not None:
        savefig(ns.output)
    show()

if __name__ == '__main__':
    
    desc = "plot the window function multipoles"
    parser = argparse.ArgumentParser(description=desc)
    
    h = 'the file holding the window multipoles'
    parser.add_argument('filename', type=str, help=h)
    
    h = 'the output file name'
    parser.add_argument('-o', '--output', type=str, help=h)

    ns = parser.parse_args()
    main(ns)