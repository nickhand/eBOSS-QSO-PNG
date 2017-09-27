from nbodykit.lab import *
from nbodykit.dataset import DataSet

import argparse
import functools
from scipy.special import legendre
from scipy.integrate import simps
from scipy.optimize import leastsq
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from matplotlib import pyplot as plt

# some global configuration variables
NMU      = 100
ZMIN     = 0.8
ZMAX     = 2.2
NZ       = 100
FKP_P0   = 6e3
FSKY_NGC = 1001.2791 * (numpy.pi/180.)**2 / (4*numpy.pi)
FSKY_SGC = 821.82659 * (numpy.pi/180.)**2 / (4*numpy.pi)

def residuals(p0, x, y, ells, *args, **kws):
    """
    Return the chi-squared residual: ``(y-model)/error``, which
    should be used in :func:`scipy.optimize.leastsq`
    
    Parameters
    ----------
    p0 : array_like
        the array of model parameters to fit
    x : array_like
        the independent variable (wavenumber in this case)
    y : array_like
        the dependent variable
    ells : list of int
        the list of multipole numbers in the full data vector
    *args : tuple of arguments
        the additional arguments to pass to :func:`model_function`
    **kws : key/value pairs
        the additional keywords to pass to :func:`model_function`
    """
    theory = []; errs = []
    transfer = kws.pop('transfer', None)
    for i, ell in enumerate(ells):
        
        if transfer is not None:
            kws['transfer'] = transfer[i]
        else:
            kws['transfer'] = None
        t, e = multipoles_theory(ell, x, p0[0], p0[1], *args, **kws)
        theory.append(t)
        errs.append(e**0.5)
        
    theory = (numpy.asarray(theory).T).ravel(order='F')
    errs = (numpy.asarray(errs).T).ravel(order='F')
    
    return (y-theory)/errs
    

def multipoles_theory(ell, k, b1, sigma, cosmo, nbar, Plin, fsky, f=1.0, ellprime=None, transfer=None):
    """
    Compute the theoretical multipoles and the covariance
    
    Parameters
    ----------
    ell : int
        the multipole number to compute the results for
    cosmo : astropy.cosmology.FlatLambdaCDM
        the astropy cosmology instance
    nbar : callable
        function that returns number density as a function of redshift
    Plin : callanle
        callable that returns the linear power spectrum as a function of wavenumber
    k : array_like
        the independent variable, wavenumber
    b1 : float
        the linear bias parameter of the model
    sigma : float
        the velocity dispersion parameter of the model
    f : float
        the growth rate value to use in the model
    ellprime : int; optional
        if not None, return the covariance between ``ell`` and ``ellprime``
    transfer : callable; optional
        a window function transfer function to apply to the theory
    
    Returns
    -------
    theory : array_like
        the theoretical multipole power spectrum
    cov : array_like
        the theoretical Gaussian covariance between ``ell`` and ``ellprime``
    """
    if ellprime is None:
        ellprime = ell
        
    # the model P(k,mu) function
    Pk = Plin(k)[:,None]
    beta = f/b1
    pkmu_model = lambda k, mu: numpy.exp(-(k*mu*sigma)**2) * b1**2 * Pk * (1 + 2*beta*mu**2 + beta**2 * mu**4)
    mus = numpy.linspace(0, 1, NMU+1)
    Pkmus = pkmu_model(k[:,None], mus[None,:])
    
    # apply a transfer to P(k,mu)?
    if transfer is not None:
        Pkmus *= transfer(k)[:,None]

    # volume of redshift shells for integral over z
    zbins = numpy.linspace(ZMIN, ZMAX, NZ+1)
    R_hi = cosmo.comoving_distance(zbins[1:]).value * cosmo.h
    R_lo = cosmo.comoving_distance(zbins[:-1]).value * cosmo.h
    dV = (4./3.)*numpy.pi*(R_hi**3 - R_lo**3) * fsky
    
    # number density
    zcen = 0.5*(zbins[:-1] + zbins[1:])
    nbar_ = nbar(zcen)
    w = 1. / (1 + nbar_*FKP_P0) # FKP weights
    
    # normalization for redshift integral
    W = ((nbar_*w)**2 * dV).sum()
    
    # k-shell volume 
    dk = numpy.diff(k).mean()
    Vk  = 4*numpy.pi*k**2*dk
    modes = 2 * Vk / (2*numpy.pi)**3
            
    # the theory P(k,mu)
    kern = (2*ell+1.)*legendre(ell)(mus)
    theory = numpy.array([simps(kern*d, x=mus) for d in Pkmus])
    
    # P(k,mu)^2 * L_ell(mu)*L_ellprime
    weights = kern * (2*ellprime+1.)*legendre(ellprime)(mus)
    power = (Pkmus[...,None] + 1./nbar_)**2
    tobin =  weights[...,None] * power[None,...]
    
    # do the sum over redshift
    Psq = ( (w*nbar_)**4 * dV * tobin).sum(axis=-1) / W**2
    Psq = Psq.mean(axis=-1)
    cov = (2*Psq/modes)
        
    return theory, numpy.squeeze(cov)

def find_bestfit(filenames, nbar_file, sample, kmin=0, kmax=0.4, window_transfer=None):
    """
    Find the best-fit model
            
    Parameters
    ----------
    filenames : list of str
        a list of file names to fit and plot
    nbar : str
        the file holding the n(z) measurement
    sample : 'NGC' or 'SGC'
        either the north galactic cap (NGC) or south galactic cap (SGC)
    kmin : float; optional
        the minimum k in the fitting range
    kmax : float; optional
        the maximum k in the fitting range
    window_transfer : str; optional
        the file holding the window function transfer 
    
    Returns
    -------
    list of tuples : 
        a list of tuples for each input filename; each tuple gives
        the DataSet holding the data (with error columns) and 
        a best-fit theory function taking ``k``, ``ell`` and ``ellprime`` 
        which returns the best-fit theory multipole and the covariance
        between ``ell`` and ``ellprime``
    """    
    sample = sample.lower()
    if sample not in ['ngc', 'sgc']:
        raise ValueError("'sample' should be one of 'NGC' or 'SGC'")
    
    # set fsky for this sample
    if sample == 'ngc':
        fsky = FSKY_NGC
    else:
        fsky = FSKY_SGC
        
    # window function transfer
    T = [None]*3
    if window_transfer is not None:
        W = numpy.loadtxt(window_transfer)
        T = [spline(W[:,0], W[:,i]) for i in range(1, W.shape[1])]    
        
    # set up the cosmology
    z = 1.5
    cosmo = cosmology.Cosmology(H0=67.7, Om0=0.31, sigma8=0.8, n_s=0.97, flat=True)
    Plin = cosmology.EHPower(cosmo, z)
    Plin_nw = cosmology.NoWiggleEHPower(cosmo, z)
    growth_rate = cosmo.growth_rate(z)
    
    # load the n(z) from file
    nbar = numpy.loadtxt(nbar_file, skiprows=3)
    nbar = spline(nbar[:,0], nbar[:,3])
        
    toret = []
    for ifile, filename in enumerate(filenames):
        
        # load results in new format
        if filename.endswith('.json'):
            r = ConvolvedFFTPower.load(filename)
            d = r.poles
            Pshot = r.attrs['shotnoise']
    
        # old plaintext results
        else:
            d = DataSet.from_plaintext(['k'], filename)
            Pshot = d.attrs['shot_noise']
            
        idx = (d['k'] >= kmin)&(d['k'] <= kmax)
        k = d['k'][idx]
        
        # subtract shot noise
        d['power_0'] -= Pshot
        
        # concatenate into one data vector
        P0 = d['power_0'][idx].real
        P2 = d['power_2'][idx].real
        P4 = d['power_4'][idx].real
        poles = numpy.hstack([P0,P2,P4])
    
        # do the chi-squared fitting
        p0 = [2.0, 5.0]
        f = lambda p0: residuals(p0, k, poles, [0, 2, 4], cosmo, nbar, Plin, fsky, f=growth_rate, transfer=T)
        coeffs, matcov = leastsq(f, p0)
        print("b1 = %.6f" %coeffs[0])
        print("sigma = %.6f" %coeffs[1])
        
        # compute the errors
        for ell in [0, 2, 4]:
            _, cov = multipoles_theory(ell, d['k'], coeffs[0], coeffs[1], cosmo, nbar, Plin, fsky, f=growth_rate, transfer=T[(ell+1)//2])
            d['error_%d' %ell] = cov**0.5
        
        def theory(b1, sigma, k, ell, ellprime=None):
            return multipoles_theory(ell, k, b1, sigma, cosmo, nbar, Plin, fsky, f=growth_rate, ellprime=ellprime, transfer=T[(ell+1)//2])
        f = functools.partial(theory, *coeffs)
        f.b1 = coeffs[0]
        f.sigma = coeffs[1]
        
        toret.append((d, f))
        
    return toret
    

    
def plot_bestfit(bestfits, normalize=False, ells=[0,2,4], labels=None, 
                output=None, ylims=None, legend_ncol=1):
    """
    Plot the best fit results
            
    Parameters
    ----------
    bestfits : list
        list of tuples holding data and theory
    normalize : bool; optional
        when plotting, normalize the multipoles
    ells : list of int; optional
        which multipoles to fit and plot
    labels : list of str; optional
        label each set of results when plotting
    output : str; optional
        save the plot to this file
    ylims : tuple; optional
        set the limits of the y axis
    """       
    # set up the cosmology
    z = 1.5
    cosmo = cosmology.Cosmology(H0=67.7, Om0=0.31, sigma8=0.8, n_s=0.97, flat=True)
    Plin = cosmology.EHPower(cosmo, z)
    Plin_nw = cosmology.NoWiggleEHPower(cosmo, z)
    growth_rate = cosmo.growth_rate(z)
         
    for i, result in enumerate(bestfits):
            
        data, theory = result
        k = data['k']
        
        # no-wiggle monopole for normalization
        norm = lambda k, b1, f: (1. + 2./3*f/b1 + 1./5*(f/b1)**2) * b1**2 * Plin_nw(k)
    
        for ell in [0,2,4]:
              
            if ell not in ells:
                continue
                
            # best-fit theory
            y = theory(k, ell)[0]
        
            # the label
            label = r"$P_%d$" %ell
            if labels is not None:
                label = labels[i] + ": %s" %label
        
            P = data['power_%d' %ell]
            yerr = data['error_%d' %ell]

            # normalize by no-wiggle
            if not normalize:
                plt.errorbar(k, k*P, k*yerr, ls='', marker='.', label=label)
                plt.semilogx(k, k*y, c='k')
            else:
                N = norm(k, theory.b1, growth_rate)
                plt.errorbar(k, P / N, yerr / N, ls='', marker='.', label=label)
                plt.plot(k, y / N, c='k')
    
    # format the axes
    ax = plt.gca()
    if len(ells) > 1 or len(bestfits) > 1:
        ax.legend(loc=0, fontsize=14, ncol=legend_ncol)
        
    ax.set_xlabel(r"$k$ [$h$/Mpc]", fontsize=16)
    if not normalize:
        ax.set_ylabel(r"$k \ P_\ell$ $[\mathrm{Mpc}/h]^2$", fontsize=16)
    else:
        ax.set_ylabel(r"$P_\ell / P_0^\mathrm{EH}$", fontsize=16)
    
    if ylims is not None:
        ax.set_ylim(ylims)
    
    if output is not None:
        plt.savefig(output)
    
    return ax

if __name__ == '__main__':
    
    desc = 'fit a simple theoretical model to the quasar multipoles'
    parser = argparse.ArgumentParser(description=desc)
    
    h = 'the name of the power spectrum file from nbodykit'
    parser.add_argument('filenames', type=str, nargs='+', help=h)
    
    h = 'the name of the file holding the n(z) for the quasars'
    parser.add_argument('nbar_file', type=str, help=h)
    
    h = 'the minimum k in the fitting range'
    parser.add_argument('--kmin', default=0, type=float, help=h)
    
    h = 'the maximum k in the fitting range'
    parser.add_argument('--kmax', default=0.2, type=float, help=h)
    
    h = 'whether to normalize the plot'
    parser.add_argument('--normalize', action='store_true', help=h)
    
    h = 'which ells to plot'
    parser.add_argument('--ells', nargs='*', default=[0,2,4], choices=[0,2,4], type=int, help=h)
    
    h = 'labels to use for more than one data file'
    parser.add_argument('--labels', nargs='*', default=None, help=h)
    
    h = 'the output file name to save to'
    parser.add_argument('-o', '--output', type=str, help=h)
    
    h = 'the y limits to apply'
    parser.add_argument('--ylims', nargs=2, type=float, help=h)
    
    ns = parser.parse_args()
    ax = run(**vars(ns))
    plt.show()
    
    
    
    
