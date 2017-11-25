from ebossqso import config
from pyRSD.rsd import QuasarSpectrum
from pyRSD.rsdfit.data import PoleCovarianceMatrix

import argparse
import numpy
from scipy.special import legendre

def _compute_covariance(config, zmin, zmax, k, ells, model, P0_FKP=0., Nmu=100, Nz=50):
    """
    The internal function to compute the covariance
    """
    # get config
    cosmo = config.cosmo
    fsky = config.fsky
    nbar = config.nbar

    # the best-fit P(k,mu)
    mus = numpy.linspace(0, 1, Nmu+1)
    Pkmu = model.power(k, mus).values
    _, mus = numpy.meshgrid(k, mus, indexing='ij')

    N1 = Pkmu.shape[0]
    N2 = len(ells)

    # volume of redshift shells for integral over z
    zbins = numpy.linspace(zmin, zmax, Nz+1)
    R_hi = cosmo.comoving_distance(zbins[1:]) * cosmo.h # in Mpc/h
    R_lo = cosmo.comoving_distance(zbins[:-1]) * cosmo.h # in Mpc/h
    dV = (4./3.)*numpy.pi*(R_hi**3 - R_lo**3)

    # compute nbar
    zcen = 0.5*(zbins[:-1] + zbins[1:])
    nbar_ = nbar(zcen)

    # weights
    w = 1. / (1 + nbar_*P0_FKP) # FKP weights

    # properly calibrate fsky
    W = ((nbar_*w)**2 * dV).sum()
    dV *= fsky
    W *= fsky

    # k-shell volume
    dk = numpy.diff(k).mean()
    Vk  = 4*numpy.pi*k**2*dk
    modes = 2 * Vk / (2*numpy.pi)**3

    # initialize the return array
    Psq = numpy.zeros((N2, N1)*2)
    modes = numpy.zeros((N2, N1)*2)

    # leg weights
    leg = numpy.array([(2*ell+1)*legendre(ell)(mus) for ell in ells])

    # P(k,mu)^2 * L_ell(mu)*L_ellprime
    weights = leg[:,None]*leg[None,:]
    power = (Pkmu[...,None] + 1./nbar_)**2
    tobin = weights[...,None] * power[None,...]

    # fill the covariance for each ell, ell_prime
    for i in range(N2):
        for j in range(i, N2):

            # do the sum over redshift first
            x = ( (w*nbar_)**4 * dV * tobin).sum(axis=-1) / W**2
            Psq[i,:,j,:] = numpy.diag(numpy.nanmean(x[i,j,:], axis=-1))
            modes[i,:,j,:] = 2 * Vk / (2*numpy.pi)**3
            if i != j:
                Psq[j,:,i,:] = Psq[i,:,j,:]
                modes[j,:,i,:] = modes[i,:,j,:]

    # reshape squared power and modes
    Psq = Psq.reshape((N1*N2,)*2)
    modes = modes.reshape((N1*N2,)*2)
    C = numpy.nan_to_num(2*Psq/modes)

    # the coordinate arrays
    k_coord = numpy.concatenate([k for i in range(len(ells))])
    ell_coord = numpy.concatenate([numpy.ones(len(k), dtype=int)*ell for ell in ells])

    return PoleCovarianceMatrix(C, k_coord, ell_coord, verify=False)


def compute_covariance(config, stats, z_eff, b1, sigma_fog, zmin, zmax, P0_FKP,
                        output, kmin=0., kmax=0.7, dk=0.005, quiet=False, rescale=1.0):
    """
    Compute the multipole covariance matrix

    Parameters
    ----------
    config : eBOSSConfig
        the configuration object
    stats : list of str
        the list of statistics to compute
    z_eff : float
        the effective redshift of the sample to evaluate the model P(k,mu)
    b1 : float
        the linear bias of the model
    sigma_fog : float
        the FOG velocity dispersion of the model
    zmin : float
        the minimum redshift of this sample
    zmax : float
        the maximum redshift of this sample
    P0_FKP : float, optional
        the FKP P0 value to use
    output : str
        the name of the file to write the output to
    kmin : float; optional
        the minimum k value of the grid used to evaluate the covariance matrix
    kmax : float; optional
        the maximum k value of the grid used to evaluate the covariance matrix
    dk : float; optional
        the k grid spacing
    lmax : int; optional
        the maximux multipole number to include
    """
    # load the model
    model = QuasarSpectrum(z=z_eff, params=config.cosmo)
    model.b1 = b1
    model.sigma_fog = sigma_fog

    if 'P0_sysfree' in stats:
        ells = [0, 2]
    else:
        ells = [int(stat[-1]) for stat in stats]

    # get C
    ks = numpy.arange(kmin, kmax, dk) + 0.5*dk
    C = _compute_covariance(config, zmin, zmax, ks, ells, model, P0_FKP=P0_FKP)

    # covariance for P0 + 2./5*P2
    if 'P0_sysfree' in stats:
        C00 = C.sel(ell1=0, ell2=0)
        C22 = C.sel(ell1=2, ell2=2)
        C02 = C.sel(ell1=0, ell2=2)
        C = C00 + 2./5 * C22 + 2 * 2./5 * C02

    # rescale
    C *= rescale

    # and save
    if not quiet:
        print("saving %s..." %output)
    C.to_plaintext(output)
