from pyRSD.rsdfit.data import PoleCovarianceMatrix
from eboss_qso.measurements import bias_model

import numpy
from scipy.special import legendre
from scipy.integrate import trapz


def _compute_covariance(config, zmin, zmax, kmin, kmax, dk, ells, model,
                        P0_FKP=0., Nmu=100, Nz=50, Nk=10):
    """
    The internal function to compute the covariance
    """
    from astropy import units, constants

    def evolve(model, z):
        model.z = z
        model.f = model.cosmo.f_z(z)
        model.sigma8_z = model.cosmo.Sigma8_z(z)
        model.b1 = bias_model(z)

    # setup
    fsky = config.fsky
    nbar = config.nbar
    cosmo = model.cosmo

    # the finer k binning
    dk_ = dk/Nk
    kedges = numpy.arange(kmin, kmax+0.5*dk_, dk_)
    kcen = 0.5*(kedges[1:] + kedges[:-1])

    # the 2D (k,mu) grid
    mus_ = numpy.linspace(0, 1, Nmu+1)
    _, mus = numpy.meshgrid(kcen, mus_, indexing='ij')

    # the output k grid
    kedges = numpy.arange(kmin, kmax+0.5*dk, dk)
    kout = 0.5*(kedges[1:] + kedges[:-1])

    N1 = len(kout)
    N2 = len(ells)

    # setup redshift binning
    zbins = numpy.linspace(zmin, zmax, Nz+1)
    zcen = 0.5*(zbins[:-1] + zbins[1:])

    # the comoving volume element
    Da = cosmo.Da_z(zcen) * cosmo.h()  # in Mpc/h
    dV = 4*numpy.pi*(1+zcen)**2 * Da**2 / cosmo.H_z(zcen) * cosmo.h()
    dV *= constants.c.to(units.km/units.second).value
    dV *= numpy.diff(zbins)

    # compute nbar
    nbar_ = nbar(zcen)

    # weights
    w = 1. / (1 + nbar_*P0_FKP)  # FKP weights

    # properly calibrate fsky
    W = ((nbar_*w)**2 * dV).sum()
    dV *= fsky
    W *= fsky

    # k-shell volume
    Vk = 4*numpy.pi/3. * (kedges[1:]**3 - kedges[:-1]**3)

    # initialize the return array
    cov = numpy.zeros((N2, N1)*2)

    # leg weights
    leg = numpy.array([(2*ell+1)*legendre(ell)(mus) for ell in ells])

    # compute Pkmu, optionally evolving in redshift
    if evolve is not None:
        Pkmu = []
        for zi in zcen:
            evolve(model, zi)
            Pkmu.append(model.power(kcen, mus_).values)
        Pkmu = numpy.asarray(Pkmu)
        Pkmu = numpy.moveaxis(Pkmu, 0, -1)
    else:
        Pkmu = model.power(kcen, mus_).values
        Pkmu = Pkmu[..., None]

    # P(k,mu)^2 * L_ell * L_ellprime
    weights = leg[:, None]*leg[None, :]
    power = (Pkmu + 1./nbar_)**2
    tobin = weights[..., None] * power[None, ...]

    # do the sum over redshift first
    x = ((w*nbar_)**4 * dV * tobin).sum(axis=-1) / W**2

    # the normalization
    norm = 4*(2*numpy.pi)**4 / (Vk**2)

    # split k into chunks to average
    N_chunks = int(len(kcen)/Nk)
    kcen_split = numpy.split(kcen, N_chunks)

    # fill the covariance for each ell, ell_prime
    for i in range(N2):
        for j in range(i, N2):

            # do the mu integral (mean is okay b/c of [0,1] domain)
            t = numpy.nanmean(x[i, j, :], axis=-1)

            # do the k averaging
            x_split = numpy.split(t, N_chunks)
            t = numpy.array([trapz(xx * kk**2, x=kk)
                             for xx, kk in zip(x_split, kcen_split)])

            cov[i, :, j, :] = norm * numpy.diag(t)
            if i != j:
                cov[j, :, i, :] = cov[i, :, j, :]

    # reshape squared power and modes
    cov = cov.reshape((N1*N2,)*2)
    cov = numpy.nan_to_num(cov)

    # the coordinate arrays
    k_coord = numpy.concatenate([kout for i in range(len(ells))])
    ell_coord = numpy.concatenate(
        [numpy.ones(len(kout), dtype=int)*ell for ell in ells])

    return PoleCovarianceMatrix(cov, k_coord, ell_coord, verify=False)


def compute_analytic_covariance(config, stats, model, zmin, zmax, P0_FKP,
                                kmin=0., kmax=0.7, dk=0.005, rescale=1.0):
    """
    Compute the multipole covariance matrix

    Parameters
    ----------
    config : eBOSSConfig
        the configuration object
    stats : list of str
        the list of statistics to compute
    model : QuasarSpectrum
        the quasar model object
    zmin : float
        the minimum redshift of this sample
    zmax : float
        the maximum redshift of this sample
    P0_FKP : float
        the FKP P0 value to use
    kmin : float; optional
        the minimum k value of the grid used to evaluate the covariance matrix
    kmax : float; optional
        the maximum k value of the grid used to evaluate the covariance matrix
    dk : float; optional
        the k grid spacing
    """
    if 'P0_sysfree' in stats:
        ells = [0, 2]
    else:
        ells = [int(stat[-1]) for stat in stats]

    # get C
    C = _compute_covariance(config, zmin, zmax, kmin,
                            kmax, dk, ells, model, P0_FKP=P0_FKP)

    # covariance for P0 + 2./5*P2
    if 'P0_sysfree' in stats:
        C00 = C.sel(ell1=0, ell2=0)
        C22 = C.sel(ell1=2, ell2=2)
        C02 = C.sel(ell1=0, ell2=2)
        C = C00 + 2./5 * C22 + 2 * 2./5 * C02

    # rescale
    C *= rescale

    return C


def compute_ezmock_covariance(version, sample, ells, kmin=0., kmax=0.7, p=None):
    """
    Compute the covariance from the EZ mocks.
    """
    from eboss_qso.measurements.results import load_ezmock_spectra

    # load the mocks
    Pell = []
    for ell in ells:
        mocks = load_ezmock_spectra(version, sample, p=p, ell=ell,
                                    subtract_shot_noise=True, average=False)
        Pell.append(mocks['power_%d' % ell].real)

    Pell = numpy.concatenate(Pell, axis=-1)
    cov = numpy.cov(Pell, rowvar=False)

    k = mocks['k'].mean(axis=0)
    k_coord = numpy.concatenate([k for i in range(len(ells))])
    ell_coord = numpy.concatenate(
        [numpy.ones(len(k), dtype=int)*ell for ell in ells])

    # create and slice to correct range
    C = PoleCovarianceMatrix(cov, k_coord, ell_coord, verify=False)
    valid_range = slice(kmin, kmax)
    C = C.sel(k1=valid_range, k2=valid_range)

    # store the number of mocks
    C.attrs['Nmock'] = len(mocks)

    return C
