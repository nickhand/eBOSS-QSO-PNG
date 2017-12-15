import argparse
from eboss_qso.fits import load_data_results, load_ezmock_results
import numpy
import collections
import emcee
from pyRSD.rsdfit.parameters import ParameterSet
from pyRSD.rsdfit.solvers.emcee_solver import ChainManager
from pyRSD.rsdfit.results import EmceeResults
from pyRSD.rsdfit.util import rsd_logging

def run_joint_mcmc_fit(kind, iterations, walkers, output, load_kwargs, joint_params=['f_nl']):
    """
    Run a joint mcmc fit to the NGC ans SGC sky regions.

    Parameters
    ----------
    kind : 'data', 'ezmock'
        the kind of data to load
    iterations : int
        the number of iterations to run
    walkers : int
        the number of walkers to use
    output : str
        the name of the output file to save
    load_kwargs : dict
        the dictionary of keywords to load the previous result for each sample
    joint_params : list of str
        the names of parameters that are only use a single value for each sky region.
    """
    # add a console logger
    rsd_logging.add_console_logger(0)

    # determine the loader
    assert kind in ['data', 'ezmock']
    if kind == 'data':
        loader = load_data_results
    elif kind == 'ezmock':
        loader = load_ezmock_results

    # load the previous results
    samples = ['N', 'S']
    drivers = {}
    for sample in samples:
        load_kwargs['sample'] = sample
        drivers[sample] = loader(**load_kwargs)
        drivers[sample].set_fit_results()

    # initialize a parameter set to hold combined data set
    theory = ParameterSet()

    # the initial state
    p0 = []
    free_names = [] # track names

    # add params that are fit jointly
    for param in joint_params:
        theory[param] = drivers['N'].theory.fit_params[param]
        p0.append(theory[param].value)
        free_names.append(param)

    # slices to set theta properly
    slices = collections.defaultdict(list)

    for tag in ['ngc', 'sgc']:
        tag_ = tag[0].upper()
        d = drivers[tag_] # the driver

        # add all free names
        for i, name in enumerate(d.theory.free_names):

            # only fit one version of joint params
            if name in joint_params:
                slices[tag_].append(joint_params.index(name))
                continue

            # add tagged version of other params
            slices[tag_].append(len(free_names))
            free_names.append(name+'_'+tag)

            par = d.theory.fit_params[name]
            p0.append(par.value)
            theory[name+'_'+tag] = par


    # the initial state
    p0 = numpy.array(p0)

    # the arrays of free parameters
    theta1 = numpy.array([drivers['N'].theory.fit_params[name].value for name in drivers['N'].theory.free_names])
    theta2 = numpy.array([drivers['S'].theory.fit_params[name].value for name in drivers['S'].theory.free_names])

    def objective(x):

        # get separate values
        theta1[:] = x[slices['N']]
        theta2[:] = x[slices['S']]

        # return
        return drivers['N'].lnprob(theta1) + drivers['S'].lnprob(theta2)


    # initialize the sampler
    ndim = len(p0)
    p0 = numpy.array([p0 + 1e-3*numpy.random.randn(ndim) for i in range(walkers)])
    sampler = emcee.EnsembleSampler(walkers, ndim, objective)

    # iterator interface allows us to tap ctrl+c and know where we are
    niters = iterations
    burnin = 0

    # do the sampling
    with ChainManager(sampler, niters, walkers, free_names, None) as manager:
        for niter, result in manager.sample(p0, None):

            # check if we need to exit due to exception/convergence
            manager.check_status()

            # update progress and test convergence
            manager.update_progress(niter)

    # re-order the chain to the sorted parameter order
    order = []
    for name in theory.free_names:
        order.append(free_names.index(name))
    sampler.chain = sampler.chain[...,order]

    # create and save the result
    result = EmceeResults(sampler, theory, burnin)
    print(result)
    result.to_npz(output)
