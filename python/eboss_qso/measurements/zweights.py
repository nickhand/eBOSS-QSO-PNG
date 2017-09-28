def bias_weight(z, cosmo):
    """
    Compute the bias redshift weight.

    .. math::

        W = D(z) [ b(z) + f(z)/3 ]
    """
    b = bias_model(z)
    D = z.map_blocks(cosmo.scale_independent_growth_factor)
    f = z.map_blocks(cosmo.scale_independent_growth_rate)

    return D*(b+f/3.)


def fnl_weight(z, p=1.6):
    """
    Compute the fnl redshift weight.

    .. math::

        W = (b(z) - p),

    with p = 1.0 or 1.6.
    """
    b = bias_model(z)
    return b - p

def bias_model(z):
    """
    Return the bias fit as a function of redshift from Laurent et al. 2017
    (1705.04718).
    """
    alpha = 0.278
    beta = 2.393
    return alpha * ( (1+z)**2 - 6.565 ) + beta
