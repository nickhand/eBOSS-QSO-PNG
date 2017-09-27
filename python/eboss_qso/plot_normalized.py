from pyRSD import pygcl
from nbodykit import files, dataset
from matplotlib import pyplot as plt

z = 1.5
cosmo = pygcl.Cosmology('boss_dr12_fidcosmo.ini', pygcl.Cosmology.EH_NoWiggle)
Plin = pygcl.LinearPS(cosmo, z)

b1 = 2.1
f = cosmo.f_z(z)
P0_kaiser = lambda k: (1. + 2./3*f/b1 + 1./5*(f/b1)**2) * b1**2 * Plin(k)
P2_kaiser = lambda k: (4./3*f/b1 + 4./7*(f/b1)**2) * b1**2 * Plin(k)
P4_kaiser = lambda k: (8./35*(f/b1)**2) * b1**2 * Plin(k) 

d = dataset.Power1dDataSet.from_nbkit(*files.Read1DPlainText('./poles_eboss_v1.6-QSO-N-eboss_v1.6_zrange_0.8_2.2.dat'))
idx = d['k'] < 0.5
k = d['k'][idx]

plt.semilogx(k, k*(d['power_0.real'][idx] - d.attrs['shot_noise']), label=r"$P_0(k)$")
plt.semilogx(k, k*d['power_2.real'][idx], label=r"$P_2(k)$")
plt.semilogx(k, k*d['power_4.real'][idx], label=r"$P_4(k)$")

plt.semilogx(k, k*P0_kaiser(k), c='k')
plt.semilogx(k, k*P2_kaiser(k), c='k')
plt.semilogx(k, k*P4_kaiser(k), c='k')

ax = plt.gca()
ax.set_xlabel(r"$k$ $[h/\mathrm{Mpc}]$", fontsize=16)
ax.set_ylabel(r"$k P_\ell(k)$ $[\mathrm{Mpc}/h]^2$", fontsize=16)
ax.legend()

plt.savefig("normalize_qso_poles.pdf")
plt.show()
