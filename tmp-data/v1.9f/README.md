# eBOSS QSO: v1.9f

Files have been downloaded from:

https://data.sdss.org/sas/ebosswork/eboss/lss/catalogs/1.9f

See notes on v1.9f here:

https://trac.sdss.org/wiki/eBOSS/QGC/September062017
https://trac.sdss.org/wiki/eBOSS/QGC/LSScatsv1.9f

## Differences in weighting schemes

### BOSS-like Weights

For the files

- eboss_v1.9f-QSO-S-eboss_v1.9f.dat.fits
- eboss_v1.9f-QSO-N-eboss_v1.9f.dat.fits

use the standard treatment of redshift failures, where the
closest neighbor in the sector has been up-weighted by +1
(the same used for BOSS). In this case, each galaxy has to be up-weighted by:

wtot = wfkp x wsys x (wnoz + wcp - 1)


### New Focal Plane Weights

The files:

- eboss_v1.9f-QSO-N-eboss_v1.9f-focal.dat.fits
- eboss_v1.9f-QSO-S-eboss_v1.9f-focal.dat.fits

contain focal plane weights to account for the redshift failure
effects. We recommend to use this kind of weighting scheme to
minimize the impact of redshift failures. In this case, each galaxy
has to be up-weighted by

w_tot = wfkp x wsys x wfocal x wcp
