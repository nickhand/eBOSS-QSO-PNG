# eBOSS QSO EZ-mocks: v1.8

The notes here are taken from:

https://trac.sdss.org/wiki/eBOSS/QGC/LSScatsv1.9f

There are several sub-versions of v1.8 of the EZ-mocks. All of these
sub-versions agree in the same footprint and only differ on how the
redshift failures and collision objects are treated.

The description of the different v1.8 sub-versions of the EZ-mocks are:

1. **v1.8e-reg**: Usual wtot = wcp + wnoz -1. This corresponds to cp + zfw
cases, if you have access to these pre-versions. This should correspond as
 well to "v1.8e" in the old notation of the EZ mocks.

2. **v1.8e-no**: The full objects are considered. No close pairs (CP) or
redshift failures (RF) are applied. This is the reference value to
test systematics. This should correspond to the "v1.8" of the EZ in the
old notation.

3. **v1.8e-fph**: cp + zf using Focal Plane Eff (template from Hector)


## NERSC Locations

- **v1.8e-reg**: /global/project/projectdirs/eboss/qso_mocks_v1.8_wnoz_wcp/zevoEZmock_QSO_v1.8e_veto_ngc.tar.gz
- **v1.8e-no**:  /global/project/projectdirs/eboss/QSO_v1.8_EZmock_fkp_corrected/zevoEZmock_?gc_nzfixed.tar.gz
- **v1.8e-fph**:  /global/project/projectdirs/eboss/EZmocks_v1.8weightfocal/
