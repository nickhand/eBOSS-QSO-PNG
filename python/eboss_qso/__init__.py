import os
import socket

# data and result directories
if socket.gethostname() == 'utley':
    # define some directories
    EBOSS_DIR = os.path.join(os.environ['THESIS_DIR'], 'eBOSS-QSO-PNG')
    EBOSS_FITS = os.path.join(EBOSS_DIR, 'fits')
    EBOSS_SPECTRA = os.path.join(EBOSS_DIR, 'measurements', 'spectra')
else:
    EBOSS_DIR = None
    EBOSS_FITS = None
    EBOSS_SPECTRA = None
