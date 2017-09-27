#! /usr/bin/env python
from lsskit.rsdfit import RSDFitRunner

SAMPLES = ['N', 'S']
ZBINS = [(0.9, 1.2), (1.2, 1.5), (1.5, 1.8), (1.8, 2.2), (0.9, 2.2)]

for sample in SAMPLES:
    for zbin in ZBINS:
        zmin, zmax = zbin

        args = (sample, zmin, zmax)
        param_file = "$THESIS_DIR/eBOSS/v1.8/Results/Params/params_v1.8-QSO-%s_%.1f-%.1f.dat" %args
        output = "$THESIS_DIR/eBOSS/v1.8/Results/Fits/QSO-%s-%.1f-%.1f" %args
        command = "rsdfit mcmc -p %s -o %s -i 1000 -w 10" %(param_file, output)
        RSDFitRunner.register(command)


if __name__ == '__main__':

    RSDFitRunner.execute()
