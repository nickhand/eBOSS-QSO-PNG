#!/bin/bash

rsync -e ssh -avzl $@ --progress nhand@cori:/global/cscratch1/sd/nhand/Research/eBOSS/measurements/ .
