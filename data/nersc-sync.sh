#!/bin/bash

rsync -e ssh -avzl $@ --progress --delete . nhand@cori:/global/cscratch1/sd/nhand/eBOSS/data/
