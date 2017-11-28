#!/bin/bash

#SBATCH -p {{ partition }}
#SBATCH -J {{ job }}.{{ cores }}
#SBATCH -o {{ output_file }}
#SBATCH -N {{ nodes }}
#SBATCH -t {{ time }}
{{ haswell_config }}

cd $PROJECT_HOME/Research/eBOSS/python

# activate environment
source /usr/common/contrib/bccp/conda-activate.sh 3.6

# install correct nbodykit version to computing nodes
bcast-pip .
bcast $LOCAL_STACK/anaconda3/envs/pyrsd-anaconda-3.6-$NERSC_HOST.tar.gz

cd $SLURM_SUBMIT_DIR

# call run-tests with desired number of cores
echo ===== Running with {{ cores }} cores =====
srun -n {{ cores }} {{ command }}
