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
bcast-pip git+git://github.com/bccp/nbodykit.git

# install correct nbodykit version to computing nodes
bcast-pip .
bcast $PROGRAMS_DIR/pkgs/$NERSC_HOST/pyRSD*

cd $SLURM_SUBMIT_DIR

# call run-tests with desired number of cores
echo ===== Running with {{ cores }} cores =====
srun -n {{ cores }} {{ command }}
