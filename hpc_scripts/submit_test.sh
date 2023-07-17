#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00 
#SBATCH --partition=devel 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=hannah.sanderson@earth.ox.ac.uk
#SBATCH --job-name="test"

module load Anaconda3
source activate $DATA/viscosity

#add a test here first
multi_run_test.sh
