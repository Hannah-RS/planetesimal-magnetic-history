#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00 
#SBATCH --partition=short 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=hannah.sanderson@earth.ox.ac.uk
#SBATCH --job-name=test
#SBATCH --output=test.out
module load Anaconda3
source activate $DATA/viscosity
#go into model directory
cd Asteroid_model
#run the model
bash multi_run.sh Results/test2/
