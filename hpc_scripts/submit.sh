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
mkdir Results/$1
#copy autoparameters file for this run
cp Run_params/$1/auto_params.csv Results/$1/auto_params.csv
#run the model
bash multi_run.sh Results/$1

