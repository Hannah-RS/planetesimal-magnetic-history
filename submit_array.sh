#! /bin/bash
#alter the SBATCH info here
#SBATCH --job-name=arrayJob
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00 


#activate your virtual environment
module load Anaconda3
source activate $DATA/viscosity
#go into model directory
cd Asteroid_model
#make directory for results
mkdir Results/$1
#make subdirectory
mkdir Results/$1/params_$SLURM_ARRAY_TASK_ID
#copy autoparameters file for this run set
cp Run_params/$1/auto_params_$SLURM_ARRAY_TASK_ID.csv Results/$1/params_$SLURM_ARRAY_TASK_ID/auto_params.csv
#run the model
bash multi_run.sh Results/$1/params_$SLURM_ARRAY_TASK_ID
