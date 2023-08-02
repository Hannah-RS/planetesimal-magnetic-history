#! /bin/bash
#SBATCH --job-name=arrayJob
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-560
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00 
#SBATCH --partition=short 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=hannah.sanderson@earth.ox.ac.uk

module load Anaconda3
source activate $DATA/viscosity
#go into model directory
cd Asteroid_model
#make directory
mkdir Results/params_$SLURM_ARRAY_TASK_ID
#copy autoparameters file for this run
cp Run_params/auto_params_$SLURM_ARRAY_TASK_ID.csv Results/params_$SLURM_ARRAY_TASK_ID/auto_params.csv
#run the model
bash multi_run.sh Results/params_$SLURM_ARRAY_TASK_ID
