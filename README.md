# learning-model
Thermal evolution model of an asteroid improving on Dodds et. al. (2021), Bryson et. al. (2019) and Nimmo (2009).

## Workflow for hpc
1. Use `create_csv.py` to set your range of parameter space you are exploring and create a folder of parameters
2. Run ` bash hpc_scripts/send_run_params.sh` to send parameters to Run_params folder on the hpc
3. Login to the hpc and move to the $DATA folder
4. Edit `submit_array2.sh` so the array=1-nfiles where nfiles is the number of csvs for all your runs
5. Submit array job using `sbatch submit_array.sh Folder` where Folder is the name of the folder for your run parameters
6. Wait for job to run
7. Run `extract.sh dir_name` where `dir_name` is your new directory to copy all results off the hpc and merge locally
8. Now enjoy analysing your data!

## Multiple runs on local machine
1. Create the directory where you wish to have the results.
2. Create auto_params.csv with run parameters in the desired folder.
3. Run "bash multi_run.sh filepath" in the main directory where filepath is the path from the model directory to where you want to save results

## Single run
1. Set `automated=False` in `parameters.py` and then change your desired parameters below. 
2. Run `solver.py`. Results will appear in Results_combined
