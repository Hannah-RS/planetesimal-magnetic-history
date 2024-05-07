# learning-model
Thermal evolution and dynamo generation model for a planetesimal. The model is described in Sanderson et. al. 2024a, "Unlocking planetesimal magnetic field histories: a refined, versatile model for thermal evolution and dynamo generation", which will be submitted to Icarus and available as a preprint soon. 


## Basic model description

## Contents of this repository

## Downloading the code
### Dependencies


## How to run the model
Describe required file structure and set up (run parameters and results folder) 
### Directory structure
```mermaid
    ---
    title: Directory structure
    ---
    flowchart TD
    Repository --- Run_parameters & Results
    Run_parameters --- a["Single_run"] & b["Multi_run"]
    b --- params_1 & params_2
    Results --- c["Single_run"] & d["Multi_run"]
    d --- e["params_1] & f"[params_2] & all_sucess_info.csv & fail_params.csv & inval_params.csv

```

### Single run
1. Set `automated=False` in `parameters.py` and then change your desired parameters below. 
2. Run `solver.py`. Results will appear in your chosen results directory

### Multiple runs one folder
1. Create the directory where you wish to have the results.
2. Create auto_params.csv with run parameters in the desired folder.
3. Run `bash multi_run.sh <filepath>` in the main directory where filepath is the path from the model directory to where you want to save results

### Multiple runs - multiple subfolders
This method is for keeping a set of parameters constant, whilst varying an individuals parameter in independent runs. It requires access to an hpc in order to create an array job. For example, the directory "300km_runs" could contain subfolders "params_1" nd "params_2" each containing an "auto_params.csv" file and could run as one directory using an array job on an hpc. Both subdirectories could be for a 300km radius planetesimal: "params_1/" could vary reference viscosity but keep all other variables constant and "params_2" could vary critical melt fraction and keep all other variables constant. This method was used in Sanderson et. al. 2024b to investigate the role of viscosity, core sulfur content and radiogenic $^{60}Fe$ on dynamo generation.
 
1. Create overall parameters directory
2. Use `create_single_csv.py` to set your range of parameter space you are exploring and create the parameter subfolders.
3. Send parameters to remote machine
4. Submit an array job to run `multi_run.sh` on each of your parameter subfolders. 
5. Copy all results off the hpc
6. Copy `all_sucess_info.csv`, `fail_params.csv` and `inval_params.csv` from Templates to results folder 
7. Use `merge_results.py <results_folder> <number_of_subfolders>` to merge all results from a group of runs into a combined data file 
8. Now enjoy analysing your data!

## Plotting the output
Npz files, example plotting notebook

## Citing the model

## Issues

## References




