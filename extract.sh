#script to copy results off the hpc and merge them together
#give the folder name on the command line
#run the command from learning-model folder 
mkdir Results_combined/$1 #make directory for results
scp -r arc-file://data/eart-astroid-evo/exet5460/Asteroid_model/Results/$1/* Results_combined/$1/ #copy from hpc
nfolders=$(ls Results_combined/$1|wc -l) #count number of folders
cp all_sucess_info.csv Results_combined/$1/all_sucess_info.csv #copy a template for final sucess results
cp fail_params.csv Results_combined/$1/fail_params.csv #copy a template for final fail results

python merge_results.py $1/ $nfolders #merge results, give folder and number of subfolders as arguments

