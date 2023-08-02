## Move results from a set of runs to their correct folders ready for extraction
cd $DATA/Asteroid_model/Results
mkdir $1 #make directory
mv params_* $1/ #move results to directory
cd ../Run_params
mv auto_params* $1/ #move run params back to their folder
cd ../.. 
mkdir Output_$1 #move output files to a folder
mv arrayJob* Output_$1/  
