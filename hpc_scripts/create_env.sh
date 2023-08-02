!# /bin/bash

# Load the version of Anaconda you need
module load Anaconda3

# Create an environment in $DATA and give it an appropriate name
export CONPREFIX=$DATA/viscosity
conda create --prefix $CONPREFIX

# Activate your environment
source activate $CONPREFIX

# Install packages
conda install numpy
conda install pandas
conda install scipy
conda install matplotlib
