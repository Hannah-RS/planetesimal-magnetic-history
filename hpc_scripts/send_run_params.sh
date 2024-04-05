#send a folder of run parameters to the hpc
cd
scp -r Documents/Code/learning-model/Run_params/$1/ arc-file://data/eart-astroid-evo/exet5460/Asteroid_model/Run_params
