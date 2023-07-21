#send a folder of run parameters to the hpc
scp -r Documents/Code/learning-model/Run_params/$1/ exet5460@gateway.arc.ox.ac.uk://data/eart-astroid-evo/exet5460/Asteroid_model/Run_params
