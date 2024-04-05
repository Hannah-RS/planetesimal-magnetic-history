# shell script to synchronise code folder with hpc
cd #change to home directory
rsync -va --exclude=*/ Documents/Code/learning-model/ arc-file://data/eart-astroid-evo/exet5460/Asteroid_model #copy everything except directories
rsync -var Documents/Code/learning-model/Templates arc-file://data/eart-astroid-evo/exet5460/Asteroid_model #copy the templates directory 
