i=0
nruns=$(awk '(NR>2)' auto_params.csv | awk -F',' '$11<0' | wc -l) #only count unrun lines

while [ $i -lt $nruns ]
do
    runi=$(awk '(NR>2)' auto_params.csv | awk -F',' '$11<0' | awk -F',' 'NR==1{print $1}') #get run number
    SECONDS=0 #start timer
    python solver.py >> output.txt #run model and write terminal output to file
    runt=$SECONDS #stop timer
    printf '%s\n' $runi $runt | paste -sd ',' >> runtime.csv #save run number and time to file
    echo Run $runi has been completed in 
    printf '%dd:%dh:%dm:%ds\n' $((runt/86400)) $((runt%86400/3600)) $((runt%3600/60)) $((runt%60))
    i=$((i+1)) 
done

