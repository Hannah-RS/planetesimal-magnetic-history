# run multiple parameter combinations
#if directory is empty except for parameters file copy csv templates from main directory
if [ $(ls $1|wc -l) -lt 1 ]; then
    echo "No auto parameters file - aborting"
else
    if [ $(ls $1|wc -l) -gt 1 ]; then
        echo "Directory already has files - appending"
    else
        echo "Directory is empty - copying templates"
        cp Templates/run_results.csv $1
    fi
    i=0
    nruns=$(awk '(NR>2)' $1/auto_params.csv | awk -F',' '$NF<-1' | wc -l) #only count unrun lines - status is in last column
    echo $nruns
    while [ $i -lt $nruns ]
    do
        runi=$(awk '(NR>2)' $1/auto_params.csv | awk -F',' '$NF<-1' | awk -F',' 'NR==1{print $1}') #get run number
        SECONDS=0 #start timer
        python solver.py $1/ >> $1/output.txt 2>&1 #run model and write terminal output to file, send error there too
        runt=$SECONDS #stop timer
        
        echo Run $runi has been completed in 
        printf '%dd:%dh:%dm:%ds\n' $((runt/86400)) $((runt%86400/3600)) $((runt%3600/60)) $((runt%60))
        i=$((i+1)) 
    done

fi


