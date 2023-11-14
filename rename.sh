# script for renaming files - don't use for renaming files upwards as you overwrite subsequent files!
#$1 directory
#$2 number to subtract from run number

#change to directory
cd Results_combined/$1

#get value of number in file 
i = 0
declare -a new_number=()
for file in run_*_diff.npz; do
    number=$(echo $file | tr -dc '0-9') #find current number
    if [ $number -gt 0 ] #don't try rename non npz files
    then
        new_number+=$(($number-$2)) #append to new number array
        #echo ${new_number[i]}
        i=$(($i+1))
        echo $new_number
       # mv run_${number}_diff.npz run_${new_number}_diff.npz
        #mv run_${number}.npz run_${new_number}.npz
    fi
done
