#!/usr/env/bin bash

# Default values for parameters

# Number of iterations for MDI to perform
num_iter=10000

# Thinning variable for MDI
thin=25

# Random seed for MDI
seed=1

# Directory to save output to
save_dir="./output/"

# Receive parameters from command line 
while getopts n:t:s:o: option; do
    case "${option}" in
        n) num_iter=${OPTARG};;
        t) thin=${OPTARG};;
        s) seed=${OPTARG};;
        o) save_dir=${OPTARG};;
    esac
done

# Find all csv files of uppercase letters and numbers
files=( *[[:upper:][:digit:]].csv )

# Find the number of files matching this criterion
max=${#files[@]}

# Loop over all files unique combinations performing
# MDI on each combination and saving to a unique file
for ((id_a=0; id_a<max; id_a++)); do
    for ((id_b = id_a+1; id_b<max; id_b++)); do

        # The tissue types (file stripped of extension)
        file_a="${files[id_a]%.*}"
        file_b="${files[id_b]%.*}"

        # The output file name
        output="${save_dir}out_${file_a}_${file_b}.csv"

        # Call MDI
        ~/Desktop/MDI/mdipp-1.0.1/mdipp N ${files[$id_a]} N ${files[$id_b]} -n ${num_iter} -t ${thin} -s ${seed} > ${output}

    done
done
