#!/usr/bin/env bash

# Default values for parameters

# Number of iterations for MDI to perform
num_iter=800

# Thinning variable for MDI
thin=1

# Directory to save output to
save_dir="./output/"

# Receive parameters from command line 
while getopts n:t:o: option; do
    case "${option}" in
        n) num_iter=${OPTARG};;
        t) thin=${OPTARG};;
        o) save_dir=${OPTARG};;
    esac
done

# Loop over all files unique combinations performing
# MDI on each combination and saving to a unique file
for ((i=1; i<11; i++)); do

    # The output file name
    output="${save_dir}out_seed_${i}.csv"

    # Call MDI
    ~/Desktop/MDI/mdipp-1.0.1/mdipp N CD14.csv N CD15.csv N CD19.csv N CD4.csv N CD8.csv N IL.csv N PLA.csv N RE.csv N TR.csv -n ${num_iter} -t ${thin} -s ${i} > ${output}

done
