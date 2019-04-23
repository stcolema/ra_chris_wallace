#!/usr/bin/env bash

# Default values for parameters

# Files to edit
# FILES="$(ls *.csv)"

# python file
PY_FILE="/home/MINTS/sdc56/Desktop/ra_chris_wallace/Data_prep/row_names_to_string.py"

# Directory to save output to
SAVE_DIR="./"

# Receive parameters from command line 
while getopts p:o: option; do
    case "${option}" in
        # f) FILES=${OPTARG};;
        p) PY_FILE=${OPTARG};;
        o) SAVE_DIR=${OPTARG};;
    esac
done

# Loop over all files unique combinations performing
# MDI on each combination and saving to a unique file
for i in *.csv; do

    # Create the output file name
	y=${i%.*}
	y=${y##*/}

    output="${SAVE_DIR}${y}_matlab.csv"

    # call="python3 ${PY_FILE} ${i} ${output}"

    # echo "${call}"

    # Call script
    python3 "${PY_FILE}" "${i}" "${output}"

done
