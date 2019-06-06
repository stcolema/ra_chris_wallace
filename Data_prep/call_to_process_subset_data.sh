#!/usr/bin/env bash

# Remember to call these from within a directory that is on the Desktop and contains 
# .csv files for the gene set data

SCALE=TRUE
CENTRE=TRUE

# Receive parameters from command line 
while getopts s:c: option; do
    case "${option}" in
        d) SCALE=${OPTARG};;
        s) CENTRE=${OPTARG};;
    esac
done

# Transpose and standardise the data
Rscript ../ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE -d . -w ./Transposed_data/ -n TRUE --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v FALSE -s ${SCALE} -c ${CENTRE} -e ".csv"

# Make a directory for meta data and move .csv recording observations lost to it
mkdir Meta_data
mv Transposed_data/Observations_lost.csv Meta_data/

# Move to directory containing transposed data and re-arrange row order, 
# fill empty probes with 0 value and record probes present
cd Transposed_data

# Fill empty probes with draws from a standard normal
Rscript /home/MINTS/sdc56/Desktop/ra_chris_wallace/Data_prep/finad_all_probe_ids.R -d . -r T

# Move back up a level and create a directory for this rearranged and filled data to it
cd ..
mkdir Filled_data
mv Transposed_data/*_filled.csv Filled_data/
mv Transposed_data/probe* Meta_data/

# Make a directory for the original downloaded data
mkdir Original_data
mv *.csv Original_data

# Rename final form of files and move to new directory
Rscript /home/MINTS/sdc56/Desktop/ra_chris_wallace/Data_prep/rename_files.R -d ./Filled_random_data/ -w ./Final_random_data/