#!/usr/bin/env bash

# Find the current data and use this as the default directory name
CURR_DATE=`date +%F`
NEW_DIR="Data_${CURR_DATE}/"

# Instruction to subset data with gene sets from Chris
SUBSET="FALSE"

# Options for subsetting data
ALLOWED_SUBSET_OPTIONS=("FALSE" "bignet" "midnet" "smallnet" "all")

# Flag indicating if subset input is one of allowed options
VIABLE_SUBSET_FLAG=0

# If this is 0 we skip the download and unzip step
DOWNLOAD=1

# Receive parameters from command line 
while getopts d:s: option; do
    case "${option}" in
        d) NEW_DIR=${OPTARG};;
        s) SUBSET=${OPTARG};;
    esac
done

# Check if the subset argument is in one of the viable options
for i in "${ALLOWED_SUBSET_OPTIONS[@]}"
do
    if [ "${i}" == "${SUBSET}" ] ; then
        VIABLE_SUBSET_FLAG=1
    fi
done

# if the subset argument is not a viable option, exit the program with an error message
if [[ ${VIABLE_SUBSET_FLAG} -eq 0 ]]; then
	echo 'SUBSET (-s) variable must be one of "FALSE", "bignet", "midnet", "smallnet" or "all".'
	exit 1
fi

# Download the data
if [[ ${DOWNLOAD} -eq 1 ]]; then
	Rscript ra_chris_wallace/Data_prep/download_data.R -d ${NEW_DIR}

	# Move to the new directory holding the downloaded data
	cd ${NEW_DIR}

	# Unzip this data
	gunzip *.gz
else 
	# Move to the new directory holding the downloaded data
	cd ${NEW_DIR}
fi

# If requested to subset, continue to subset (currently not supported)
if [[ "${SUBSET}" != "FALSE" ]]; then
	Rscript /home/MINTS/sdc56/Desktop/ra_chris_wallace/Data_prep/gene_subsetter_original_data.R -a TRUE -d ./ -e .txt -w ./Subsetted_data/ -g "${SUBSET}" --gene_data /home/MINTS/sdc56/Desktop/ra_chris_wallace/Data/stephen-genesets.RData -p /home/MINTS/sdc56/Desktop/ra_chris_wallace/Data/full_probe_key.csv -t TRUE

	if [[ "${SUBSET}" == "all" ]]; then
		NEW_SUB_DIRS=("bignet" "midnet" "smallnet")
	else
		NEW_SUB_DIRS=(${SUBSET})
	fi
fi

exit 1

# Transpose the data and remove NAs
Rscript ../ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE -d . -w ./Transposed_data/ -n TRUE --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v FALSE

# Make a directory for meta data and move .csv recording observations lost to it
mkdir Meta_data
mv Transposed_data/Observations_lost.csv Meta_data/

# Move to directory containing transposed data and re-arrange row order, 
# fill empty probes with 0 value and record probes present
cd Transposed_data
Rscript /home/MINTS/sdc56/Desktop/ra_chris_wallace/Data_prep/finad_all_probe_ids.R -d .

# Move back up a level and create a directory for this rearranged and filled data to it
cd ..
mkdir Filled_data
mv Transposed_data/*_filled.csv Filled_data/
mv Transposed_data/probe* Meta_data/

# Make a directory for the original downloaded data
mkdir Original_data
mv *.txt Original_data
