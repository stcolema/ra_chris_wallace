#!/usr/bin/env bash
Rscript ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE -d ./MDI/Data/Original\ data/ -w ./VSN_NA_min_data/ -n TRUE --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v TRUE

