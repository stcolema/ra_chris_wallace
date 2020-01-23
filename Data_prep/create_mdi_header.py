#!/usr/bin/env python

"""
Author: Stephen D. Coleman (the D is for David)
Dedication: This script goes out to all of the boys
Description: This creates a header for MDI MATLAB output in keeping with the
  format used by Rich Savage
"""
# Import statements
from sys import argv
import os
import glob
from itertools import combinations
from getopt import getopt

def create_header(input_file_name:str, n_dataset:int,
                  sep:str = ","):
    """Create a string of column names
    input_file_name: a string, the name of the file to take the gene names from
    n_dataset: an int, the number of datasets used in the analysis
    sep: a string, the separator between items in the output file (default is comma)
    """
    # The range of datasets (Python is 0 indexed)
    dataset_index = range(1, n_dataset + 1)

    # Define the strings for the mass parameter columns
    mass_params = ["MassParameter_" + str(i) for i in dataset_index]

    # Create the combinations of datasets relevant to the phi parameters
    phi_param_comb = combinations(dataset_index, 2)

    # Define the phi parameter headers, these are of the style Phi13
    phi_params = ["Phi_" + "".join(map(str, x)) for x in phi_param_comb]

    # Find the appropriate gene / sample names from one of the input files
    with open(input_file_name) as f:
        row_names = [row.split(",")[0].replace('"', '') for row in f]

    # the first entry is a blank "", drop this (the entry in a csv file above
    # the row names and beside the column names, cell(0,0))
    row_names = row_names[1:]

    # Define the list of dataset names
    dataset_names = ["Dataset" + str(i) for i in dataset_index]

    # For each object in each dataset create a header for the allocations
    alloc_names = [dt + "_" + entry for dt in dataset_names for entry in row_names]

    # The headers as a list of strings
    header_list = mass_params + phi_params + alloc_names

    # The headers as a string ready to write to a csv file
    header_string = sep.join(header_list)  + "\n"

    return(header_string)

if __name__ == "__main__":

    out_file_name = argv[1]
    input_file_name = argv[2]
    n_dataset = int(argv[3])
    sep = ","

    if(len(argv) > 4):
        sep = argv[4]

    header = create_header(input_file_name, n_dataset, sep)

    with open(out_file_name, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(header.rstrip('\r\n') + '\n' + content)
