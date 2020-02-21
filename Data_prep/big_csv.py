#!/usr/bin/env python

"""
Author: Stephen D. Coleman (the D is for David)
Dedication: This script goes out to all of the boys
Description: this is a script to take the last lines from csv files and compile
    them into one large file. This is part of the CONSENSUS INFERENCE pipeline.
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


def find_ext_files(dir_path:str = ".",
                   extension:str = ".csv"):

    EXT = "*" + extension

    for path, subdir, files in os.walk(dir_path):

        all_csv_files = [file for file in glob.glob(os.path.join(path, EXT))]
        break

    return(all_csv_files)




def big_csv(input_dir: str,
            out_file_name: str,
            line_to_keep: int = 2,
            header: str = "True",
            extension: str = ".csv",
            input_file_name: str = ".",
            n_dataset: int = 0):
    """Replaces all sep with itself twice (for matlab version of MDI)
    """

    input_files = find_ext_files(input_dir, extension)

    # Open the inputted files
    with open(out_file_name, 'w+') as out_f:

        # If no header is present, create one
        if(header.lower() in ["false", "f"]):
            header_string = create_header(input_file_name, n_dataset, sep = ",")
            out_f.write(header_string)

        for i, file_name in enumerate(input_files):

            with open(file_name, 'r') as f:

                # iterate over the lines in f
                for j, line in enumerate(f):

                    if header.lower() in ["true", "t"]:

                        # if the header line, write directly to file with no changes
                        if j == 0:

                            if i == 0:

                                out_f.write(line)

                            continue


                    if j == line_to_keep:
                        # Else find the first entry and add the appendage
                        out_f.write(line)

                    if j != line_to_keep:
                        continue



if __name__ == "__main__":

    dir_path = argv[1]
    out_file_name = argv[2]
    ext = ".csv"
    line_to_keep = 2
    header = "True"
    input_file_name = ""
    n_dataset = None

    if(len(argv) > 3):
        line_to_keep = int(argv[3])

    if(len(argv) > 4):
        header = argv[4]

    if(len(argv) > 5):
        ext = argv[5]

    if(len(argv) > 6):
        input_file_name = argv[6]

    if(len(argv) > 7):
        n_dataset = int(argv[7])

    if(header.lower() not in ['true', 't', 'false', 'f']):
        raise ValueError("header value is not true or false. Please check.")

    # find_ext_files(dir_path, ext)
    big_csv(dir_path, out_file_name, line_to_keep, header, ext, input_file_name, n_dataset)
