#!/usr/bin/env python

"""
Author: Stephen D. Coleman (the D is for David)
Studnet nr: 940309160050
Dedication: This script goes out to all of the boys
Description: this is a script to take all lines from csv files 
"""
# Import statements
from sys import argv
import os
import glob

def find_ext_files(dir_path:str = ".",
                   extension:str = ".csv"):
    
    EXT = "*" + extension

    for path, subdir, files in os.walk(dir_path):

        all_csv_files = [file
                         for file in glob.glob(os.path.join(path, EXT))]
        break

    return(all_csv_files)




def big_csv(input_dir: str,
            out_file_name: str,
            line_to_keep: int = 2,
            extension: str = ".csv"):
    """Replaces all sep with itself twice (for matlab version of MDI)
    """

    input_files = find_ext_files(input_dir, extension)

    # Open the inputted files
    with open(out_file_name, 'w+') as out_f:
        for i, file_name in enumerate(input_files):
            with open(file_name, 'r') as f:

                # iterate over the lines in f
                for j, line in enumerate(f):

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

    if(len(argv) > 3):
        line_to_keep = argv[3]

    if(len(argv) > 4):
        ext = argv[4]

    # print(dir_path)

    # find_ext_files(dir_path, ext)
    big_csv(dir_path, out_file_name, line_to_keep, ext)