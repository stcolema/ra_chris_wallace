#!/usr/bin/env python

"""
Author: Stephen D. Coleman (the D is for David)
Studnet nr: 940309160050
Dedication: This script goes out to all of the boys
Description: this is a script to generate csv files with row names appended 
by an X (and thus autmoatically a string)
"""
# Import statements
from sys import argv


def row_names_to_string(file_name: str,
                        out_file_name: str,
                        appendage: str = "X",
                        sep: str = ","):
    """Replaces all sep with itself twice (for matlab version of MDI)
    """

    # Open the inputted files
    with open(file_name, 'r') as f:
        with open(out_file_name, 'w+') as out_f:

            # iterate over the lines in f
            for i, line in enumerate(f):

                # if the header line, write directly to file with no changes
                if i == 0:
                    out_f.write(line)
                    continue

                # Else find the first entry and add the appendage
                row_name_index = line.find(sep)
                new_row_name = line[0:row_name_index] + appendage
                new_line = new_row_name + line[row_name_index:]
                out_f.write(new_line)

if __name__ == "__main__":
    file_name = argv[1]
    out_file_name = argv[2]
    sep = ","
    app = "X"

    if(len(argv) > 3):
        app = argv[3]

        if(len(argv) > 4):
            sep = argv[4]

    row_names_to_string(file_name, out_file_name, app, sep)