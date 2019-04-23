#!/usr/bin/env python

"""
Author: Stephen D. Coleman (the D is for David)
Studnet nr: 940309160050
Dedication: This script goes out to all of the boys
Description: this is a script to generate csv files with empty columns between
every filled column (required by matlab MDI)
"""
# Import statements
from sys import argv


def double_separator(file_name: str,
                     out_file_name: str,
                     sep: str = ","):
    """Replaces all sep with itself twice (for matlab version of MDI)
    """
    with open(file_name, 'r') as f:
        with open(out_file_name, 'x') as out_f:
            for line in f:
                new_line = line.replace(sep, 2 * sep)
                out_f.write(new_line)

if __name__ == "__main__":
    file_name = argv[1]
    out_file_name = argv[2]
    sep = ","

    if(len(argv) > 3):
        sep = argv[3]

    double_separator(file_name, out_file_name, sep)