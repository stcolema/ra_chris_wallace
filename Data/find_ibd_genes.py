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
from itertools import compress
from typing import List, Dict
import csv


def find_genes(file_name: str, key_str: str = "               ") -> List[str]:

    # Declare lsit ot hold genes
    genes = []

    # Open the file containing the gene set and iterate over lines
    with open(file_name, 'r') as f:
        for line in f:

            # Find the point at which (RefSeq) occurs, the string
            # preceding the gene names in each line
            first_ind = line.find(key_str) 

            # If this is not found skip to the next line
            if(first_ind == - 1):
                continue

            # Find returns the point at which the substring begins,
            # we want its end, so add its length
            first_ind += len(key_str)

            # The gene names are enclosed by a semi-colon if there is 
            # any further description.
            last_ind = line.find(";")

            # If no semi-colon is present we want the remainder of 
            # the string
            if(last_ind == -1):

                # Strip and split the section of the line containing gene names
                genes += line[first_ind:].strip().split(", ")
            else:
                genes += line[first_ind:last_ind].strip().split(", ")
        

    return(genes)

if __name__ == "__main__":

    in_file = argv[1]
    out_file = argv[2]

    my_genes = find_genes(in_file)
    # my_genes = find_genes("/home/MINTS/sdc56/Desktop/kegg_ibd_genes.txt")
    #my_genes = find_genes("/home/MINTS/sdc56/Desktop/kegg_ibd_genes", key_str = "(RefSeq)")

    # with open("/home/MINTS/sdc56/Desktop/ste_out.csv", "w") as outfile:
    with open(out_file, "w") as outfile:
        outfile.write("Genes\n")
        for entries in my_genes:
            outfile.write(entries)
            outfile.write("\n")
    
    print(my_genes)