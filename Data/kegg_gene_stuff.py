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



def find_gene_sets(set_names: List[str],
    file_name: str,
    out_file_name: str) -> Dict[str, List[str]]:

    gene_dict = {}
    with open(file_name, 'r') as f:
        with open(out_file_name, 'x') as o:
            for line in f:


                line_of_i = line.split()

                x = [True if p in line_of_i[0] else False for p in pathways]

                path_to_drop = list(compress(pathways, x))

                if not any(x):
                    continue


                pathways.remove(path_to_drop[0])

                gene_dict.update( {line_of_i[0] : line_of_i[2:]} )


    return(gene_dict)

if __name__ == "__main__":

    file_name = "/home/MINTS/sdc56/Desktop/kegg_msigdb_2.txt"
    out_file_name = "/home/MINTS/sdc56/Desktop/kegg_genes.csv"
    gen_url = "http://www.broadinstitute.org/gsea/msigdb/cards/"
    gene_dict = {}
    pathways = ['NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY',
                'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY'
                ]

    my_genes = find_gene_sets(pathways, file_name, out_file_name)

    print(my_genes.keys())