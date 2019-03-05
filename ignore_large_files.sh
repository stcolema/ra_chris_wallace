#!/usr/env/bin bash

# Default values for parameters

# Filelocation to add filenames to
file=".gitignore"

# Size threshold to find files exceeding
size="50M"

# Receive parameters from command line 
while getopts f:s: option; do
    case "${option}" in
        f) file=${OPTARG};;
        s) size=${OPTARG};;
    esac
done

# Add files to git ignore (the grep aspect does not work)
find . -size +${size} | grep -qxF ${file} || find . -size +${size} | cat >> ${file}

