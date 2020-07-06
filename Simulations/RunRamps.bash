#!/bin/bash
# This script finds ramps in all files of a directory
# Input is a directory path

while read FN
do
    python3 ExtRamp.py -i $1/${FN} -o $1/ramps_${FN} -u data/GRCh38_longestIsoforms.fa
done < <(find $1 -maxdepth 1 -type f -name '*.fasta' -printf "%f\n")
