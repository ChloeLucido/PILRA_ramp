#! /usr/bin/env python
#This script will take the PILRA file and make files for 10%, 20%, etc. for both synonymous and nonsynonymous codon changes
#Input - PILRA fasta file (Located in data/pilra.fasta)

import sys
from Bio.Seq import Seq
import os
if not os.path.exists('data/pilra_permutations'):
        os.makedirs('data/pilra_permutations')
infile = open(sys.argv[1],"r")

header = infile.readline().strip()
seq = infile.readline().upper()
seq = seq[:len(seq) - len(seq) % 3]
infile.close()

######### Permute the random genes 10% at a time

length = len(seq)
tenth = len(seq) // 10
currentIndex = 1
synFile = open("data/pilra_permutations/" + str(currentIndex * 10) + "%_synonymous_permutations.fasta", "w")
nonFile = open("data/pilra_permutations/" + str(currentIndex * 10) + "%_non-synonymous_permutations.fasta", "w")

for n in range(len(seq)):
    if (n // tenth) + 1 > currentIndex and (n // tenth) + 1 < 11: 
        currentIndex += 1
        synFile.close()
        nonFile.close()
        synFile = open("data/pilra_permutations/" + str(currentIndex * 10) + "%_synonymous_permutations.fasta", "w")
        nonFile = open("data/pilra_permutations/" + str(currentIndex * 10) + "%_non-synonymous_permutations.fasta", "w")

    nucs = ['A','T','C','G']
    nucs.remove(seq[n])
    for d in nucs:
        if (n % 3) == 0:
            codon = seq[n: n + 3]
            newCodon = d + seq[n+1: n+3]
        elif (n % 3) == 1:
            if n + 2 >= length:
                newCodon = seq[n-1] + d + seq[n + 1:]
            else:
                newCodon = seq[n - 1] + d + seq[n+2]
        elif (n % 3) == 2:
            if n + 1 >= length:
                newCodon = seq[n-2: n] + d
            else:
                newCodon = seq[n - 2: n] + d
            
        aa1 = Seq(codon).translate()
        aa2 = Seq(newCodon).translate()
        if str(aa1) == str(aa2):
            synFile.write(">ChangedPos_" + str(n) + "_from_" + seq[n] + "_to_" + d + header + "\n")
            synFile.write(seq[:n] + d + seq[n+1:] + "\n")
        else:
            nonFile.write(">ChangedPos_" + str(n) + "_from_" + seq[n] + "_to_" + d + header + "\n")
            nonFile.write(seq[:n] + d + seq[n+1:] + "\n")

