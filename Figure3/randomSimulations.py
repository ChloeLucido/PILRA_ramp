#! /usr/bin/env python
#This script will take the hg38 ref file and make files for 10%, 20%, etc. for both synonymous and nonsynonymous codon changes
#Input 1 - hg38 fasta file
#Input 2 - an output directory

import sys
from Bio.Seq import Seq
import random
from random import sample

random.seed(0)

infile = open(sys.argv[1],"r")

genes = []

header = ""
for line in infile:
	line = line.strip()
	if line[0] == ">":
		header = line
	else:
		line = line.upper()
		genes.append([header,line])

randGenes = sample(genes,1000)

####Open all of the files
for x in range(10,101,10):
	synFile = open(sys.argv[2] + "/" + str(x) + "%_synonymous_permutations.fasta", "w")
	nonFile = open(sys.argv[2] + "/" + str(x) + "%_non-synonymous_permutations.fasta", "w")
	synFile = open(sys.argv[2] + "/" + str(x) + "%_synonymous_permutations.fasta", "w")
	nonFile = open(sys.argv[2] + "/" + str(x) + "%_non-synonymous_permutations.fasta", "w")

#########Ok now lets start permuting our random genes 10% at a time
for gene in genes:
	header = gene[0]
	seq = gene[1]
	length = len(seq)
	tenth = len(seq) // 10
	currentIndex = 1
	
	synFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_synonymous_permutations.fasta", "a")
	nonFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_non-synonymous_permutations.fasta", "a")
	for n in range(len(seq)):
		if (n // tenth) + 1 > currentIndex and (n // tenth) + 1 < 11: 
			currentIndex += 1
			synFile.close()
			nonFile.close()
			synFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_synonymous_permutations.fasta", "a")
			nonFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_non-synonymous_permutations.fasta", "a")
	
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



