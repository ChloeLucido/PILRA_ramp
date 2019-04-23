#! /usr/bin/env python
#This script will take the PILRA file and make files for 10%, 20%, etc. for both synonymous and nonsynonymous codon changes
#Input 1 - PILRA fasta file
#Input 2 - an output directory

import sys
from Bio.Seq import Seq

infile = open(sys.argv[1],"r")

header = infile.readline().strip()
seq = infile.readline()
seq = seq[:len(seq) - len(seq) % 3]
infile.close()

#########Ok now lets start permuting our random genes 10% at a time

length = len(seq)
tenth = len(seq) // 10
currentIndex = 1
synFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_synonymous_permutations.fasta", "w")
nonFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_non-synonymous_permutations.fasta", "w")

print("Seq: ", seq)
print("Len: ", len(seq))
for n in range(len(seq)):
	if (n // tenth) + 1 > currentIndex and (n // tenth) + 1 < 11: 
		currentIndex += 1
		synFile.close()
		nonFile.close()
		synFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_synonymous_permutations.fasta", "w")
		nonFile = open(sys.argv[2] + "/" + str(currentIndex * 10) + "%_non-synonymous_permutations.fasta", "w")

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
			
		#print("Codon: ", codon)
		#print("NewCodon: ", newCodon)
		aa1 = Seq(codon).translate()
		aa2 = Seq(newCodon).translate()
		#print("aa1: ", aa1)
		#print("aa2: ", aa2)
		if str(aa1) == str(aa2):
			print("Synonymous")
			synFile.write(">ChangedPos_" + str(n) + "_from_" + seq[n] + "_to_" + d + header + "\n")
			synFile.write(seq[:n] + d + seq[n+1:] + "\n")
		else:
			print("Non synonymous")
			nonFile.write(">ChangedPos_" + str(n) + "_from_" + seq[n] + "_to_" + d + header + "\n")
			nonFile.write(seq[:n] + d + seq[n+1:] + "\n")



