#! /usr/bin/env python
import sys
import os
#This script will calculate the ramps changed. 
#input is ../data/grch_longest_isoforms_ramps.fa

###Make a dictionary of the original ramps
original = open(sys.argv[1],"r")
originalRamps = {}
geneName = ""
for line in original:
	if line[0] == ">":
		geneName = line.split("gene=")[1].split(";")[0]
	else:
		originalRamps[geneName] = line.strip()


allSequences = {}
filesToLines = {}
synLocToPercent = {}
nonSynLocToPercent = {}
for filename in os.listdir("."):
	if ".fasta" not in filename:
		continue
	if "ramps" not in filename:
		seqs = []
		infile = open(filename,"r")
		numLines = 0
		for line in infile:
			if line[0] != ">":
				seqs.append(line.strip())
				numLines += 1
		filesToLines[filename] = numLines
		infile.close()
		allSequences[filename] = seqs
		continue

print("dic: ", filesToLines)

for filename in os.listdir("."):

	if "ramps" not in filename or "fasta" not in filename:
		continue
	print(filename)
	infile = open(filename,"r")
	numSame = 0
	numChanged = 0
	numDestroyed = 0
	numCreated = 0
	for line in infile:
		if line[0] == ">":
			#Get the original ramp
			geneName = line.split("gene=")[1].split(";")[0]
			if geneName not in originalRamps:
				print(geneName + " does not have a ramp. ")
				originalRamp = ""
			else:
				originalRamp = originalRamps[geneName]
			continue
		if len(line) == len(originalRamp):
			numSame += 1
		elif:
			len(originalRamp) == 0:
				numCreated += 1
		else:
			numChanged += 1

	print("Number of unchanged ramps: ", numSame)
	originalName = filename.replace("ramps_","")
	location = int(originalName.split("%")[0])
	total = filesToLines[originalName]
	if numSame + numChanged < total:
		numSame += total - numSame - numChanged
	percentChanged = 1 - (float(numSame) / total)
	print("Percent of ramps changed: ", percentChanged)
	if "non-synonymous" in filename:
		nonSynLocToPercent[location] = percentChanged
	else:
		synLocToPercent[location] = percentChanged
	
	infile.close()

print("Syn Loc to percent: ", synLocToPercent)
print("Non syn loc to percent: ", nonSynLocToPercent)


outfile = open("super_permutations_data.tsv", "w")
outfile.write("Location\tPercent\tMutation\n")
for x in range(10,101,10):
	outfile.write(str(x) + "\t" + str(synLocToPercent[x]) + "\t" +  "Synonymous\n")
			
			
for x in range(10,101,10):
	outfile.write(str(x) + "\t" + str(nonSynLocToPercent[x]) + "\tNon-synonymous\n")

outfile.close()

