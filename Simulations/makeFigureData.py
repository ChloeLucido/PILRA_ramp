#! /usr/bin/env python
import sys
import os
#input is ../reference/grch_longest_isoforms_ramps.fa

###Make a dictionary of the original ramps
original = open("data/GRCh38_longestIsoforms_ramps.fa","r")
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
for filename in os.listdir(sys.argv[1]):
	if ".fasta" not in filename:
		continue
	if "ramps" not in filename:
		headers = []
		infile = open(sys.argv[1] + "/" + filename,"r")
		numLines = 0
		for line in infile:
			if line[0] != ">":
				numLines += 1
			else:
				headers.append(line.strip())
		filesToLines[filename] = numLines
		infile.close()
		allSequences[filename] = headers
		continue

for filename in os.listdir(sys.argv[1]):

	if "ramps" not in filename or "fasta" not in filename:
		continue
	infile = open(sys.argv[1] + "/" + filename,"r")
	numSame = 0
	numChanged = 0
	numCreated = 0
	mutantRamps = []
	for line in infile:
		if line[0] == ">":
			#Get the original ramp
			geneName = line.split("gene=")[1].split(";")[0]
			if geneName not in originalRamps:
				originalRamp = ""
			else:
				originalRamp = originalRamps[geneName]
			mutantRamps.append(line.strip())
			continue
		if len(originalRamp) == 0:
				numCreated += 1
		elif len(line.strip()) != len(originalRamp.strip()):
			numChanged += 1

	numDestroyed = 0
	originalName = filename.replace("ramps_","")
	for h in allSequences[originalName]:
		gene = h.split("gene=")[1].split(";")[0]
		if gene in originalRamps and h.strip() not in mutantRamps:
			numDestroyed += 1
	location = int(originalName.split("%")[0])
	total = filesToLines[originalName]
	numSame = total - numChanged - numCreated - numDestroyed
	percentChanged = float(numChanged) / total
	percentCreated = float(numCreated) / total
	percentDestroyed = float(numDestroyed) / total

	if "non-synonymous" in filename:
		nonSynLocToPercent[location] = [percentChanged,percentCreated,percentDestroyed]
	else:
		synLocToPercent[location] = [percentChanged,percentCreated,percentDestroyed] 
	
	infile.close()


outfile = open(sys.argv[2], "w")
outfile.write("Location\tChanged\tCreated\tDestroyed\tMutation\n")
for x in range(10,101,10):
	percents = synLocToPercent[x]
	outfile.write(str(x) + "\t" + str(percents[0]) + "\t" + str(percents[1]) + "\t" + str(percents[2]) + "\t"  "Synonymous\n")
			
			
for x in range(10,101,10):
	percents = nonSynLocToPercent[x]
	outfile.write(str(x) + "\t" + str(percents[0]) + "\t" + str(percents[1]) + "\t" + str(percents[2]) + "\t"  "Non-Synonymous\n")

outfile.close()

