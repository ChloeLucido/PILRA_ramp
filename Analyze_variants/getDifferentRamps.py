#! /usr/bin/env python
import sys

regGenes = open(sys.argv[1])
ramps = open(sys.argv[2])
output = open(sys.argv[3],'w')

normalToDif = {}
difToNormal = {}
lastNormal = ""
header = regGenes.readline()
while header !="":
    seq = regGenes.readline()
    if header[1:4] =='ID=':
        #if not header in normalToDif:
        #    normalToDif[header] = set()
        lastNormal = header
    else:
        difToNormal[header] = lastNormal
        if lastNormal in normalToDif:
            print (lastNormal + " is already in dictionary!!")
        normalToDif[lastNormal] = header
    header = regGenes.readline()
regGenes.close()
normalToRamp = {}
header = ramps.readline()
while header != "":
    seq = ramps.readline()
    if header[1:4] =='ID=':
        normalToRamp[header] = seq

    header = ramps.readline()

ramps.close()
ramps = open(sys.argv[2])

header = ramps.readline()
allUsedNormalHeaders = set()
while header != "":
    seq = ramps.readline()
    if header[1:4] !='ID=':
        normalHeader = difToNormal[header]
        allUsedNormalHeaders.add(normalHeader)
        if normalHeader in normalToRamp:
            if len(seq) != len(normalToRamp[normalHeader]):
                output.write(normalHeader)
                output.write(normalToRamp[difToNormal[header]])
                output.write(header)
                output.write(seq)
        else:
                output.write(normalHeader)
                output.write("NONE\n")
                output.write(header)
                output.write(seq)
    header = ramps.readline()
for head in normalToRamp:
    if not head in allUsedNormalHeaders:
        output.write(head)
        output.write(normalToRamp[head])
        output.write(normalToDif[head] + "NONE\n")
ramps.close()
output.close()



