#!/usr/bin/env python
'''
First, all known AD SNPs are added to an id set and the finished targets set.
Next, each line of the ld file with all the scores is parsed to find any line with one of the AD snp IDs.
If a SNP id is found, both snps are added to the target set if the LD is >= 0.8
The target set results in all AD SNPs and any SNP in high LD with any of them. 
Target set is written to target_snps.txt
'''
print("Extracting SNP ids in high LD with AD SNPs")
id_file = open("AD_snp_ids.txt", "r") 
ld_scores = open("gwas.ld", "r") 
out = open("AD_LD_snps_with_scores.txt", "w")
target = open("target_snps.txt", "w")
ids = set()
targets = set()
for line in id_file:
	ids.add(line.strip())
	targets.add(line.strip())
for line in ld_scores:
	fields = line.strip().split()
	if fields[2] in ids and float(fields[6]) >= 0.8:
		out.write('\t'.join(fields) + '\n')
		targets.add(fields[2])
		targets.add(fields[5])
	if fields[5] in ids and float(fields[6]) >= 0.8:
		out.write('\t'.join(fields) + '\n')
		targets.add(fields[2])
		targets.add(fields[5])
print("Done.")
for t in targets:
	target.write(t+'\n')
out.close()
id_file.close()
ld_scores.close()
target.close()
