#!/usr/bin/python

import sys

scaffolds = []
bed = []

with open(sys.argv[1], "r") as handle: #bed file with categories
    #is the xipho_transcripts.bed which is a 
    #mod of xipho_transcript_info_REDO.tsv
    m = handle.readline()  #reads header so when you iterate you won't include header
    for line in handle:
        line = line.strip().split('\t')
        line[1] = int(line[1])  #beginnning bp
        line[2] = int(line[2])  #end bp
        if line[0] not in scaffolds:  #line[0] is scaffold name
            scaffolds.append(line[0])
        bed.append(line)

raw_snps = []
with open(sys.argv[2], "r") as handle: #SNP header file in long format, split into scaffold/pos
#read everything into memory 
    for line in handle:
        line = line.strip().split('\t')
        line[1] = int(line[1])  #bp locus on scaffold, class as int
        raw_snps.append(line)

#good_snp = []
for i in scaffolds:
    temp_bed = [j for j in bed if j[0] == i] #all the genes on that scaffold
    temp_snp = [k for k in raw_snps if k[0] == i] #all the SNPs on that scaffold
    for k in temp_bed:
#        print(len(k))
#        print(len(temp_snp[0]))
        good_snp = [":".join(str(elem) for elem in l) for l in temp_snp if l[1] >= k[1] and l[1] <= k[2]] #SNPs that are greater than start, smaller than end of each gene
        
        if good_snp: 
            
            print(*good_snp, sep="\n") #May produce multiple SNP calls if SNP is located on multiple trascripts 


#print out SNP IDs
