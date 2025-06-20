# -*- coding: utf-8 -*-
#!/usr/bin/python

"""
Created on Tue Nov 15 13:53:41 2022

@author: Shady
"""

import sys
import pandas as pd
#import os

"""
HOW TO USE THIS SCRIPT
from the command line, run the following: 
python3 SNPonly_forHigherComputation_v4.py xipho_transcript_info_REDO.tsv [long version of population specific SNPnames]

argv[1] #is the xipho_transcript file
argv[2] #is the population specific long_par1..._SNPnames.tsv file

"""
#must import os to use the following line:
#os.chdir("C:/Users/Shady/Documents/SchumerMtNucData") #MAY NOT NEED, DEF DON'T NEED FOR OVIS

transfn = sys.argv[1] #"head500_xipho_transcript_info_REDO.tsv" #will be replaced with sys.argv[1]
popfn = sys.argv[2] #"head500_long_par1_CALL_SNPnames.tsv" #will be replaced w/ sys.argv[2]

transcriptDf = pd.read_table(transfn, sep = '\t') #make df from xipho_transcript file
popDf = pd.read_table(popfn, sep = '\t', names = ["SNP_scaffold", "locus"]) #this is reading in the first line as the header

print(transcriptDf) # will print the whole df
print(transcriptDf.loc[1,'seqname']) #this is how you print a particular value
print(transcriptDf.loc[[1]]) #will print all cols for a row

numPop = len(popDf)
numTrans = len(transcriptDf)

#iterate through SNPs and compare to transcript file
matches = []

for i in range(0, numPop):
    bp = popDf.loc[i, 'locus'] #get bp from population specific file
    scaf = popDf.loc[i, 'SNP_scaffold'] #get scaffold name from pop specific file
    for j in range(0, numTrans):
        scaffold = transcriptDf.loc[j, 'seqname'] #get scaffold name from transcript file for comparison
        if scaffold == scaf:
            b = transcriptDf.loc[j, 'begin']
            e = transcriptDf.loc[j, 'end']
            #does the SNP fall within the transcript start/end?
            if bp >= b and bp <= e: 
                group = [scaf, bp, transcriptDf.loc[j, 'attribute'], transcriptDf.loc[j, 'full_geneName'], b, e, transcriptDf.loc[j, 'OGfn'], transcriptDf.loc[j, 'class']]
                matches.append(group)
            else:
                next
        else:
            next
    if len(matches) != i + 1:
        group = [scaf, bp, 'NA', 'NA', 'NA', 'NA', 'NA', 'intergenic']
        matches.append(group)
    print(i)
                
#write the results to file
finColNames = 'SNP_scaffold' + '\t' + 'locus' + '\t' + 'attribute' + '\t' + 'full_geneName' + '\t' + 'begin' + '\t' + 'end' + '\t' + 'OGfn' + '\t' + 'class' + '\t' + 'col_index'
outfn = popfn.split(sep='_')[-2] + "_SNPinfo.tsv"

with open(outfn, "w") as file:
    #write header
    file.write(finColNames)
    file.write('\n')
    
    #iteratively write the matches to the file
    for line in matches:
        for item in line:
            file.write(str(item))
            file.write('\t')
        file.write('\n')
