#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:52:33 2023

@author: kusters
"""

#import sys

#THIS CALCULATES MT ANCESTRY PER INDIVIDUAL AND DETERMINES MT HAPLOTYPE

#par1fn = sys.argv[1]
par1fn = "ancestry-probs-par1_allchrs_CHAF_mito_Feb2021.tsv"
#par2fn = sys.argv[2]
par2fn = "ancestry-probs-par2_allchrs_CHAF_mito_Feb2021.tsv"

outfn1 = par1fn.split(sep = '.')
outfn1 = outfn1[0] + "FILTERED.tsv"

outfn2 = par2fn.split(sep = '.')
outfn2 = outfn2[0] + "FILTERED.tsv"
outfn3 = "CHAF_ancestryResults.tsv"

# anc = pd.read_table(ancestryfn, header = 0, nrows = 5)

with open(par1fn, "r") as anc1, open(par2fn, "r") as anc2, open("CHAF_MtAncestryCall.tsv", "w") as outf1:#, open(outfn1, "w") as file1:#, open(outfn2, "w") as file2, open (outfn3, "w") as resultFile:
    head = anc1.readline() #skip header line for reading
    #file1.write(head) #write header of OG file to output file
    
    head2 = anc2.readline() #do same for ancestry file for par2
    #file2.write(head2)
    
    outf1.write("FishID\tPar1 Avg Mt Ancestry\tPar2 Avg Mt Ancestry\tMajor Mt Ancestry\tConclusion\n")
    
    #resultFile.write("FishID\t Percent Xbir Ancestry\n")
    
    # goodFish = []
    count = 0 #limit num rows read for testing
    excluded1 = []
    excluded2 = []
    
    popAvgAnc = 0
    
    
    #iterate through ancestry file line by line
    for line1, line2 in zip(anc1, anc2):
        count = count + 1
        mtancestry = 0
        hetTotal = 0
        par2Total = 0
        
        #wholeline = anc1.readline()
        wholeline = line1
        stripped = wholeline.strip()
        modline = stripped.split(sep = '\t')
        
        #wholeline2 = anc2.readline()
        wholeline2 = line2
        stripped2 = wholeline2.strip()
        modline2 = stripped2.split(sep = '\t')
        
        # if (str(modline[0]) != "HUEX-XI-19-80.R1.fastq"): #testing bc this fish is being plotted incorrectly. nee to remove this if/else (and break at end) when done
        #     continue
        # else:
        
        anc1Probs = [] #collect probabilities for anc1
        anc2Probs = [] #collect probs for anc2
        hetProbs = [] #collect probs for heterozygous loci
        genomAnc = 0
        
        # resultFile.write(str(modline[0])) #write fishID
        # resultFile.write("\t")
        
        outf1.write(str(modline[0])) #writes fish ID
        outf1.write("\t")
        
        mtsnpcount = 0
        
        for i in range(1, len(modline)): #skip i = 0 bc it's fish ID, str
            if (i < 628790 - 1): #the SNP number that mtSNPs start on
                continue
            else:
                if (float(modline[i]) < 0.8 and float(modline[i]) > 0.2 ): #skips any par1 posterior probability that is between 0.2 and 0.8
                    continue
                if (float(modline2[i]) < 0.8 and float(modline2[i]) > 0.2 ): #skips any par2 posterior probability that is between 0.2 and 0.8
                    continue
                mtsnpcount = mtsnpcount + 1
                #mtancestry = mtancestry + modline[i]
                het = 1 - (float(modline[i]) + float(modline2[i]))
                hetTotal = hetTotal + het
                par2Total = par2Total + float(modline2[i]) #+ (het / 2)
                # if (het > float(modline[i]) and het > float(modline2[i])):
                #     print(modline[0], " has heterozygous ancestry at mt SNP ", mtsnpcount, " for total of ", het)
                locInher = float(modline[i]) #+ (het / 2)
                mtancestry = mtancestry + locInher      
        if (mtsnpcount == 0):
            outf1.write("NA")  #writes par1 ancestry
            outf1.write("\t")
            outf1.write("NA") #writes par2 ancestry
            outf1.write("\t")
            outf1.write("NA") #writes major ancestry
            outf1.write("\t")
            outf1.write("No mtSNPs\n") #writes conclusion
            continue
            #print(modline[0])
            #sys.exit(1)
        mtancestry = mtancestry / mtsnpcount
        hetTotal = hetTotal / mtsnpcount
        par2Total = par2Total / mtsnpcount
        
        outf1.write(str(mtancestry))  #writes par1 ancestry
        outf1.write("\t")
        outf1.write(str(par2Total))
        outf1.write("\t")
        
        print(str(mtancestry))
        
        #if(mtsnpcount )
        
        if (mtancestry > 0.8):
            #print("went to > 0.8)")
            outf1.write(str(mtancestry)) #writes major mt ancestry
            outf1.write("\t")
            outf1.write("Par1\n") #for HC, should be Xcor bc par 1 is Xbir
        elif (hetTotal > mtancestry and hetTotal > par2Total):
            outf1.write(str(hetTotal))
            outf1.write("\t")
            outf1.write("Heterozygous\n")
        elif (par2Total > 0.8):
            outf1.write(str(par2Total))
            outf1.write("\t")
            outf1.write("Par2\n")
        else:
            outf1.write(str(mtancestry))
            outf1.write("+")
            outf1.write(str(hetTotal))
            outf1.write("+")
            outf1.write(str(par2Total))
            outf1.write("\t")
            outf1.write("No assignment\n")
        #break