#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:22:37 2022
Made to filter ancestry files for % ancestry
any individual fish having pure ancestry >20% will be removed

@author: kusters

HOW TO USE THIS:
python3 [name of ancestry file for par1] [fn for ancestry par2]

"""

# import pandas as pd
import sys

#par1fn = sys.argv[1]
par1fn = "ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv"
#par2fn = sys.argv[2]
par2fn = "ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv"

outfn1 = par1fn.split(sep = '.')
outfn1 = outfn1[0] + "FILTERED.tsv"

outfn2 = par2fn.split(sep = '.')
outfn2 = outfn2[0] + "FILTERED.tsv"
outfn3 = "HC_ancestryResults.tsv"

nNucSNPs = 689968 #col number of first mtSNP. 689968 for HC, 629432 for CALL, 628790 for CHAF

mtSNPnum = 689968 #first mtSNP. anything past this will NOT be included in ancestry determination

# anc = pd.read_table(ancestryfn, header = 0, nrows = 5)

with open(par1fn, "r") as anc1, open(par2fn, "r") as anc2, open(outfn1, "w") as file1, open(outfn2, "w") as file2, open (outfn3, "w") as resultFile:
    head = anc1.readline() #skip header line for reading
    file1.write(head) #write header of OG file to output file
    
    head2 = anc2.readline() #do same for ancestry file for par2
    file2.write(head2)
    
    resultFile.write("FishID\t Percent Par1 Ancestry\tPercent Par2 Ancestry\n")
    
    # goodFish = []
    count = 0 #limit num rows read for testing
    excluded1 = []
    excluded2 = []
    
    popAvgAnc = 0
    
    #iterate through ancestry file line by line
    for line1, line2 in zip(anc1, anc2):
        count = count + 1
        #wholeline = anc1.readline()
        wholeline = line1
        stripped = wholeline.strip()
        modline = stripped.split(sep = '\t')
        
        #wholeline2 = anc2.readline()
        wholeline2 = line2
        stripped2 = wholeline2.strip()
        modline2 = stripped2.split(sep = '\t')
        
        anc1Probs = [] #collect probabilities for anc1
        anc2Probs = [] #collect probs for anc2
        hetProbs = [] #collect probs for heterozygous loci
        genomAnc = 0
        par2Anc = 0
        numNucSNPs = 0
        
        resultFile.write(str(modline[0])) #write fishID
        resultFile.write("\t")
        
        for i in range(1, len(modline)): #skip i = 0 bc it's fish ID, str
            # anc1Probs.append(float(modline[i]))
            # anc2Probs.append(float(modline2[i]))
            if (i <= nNucSNPs - 1): 
                if (float(modline[i]) < 0.8 and float(modline[i]) > 0.2 ): #if par1 is not good
                    continue #then skip this whole SNP
                if (float(modline2[i]) < 0.8 and float(modline2[i]) > 0.2): #if par2 is not good
                    continue #then skip this whole SNP
                het = 1 - (float(modline[i]) + float(modline2[i])) 
                # hetProbs.append(het)
                locInher = float(modline[i]) + (het / 2)
                genomAnc = genomAnc + locInher
                numNucSNPs = numNucSNPs + 1
                
                par2 = float(modline2[i]) + (het / 2)
                par2Anc = par2Anc + par2
            else:
                continue
        if (numNucSNPs == 0):
            excluded1.append([modline[0], 0])
            resultFile.write("NA")  #writes par1 ancestry
            resultFile.write("\t")
            resultFile.write("NA\n") #writes par2 ancestry
            continue
        
        genomAnc = genomAnc / (numNucSNPs)
        par2Anc = par2Anc / (numNucSNPs)
        print(len(modline) - 1)
        
        resultFile.write(str(genomAnc)) #write par 1 % ancestry 
        resultFile.write("\t")
        resultFile.write(str(par2Anc)) #write par 2 $ ancestry
        resultFile.write("\n")
        
        popAvgAnc = popAvgAnc + genomAnc
            
        if (genomAnc <= 0.8 and genomAnc >= 0.2):
            # goodFish.append(line)
            file1.write(wholeline) #write non-filtered individual to output file
            file2.write(wholeline2)
        else:
            excluded1.append([modline[0], genomAnc])
            excluded2.append([modline2[0], genomAnc])
            if modline[0] != modline2[0]:
                print("ERROR: not same fish. " + str(modline[0]) + " is not the same as " + str(modline2[0]))
        
        # for i in modline:
        #     if (i != modline[0]):
        #         probs.append(float(i)) #modline[0] is a str (is fish ID), need to not have it
        
        # anc2Probs = []
        # for j in modline2:
        #     if (j != modline2[0]):
        #         anc2Probs.append(float(j))
        
        # anc1Sum = sum(anc1Probs) / (len(modline) - 1) #get avg ancestry probability
        # anc2Sum = sum(anc2Probs) / (len(modline) - 1) #avg ancestry for anc2
        # hetSum = sum(hetProbs) / (len(modline) - 1)
        # print(str(anc1Sum)) #for testing
        
        
        # if (anc1Sum <= 0.8 and anc1Sum >= 0.2):
        #     goodFish.append(line)
        #     file1.write(line) #write non-filtered individual to output file
        # else:
        #     excluded.append(modline[0])
        
        # if (count == 110): #limit num rows read for testing
        #     break

    popAvgAnc = popAvgAnc / count
    resultFile.write("\n\n\n")
    resultFile.write("Average ancestry for this population is: ")
    resultFile.write(str(popAvgAnc))



print(str(len(excluded1)) + " were removed from par1 due to ancestry")
print("They were: " + str(excluded1))
    

 #PROBABLY ACTUALLY NEED TO USE PANDAS DF, NOT THIS
