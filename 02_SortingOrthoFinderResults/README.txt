This directory contains how I split the Xiphophorus genes into non-n-mt, non-interacting n-mt, and interacting n-mt gene lists. 

First, OrthoFinder was ran containing proteins from Xiphophorus, Homo sapiens, and Mus musculus. 

Then, I used MitoCarta 3.0 to get the n-mt gene lists (both interacting and not). These were converted into uniprot IDs based on a blast and subsequent categorization by Dan's perl script and my postMCtoUP.R script. 

The AllNmt... and InteractingNmt... gene lists were used to pull xiphophorus protein containing orthogroup files into either interacting or non-interacting (those remaining from the AllNmt list after the interacting ones were moved to a new dir)

The selectOGs.sh was used to move all of the files. 

The Results_Jun10 directory was from my original run of OrthoFinder in 2021. 

All the rest of this directory was created/ran in Aug 2024. 
