#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=00:20:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=selectNmts

#get only orthogroups that contain xiphophorus files bc the '\.t' is found in xipho fasta headers
echo "$(grep -nr '\.t' Results_Jun10/Orthogroup_Sequences)" > xiphophorus_orthogroups.txt

#get the filenames from previous output
cut -d ":" -f -1 xiphophorus_orthogroups.txt > xiphophorusfn.txt
uniq xiphophorusfn.txt > uniq_xiphophorusfn.txt

#copy the xiphophorus files to a new dir
mkdir xiphophorus_groups

echo "moving files to xiphophorus_groups..."

for file in $(<uniq_xiphophorusfn.txt)
do cp "$file" xiphophorus_groups
done

echo "DONE moving files to xiphophorus_groups"

#select n-mt genes based off of file containing mitocarta gene names and move to new dir
echo "$(grep -nrf AllNmt_UniprotID_List.txt xiphophorus_groups)" > nmtxiphophorus.txt
cut -d ":" -f -1 nmtxiphophorus.txt > nmtxiphophorusfn.txt
uniq nmtxiphophorusfn.txt > uniq_nmtxiphophorusfn.txt

echo "adding files to xipho_nmt..."

#mkdir /home/skuster@colostate.edu/projects/mtnuc/xiphophorus/orthoTest/OrthoFinder/divideOGs/xipho_nmt
mkdir xipho_nmt 

for file in $(<uniq_nmtxiphophorusfn.txt)
do mv "$file" xipho_nmt  #/home/skuster@colostate.edu/projects/mtnuc/xiphophorus/orthoTest/OrthoFinder/divideOGs/xipho_nmt
done #this has separated non-n-mt (xiphophorus_groups) from n-mt (xipho_nmt)

echo "DONE adding files to xipho_nmt"

#select direct n-mt genes and add to new dir
echo "$(grep -nrf InteractingNmt_UniprotID_List.txt xipho_nmt)" > directNmts.txt
cut -d ":" -f -1 directNmts.txt > directNmtsxiphophorusfn.txt
uniq directNmtsxiphophorusfn.txt > uniq_directNmtsxiphophorusfn.txt

echo "adding files to xipho_direct_nmt"

#mkdir /home/skuster@colostate.edu/projects/mtnuc/xiphophorus/orthoTest/OrthoFinder/divideOGs/xipho_direct_nmt
mkdir xipho_direct_nmt

for file in $(<uniq_directNmtsxiphophorusfn.txt)
do mv "$file" xipho_direct_nmt #/home/skuster@colostate.edu/projects/mtnuc/xiphophorus/orthoTest/OrthoFinder/divideOGs/xipho_direct_nmt
done #this has separated non-direct n-mt (xipho_nmt) from direct interacting n-mt (xipho_direct_nmt)

echo "DONE adding files to xipho_direct_nmt"

