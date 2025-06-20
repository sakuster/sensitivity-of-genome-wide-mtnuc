#!/bin/bash

ind1=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 1-10000)
cut -f $ind1 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset01.tsv
cut -f $ind1 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset01.tsv

ind2=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 10001-20000)
cut -f $ind2 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset02.tsv
cut -f $ind2 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset02.tsv

ind3=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 20001-30000)
cut -f $ind3 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset03.tsv
cut -f $ind3 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset03.tsv

ind4=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 30001-40000)
cut -f $ind4 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset04.tsv
cut -f $ind4 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset04.tsv

ind5=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 40001-50000)
cut -f $ind5 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset05.tsv
cut -f $ind5 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset05.tsv

ind6=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 50001-60000)
cut -f $ind6 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset06.tsv
cut -f $ind6 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset06.tsv

ind7=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 60001-70000)
cut -f $ind7 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset07.tsv
cut -f $ind7 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset07.tsv

ind8=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 70001-80000)
cut -f $ind8 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset08.tsv
cut -f $ind8 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset08.tsv

ind9=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 80001-90000)
cut -f $ind9 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset09.tsv
cut -f $ind9 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset09.tsv

ind10=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 90001-100000)
cut -f $ind10 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset10.tsv
cut -f $ind10 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset10.tsv

ind11=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 100001-110000)
cut -f $ind11 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset11.tsv
cut -f $ind11 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset11.tsv

ind12=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 110001-120000)
cut -f $ind12 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset12.tsv
cut -f $ind12 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset12.tsv

ind13=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 120001-130000)
cut -f $ind13 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset13.tsv
cut -f $ind13 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset13.tsv

ind14=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 130001-140000)
cut -f $ind14 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset14.tsv
cut -f $ind14 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset14.tsv

ind15=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 140001-150000)
cut -f $ind15 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset15.tsv
cut -f $ind15 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset15.tsv

ind16=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 150001-160000)
cut -f $ind16 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset16.tsv
cut -f $ind16 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset16.tsv

ind17=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 160001-170000)
cut -f $ind17 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset17.tsv
cut -f $ind17 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset17.tsv

ind18=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 170001-180000)
cut -f $ind18 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset18.tsv
cut -f $ind18 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset18.tsv

ind19=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 180001-190000)
cut -f $ind19 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset19.tsv
cut -f $ind19 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset19.tsv

ind20=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 190001-200000)
cut -f $ind20 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset20.tsv
cut -f $ind20 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset20.tsv

ind21=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 200001-210000)
cut -f $ind21 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset21.tsv
cut -f $ind21 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset21.tsv

ind22=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 210001-220000)
cut -f $ind22 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset22.tsv
cut -f $ind22 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset22.tsv

ind23=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 220001-230000)
cut -f $ind23 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset23.tsv
cut -f $ind23 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset23.tsv

ind24=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 230001-240000)
cut -f $ind24 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset24.tsv
cut -f $ind24 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset24.tsv

ind25=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 240001-250000)
cut -f $ind25 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset25.tsv
cut -f $ind25 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset25.tsv

ind26=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 250001-260000)
cut -f $ind26 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset26.tsv
cut -f $ind26 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset26.tsv

ind27=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 260001-270000)
cut -f $ind27 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset27.tsv
cut -f $ind27 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset27.tsv

ind28=$(cat CALL_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 270001-287834)
cut -f $ind28 ancestry-probs-par1_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par1CALLdata_subset28.tsv
cut -f $ind28 ancestry-probs-par2_allchrs_CALL_alladults_hmmFebruary2021.tsv > forLD_nonnmts_par2CALLdata_subset28.tsv
