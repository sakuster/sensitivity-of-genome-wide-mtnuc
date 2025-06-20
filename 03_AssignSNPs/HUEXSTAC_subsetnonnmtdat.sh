#!/bin/bash

ind1=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 1-10000)
cut -f $ind1 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset01.tsv
cut -f $ind1 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset01.tsv

ind2=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 10001-20000)
cut -f $ind2 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset02.tsv
cut -f $ind2 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset02.tsv

ind3=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 20001-30000)
cut -f $ind3 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset03.tsv
cut -f $ind3 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset03.tsv

ind4=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 30001-40000)
cut -f $ind4 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset04.tsv
cut -f $ind4 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset04.tsv

ind5=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 40001-50000)
cut -f $ind5 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset05.tsv
cut -f $ind5 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset05.tsv

ind6=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 50001-60000)
cut -f $ind6 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset06.tsv
cut -f $ind6 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset06.tsv

ind7=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 60001-70000)
cut -f $ind7 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset07.tsv
cut -f $ind7 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset07.tsv

ind8=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 70001-80000)
cut -f $ind8 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset08.tsv
cut -f $ind8 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset08.tsv

ind9=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 80001-90000)
cut -f $ind9 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset09.tsv
cut -f $ind9 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset09.tsv

ind10=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 90001-100000)
cut -f $ind10 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset10.tsv
cut -f $ind10 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset10.tsv

ind11=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 100001-110000)
cut -f $ind11 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset11.tsv
cut -f $ind11 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset11.tsv

ind12=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 110001-120000)
cut -f $ind12 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset12.tsv
cut -f $ind12 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset12.tsv

ind13=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 120001-130000)
cut -f $ind13 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset13.tsv
cut -f $ind13 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset13.tsv

ind14=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 130001-140000)
cut -f $ind14 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset14.tsv
cut -f $ind14 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset14.tsv

ind15=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 140001-150000)
cut -f $ind15 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset15.tsv
cut -f $ind15 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset15.tsv

ind16=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 150001-160000)
cut -f $ind16 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset16.tsv
cut -f $ind16 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset16.tsv

ind17=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 160001-170000)
cut -f $ind17 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset17.tsv
cut -f $ind17 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset17.tsv

ind18=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 170001-180000)
cut -f $ind18 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset18.tsv
cut -f $ind18 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset18.tsv

ind19=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 180001-190000)
cut -f $ind19 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset19.tsv
cut -f $ind19 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset19.tsv

ind20=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 190001-200000)
cut -f $ind20 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset20.tsv
cut -f $ind20 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset20.tsv

ind21=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 200001-210000)
cut -f $ind21 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset21.tsv
cut -f $ind21 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset21.tsv

ind22=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 210001-220000)
cut -f $ind22 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset22.tsv
cut -f $ind22 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset22.tsv

ind23=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 220001-230000)
cut -f $ind23 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset23.tsv
cut -f $ind23 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset23.tsv

ind24=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 230001-240000)
cut -f $ind24 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset24.tsv
cut -f $ind24 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset24.tsv

ind25=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 240001-250000)
cut -f $ind25 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset25.tsv
cut -f $ind25 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset25.tsv

ind26=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 250001-260000)
cut -f $ind26 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset26.tsv
cut -f $ind26 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset26.tsv

ind27=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 260001-270000)
cut -f $ind27 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset27.tsv
cut -f $ind27 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset27.tsv

ind28=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 270001-280000)
cut -f $ind28 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset28.tsv
cut -f $ind28 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset28.tsv

ind29=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 280001-290000)
cut -f $ind29 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset29.tsv
cut -f $ind29 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset29.tsv

ind30=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 290001-300000)
cut -f $ind30 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset30.tsv
cut -f $ind30 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset30.tsv

ind31=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 300001-310000)
cut -f $ind31 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset31.tsv
cut -f $ind31 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset31.tsv

ind32=$(cat HUEXSTAC_index_nonnmts_SNPs.txt | tr '\n' ',' | cut -d ',' -f 310001-319679)
cut -f $ind32 ancestry-probs-par1_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par1HUEXSTACdata_subset32.tsv
cut -f $ind32 ancestry-probs-par2_allchrs_thinned_HUEX_STAC_cortezi_cluster_March2021.tsv > forLD_nonnmts_par2HUEXSTACdata_subset32.tsv
