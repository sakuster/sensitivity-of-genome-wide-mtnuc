##dN/dS analysis

library(ape)
library(seqinr)
library(ggridges)
library(dplyr)
library(ggplot2)
library(car)
library(ggrepel)

## for HC, ndufs5 pos 31 should have an A for xcor and a M for xbir
#g8054.t1
#ScyDAA6.1934.HRSCAF.2318.2066737
#subs[433844,]
#geneID = 9133


#ndufa13
#g3279.t1
#geneID = 11916
#pos 140 xbir = L, xcor = F
# ScyDAA6.2393.HRSCAF.2888.12131228
#  ScyDAA6.2393.HRSCAF.2888.12131624 
#  ScyDAA6.2393.HRSCAF.2888.12132291 
#  ScyDAA6.2393.HRSCAF.2888.12133227 
# looking for ScyDAA6-2393-HRSCAF-2888:12133015, 16, 17 -- does not exist in sub OR HC anc data

population <- "HUEXSTAC"

#read in substitution data for xcor & xbir
path <- "/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/"

#read in substitution data
subs <- read.csv(paste0(path, "xbir_xcor_10x_genomes_identify_pairwise_dxy.txt"), #read.csv(paste0(path, "Xbirchmanni_Xcortezi_ancestry_informative_sites_addmito_filterF1 copy.txt"),
                 header = F, quote = "", sep = "\t") %>%
  mutate(subType = "", xbirAA = "", xcorAA = "", AAresid = NA, codonPos = NA) %>%
  rename(SNP_scaffold = V1, locus = V2, Xbir = V3, Xcor = V4) 


#read in gtf file
gtf <- read.csv("/Users/kusters/Documents/XiphoStartingOver_Aug2024/00_StartingDataFiles/xiphophorus_birchmanni_10x_12Sep2018_yDAA6.gtf", 
                header = F, quote = "", sep = "\t")
colnames(gtf) <- c("ScaffoldID", "annotSoftware", "feature", "start", "end", "score", "strand", "frame", "attribute", "geneName", "description")

#read in allele freq data
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/HC_allelefreq/10-21-25_HCallelefreq/v4_HCalleleFreq")

dir <- read.table("nucallelefreq_forLD_directnmts_par1HUEXSTACdata.tsv", header = TRUE, sep = "\t") %>%
  mutate(class = "direct_n-mt")
indir <- read.table("nucallelefreq_forLD_indirectnmts_par1HUEXSTACdata.tsv", header = TRUE, sep = "\t") %>%
  mutate(class = "indirect_n-mt")
non <- read.table("nucallelefreq_forLD_allnonnmts_par1HUEXSTACdata.tsv", header = TRUE, sep = "\t") %>%
  mutate(class = "non_n-mt")

#join data and filter to make sure only well supported AIMs are present
dat <- rbind(dir, indir, non) %>%
  filter(par2samples > 50)

#read in AIM info
snp <- read.table(paste0("../../", population, "_SNPinfo.tsv"), header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) 

snp$ID <- gsub("-", ".", snp$ID)

#join AIMs to allele frequency data
joined_SNPs <- right_join(snp, dat %>% mutate(ID = ScaffoldID), 
                          by = "ID")

#join AIMs, allele freqs, to substitution data, keeping only the subs in joined_SNPs
goodDat <- left_join(joined_SNPs, subs, by = c("SNP_scaffold", "locus"))

#read in fasta file of CDS transcripts (nucleotide)
cds <- read.fasta(paste0(path, "X_birchmanni_yDAA6_CDS.fasta"),
                  seqtype = "DNA")

#identify non-synonymous and synonymous variants

#iterate through each gene that has AIMs 
for (geneID in 1:length(unique(goodDat$attribute))) {
  #SAK testing
  # if (geneID < 6000) {
  #   next
  # } else if (geneID > 7000) {
  #   break
  # }
  
  #get attribute ID (transcript ID) and the CDS annotation for that gene
  geneOfInterest <- unique(goodDat$attribute)[geneID]
  gene <- gtf %>% filter(grepl(geneOfInterest, attribute))
  
  #get all of the AIMs that are in the given gene
  gene_snpID <- snp %>% 
    filter(attribute == geneOfInterest) %>%
    select(SNP_scaffold, locus, ID) 
  
  #iterate through each AIM for a given gene
  for (k in 1:length(gene_snpID$ID)) {
    #SAK testing
    #k <- 1
    # if (geneOfInterest == "g11454.t1") { #} | geneOfInterest == "g11462.t1") {
    #   next
    # }
    
    #determine which CDS the AIM is in
    cdsIndex <- NA
    for (i in 1:length(gene$feature)) {
      if (gene_snpID$locus[k] >= gene$start[i] & gene_snpID$locus[k] <= gene$end[i]) {
        print(paste("the AIM is in CDS number", as.character(i)))
        cdsIndex <- i
        if (cdsIndex == 1) { #if the AIM is in the start codon (if +, stop if -), change the cdsIndex to the 1st CDS, which includes the start codon (no need to do this for stop codon bc the last CDS will be iterated through before stop codon (if +, start if -_))
          cdsIndex <- 2
        }
        break #breaks for loop through CDS rows once a match is found. assumes AIM only falls in one exon, which it should
      } else {
        print(paste("the AIM is not in cds number", as.character(i)))
      }
    }
    
    #if the AIM was not found in a CDS, note it in goodDat and go to next AIM
    if(is.na(cdsIndex)) {
      index1 <- which(goodDat$ID == gene_snpID$ID[k])
      goodDat$subType[index1] <- NA
      next
    }
    
    #calculate from start codon what bp in gene it is
    startCodon <- gene %>% filter(feature == "start_codon")
    cdsPositions <- gene %>% filter(feature == "CDS")
    if (length(unique(gene$strand)) > 1) { #checks that all rows have the same strandedness
      stop(paste("FATAL ERROR: not all exons have the same strandedness. This message was caused by the AIM"), as.character(gene_snpID$ID[k]))
    }
    #if gene is on the - strand, start codon is last CDS and bp loci decrease instead of increasing
    if (gene$strand[2] == "-") { 
      startCodonStart <- startCodon$end #so the end bp is the beginning of the gene
      CDSBeg <- gene[cdsIndex, "end"] #get the beginning bp of the CDS the AIM was found in
      
      #calculate the length of the previous exons for a - strand gene
      prevExonLength <- 0
      if (cdsIndex == (length(gene$feature) - 1) ) { #if the AIM is in the 1st CDS (2nd to last row for a - encoded gene)
        #the # bps from the start codon to the substitution = the bp of the start codon - the locus of the AIM + 1 (for difference vs location)
        diffFromStart <- (CDSBeg - gene_snpID$locus[k]) + 1 
      } else {
        for (j in (cdsIndex + 1):(length(gene$feature) - 1)) { #need to remove start / stop codon -- make sure that it always goes stop, cds, start if - strand
          #calculate summed length of all previous exons, which is the sum of the difference in start & stop for each CDS
          prevExonLength <- prevExonLength + gene[j,"end"] - gene[j,"start"] + 1 #+1 bc of counting vs position number
        }
        #diffFromStart will now be the total of the previous exon length plus the difference from the CDS start and the AIM locus
        diffFromStart <- prevExonLength + (CDSBeg - gene_snpID$locus[k]) + 1 
      }
      
      #if gene is on + strand, start codon is first CDS and bp loci increase
    } else if (gene$strand[2] == "+") {
      startCodonStart <- startCodon$start #get first bp of gene (1st bp of start codon)
      CDSBeg <- gene[cdsIndex, "start"] #get the first bp of the CDS that the AIM is in
      
      prevExonLength <- 0
      if (cdsIndex == 2) { #if the AIM is in the 1st CDS (2nd row for a gene)
        diffFromStart <- (gene_snpID$locus[k] - CDSBeg) + 1 #then the diff from start is just the difference (plus one) between CDS start and the AIM locus
      } else { 
        for (j in 2:(cdsIndex - 1)) { #excluding stop codon
          #calculate previous exon length, which is the sum of the difference in each previous CDS's beginning and end
          prevExonLength <- prevExonLength + gene[j,"end"] - gene[j,"start"] + 1 #+1 bc of the weird subtraction vs position thing
        }
        #then calculate difference from gene start by adding the total of the previous exons' lengths to the difference in the CDS beginning and the AIM locus
        diffFromStart <- prevExonLength + (gene_snpID$locus[k] - CDSBeg) + 1
      }
      
    }
    
    #figure out what codon position the AIM is based on the difference from start
    #and pull the codon including the AIM from the transcript fasta
    if (diffFromStart %% 3 == 1) {
      aaPos <- 1
      aimFasta <- cds[[geneOfInterest]][diffFromStart:(diffFromStart + 2)]
    } else if (diffFromStart %% 3 == 2) {
      aaPos <- 2
      aimFasta <- cds[[geneOfInterest]][(diffFromStart - 1):(diffFromStart + 1)]
    } else if (diffFromStart %% 3 == 0) {
      aaPos <- 3
      aimFasta <- cds[[geneOfInterest]][(diffFromStart - 2):diffFromStart]
    }
    
    xbir <- ""
    xcor <- ""
    
    #find the row number for the row that has the AIM ID in it
    indexGoodDat <- which(goodDat$ID == gene_snpID$ID[k])
    
    #translate the codon for xbir (the genome read in)
    xbir <- seqinr::translate(aimFasta, frame = 0)
    
    #check that the seqs are what they should be based on data in goodDat, so pull the substitution info in goodDat for this AIM
    xbirSub <- toupper(aimFasta[aaPos])
    xcorSub <- toupper(aimFasta[aaPos])
    if (gene$strand[2] == "-") {
      checkSub <- toupper(comp(goodDat[which(goodDat$ID == gene_snpID$ID[k]), "Xbir"])) #must r/c if on minus strand
    } else if (gene$strand[2] == "+") {
      checkSub <- toupper(goodDat[which(goodDat$ID == gene_snpID$ID[k]), "Xbir"])
    }
    
    #if the searched for AIM cannot be found in the data
    if (length(checkSub) == 0) {
      goodDat$xbirAA[indexGoodDat] <- "unk"
      goodDat$xcorAA[indexGoodDat] <- "unk"
      goodDat$AAresid[indexGoodDat] <- NA
      goodDat$codonPos[indexGoodDat] <- NA
      goodDat$subType[indexGoodDat] <- "unk"
      next #and then move to the next AIM
    } else if (is.na(checkSub)) {     #if there is not substitution data for the AIM found in a CDS (and is in the allele freq data), it is missing from the substitution data, in which case, you should note that in an easy-to-find-way later
      goodDat$xbirAA[indexGoodDat] <- "unk"
      goodDat$xcorAA[indexGoodDat] <- "unk"
      goodDat$AAresid[indexGoodDat] <- NA
      goodDat$codonPos[indexGoodDat] <- NA
      goodDat$subType[indexGoodDat] <- "unk"
      next #and then move to the next AIM
    }
    
    #otherwise, check that your xbir nt for the substitution is correct
    if (xbirSub != checkSub) {
      #after testing & chatting with molly, these will be treated also as "unk", as this is likely error in using non-filtered substitution data
      goodDat$xbirAA[indexGoodDat] <- "unk"
      goodDat$xcorAA[indexGoodDat] <- "unk"
      goodDat$AAresid[indexGoodDat] <- NA
      goodDat$codonPos[indexGoodDat] <- NA
      goodDat$subType[indexGoodDat] <- "unk"
      next #and then move to the next AIM
      
      # stop(paste("FATAL ERROR: it appears you do not have the right gene position. AIM", 
      #            as.character(gene_snpID$ID[k]), 
      #            "in",
      #            as.character(geneOfInterest), 
      #            "should have", 
      #            as.character(goodDat[which(goodDat$ID == gene_snpID$ID[k]), "Xbir"]),
      #            "but you have pulled",
      #            as.character(xbirSub)))
    }
    
    #substitute xcor seq
    xcorseq <- aimFasta
    if (gene$strand[2] == "-") {
      xcorseq[aaPos] <- comp(goodDat$Xcor[indexGoodDat]) #must r/c if on minus strand
    } else if (gene$strand[2] == "+") {
      xcorseq[aaPos] <- goodDat$Xcor[indexGoodDat]
    }
    
    
    xcor <- seqinr::translate(xcorseq, frame = 0)
    
    if (xbir == xcor) {
      goodDat$subType[indexGoodDat] <- "S"
    } else {
      goodDat$subType[indexGoodDat] <- "NS"
    }
    
    goodDat$xbirAA[indexGoodDat] <- xbir
    goodDat$xcorAA[indexGoodDat] <- xcor
    goodDat$AAresid[indexGoodDat] <- ceiling(diffFromStart / 3)
    goodDat$codonPos[indexGoodDat] <- aaPos
    
  } #end loop through AIMs per gene
  
} #end loop through goodDat genes

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")
population <- "HUEXSTAC"
goodDat <- read.table(paste0(population, "goodDat-w-pcor.tsv"),
                       sep = "\t", header = T, quote = "")

#create data frame with # AIMs/gene, # syn subs / gene, # non-syn subs / gene
aimCount <- goodDat %>%
  group_by(attribute) %>%
  summarise(numAIMs.gene = n_distinct(ID))

classInfo <- goodDat %>%
  distinct(attribute, class.x) %>%
  select(attribute, class.x)

avgLD <- goodDat %>%
  group_by(attribute, subType) %>%
  summarise(avgAFreq = mean(par2allele_freq),
            numFishSamples = sum(par2samples) / n(),
            numAIMs.gene.sub = n_distinct(ID)) %>%
  left_join(aimCount, by = "attribute") %>%
  left_join(classInfo, by = "attribute") #add class

goodDat$subType[is.na(goodDat$subType)] <- "non-cds" 
goodDat %>% #goodDat has all sub types, including AIMs outside of CDS
  group_by(class.x, subType) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(pct_within_class = n / sum(n))


#calculate CDS length for each gene
geneInfo <- data.frame(attribute = names(cds), 
                       geneLength = NA)
for (i in 1:length(geneInfo$attribute)) {
  geneInfo$geneLength[i] = length(cds[[i]]) 
} 

joined_avg <- left_join(avgLD, geneInfo, by = "attribute")

joined_avg$subType[is.na(joined_avg$subType)] <- "non-cds" 

joined_avg$subType <- factor(joined_avg$subType,
                             levels = c("unk", "non-cds", "S", "NS"))

#write data to output

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

write.table(goodDat, paste0(population, "goodDat-w-pcor.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

write.table(joined_avg, paste0(population, "joined-avg-w-pcor.tsv"),
            sep = "\t", row.names = F, col.names = T, quote = F)

#read and join in LD vals for D' and r^2

# #D'
# setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/d-prime/", population, "/"))
# 
# Ddir <- read.table(paste0(population, "directnmt_avgd-prime.tsv"), header = TRUE, sep = "\t") %>% 
#   mutate(class = "direct_n-mt")
# Dindir <- read.table(paste0(population, "indirectnmt_avgd-prime.tsv"), header = TRUE, sep = "\t") %>% 
#   mutate(class = "indirect_n-mt")
# Dnon <- read.table(paste0(population, "nonnmt_avgd-prime.tsv"), header = TRUE, sep = "\t") %>% 
#   mutate(class = "non_n-mt")
# 
# Ddat <- rbind(Ddir, Dindir, Dnon)
# 
# #r^2
# setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/r-squared/", population, "/"))
# 
# rdir <- read.table(paste0(population, "directnmt_avgr-squared.tsv"), header = TRUE, sep = "\t") %>% 
#   mutate(class = "direct_n-mt")
# rindir <- read.table(paste0(population, "indirectnmt_avgr-squared.tsv"), header = TRUE, sep = "\t") %>% 
#   mutate(class = "indirect_n-mt")
# rnon <- read.table(paste0(population, "nonnmt_avgr-squared.tsv"), header = TRUE, sep = "\t") %>% 
#   mutate(class = "non_n-mt")
# 
# rdat <- rbind(rdir, rindir, rnon)

# rjoined_SNPs <- right_join(snp, rdat %>% mutate(ID = ScaffoldID) %>% select(-locus), 
#                            by = "ID")


# Make non-syn dataframe w/ all LD vals

nonSynGenes <- joined_avg %>%
  filter(subType == "NS") %>%
  select(attribute)

synSubs <- joined_avg %>%
  anti_join(nonSynGenes, by = "attribute") %>%
  filter(subType == "S") %>%
  mutate(NS = 0) 
synSubs$numAIMs.gene.sub <- 0

quantDat <- joined_avg %>%
  filter(subType == "NS") %>%
  mutate(NS = numAIMs.gene.sub / geneLength) 
quantDat <- rbind(quantDat, synSubs) 
quantDat$class.x <- factor(quantDat$class.x,
                           levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt")) ###NEED TO CHANGE TO MAKE NON-N-MT THE REF

#write final output df with all LD values and non-syn calculations
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

write.table(quantDat, paste0(population, "quantDat-all-LD.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)


