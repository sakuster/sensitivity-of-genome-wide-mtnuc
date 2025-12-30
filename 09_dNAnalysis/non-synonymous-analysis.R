##dN/dS analysis

library(ape)
library(seqinr)
library(ggridges)
library(dplyr)
library(ggplot2)
library(car)
library(ggrepel)

population <- "CALL"

#read in substitution data for xmal & xbir
path <- "/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/"
subs <- read.csv(paste0(path,"Xbirchmanni10xgenome_addmito_ancestry_informative_sites_filterF1"),
                 header = F, quote = "", sep = "\t") %>%
mutate(subType = "", xbirAA = "", xmalAA = "", AAresid = NA, codonPos = NA) %>%
rename(SNP_scaffold = V1, locus = V2, Xbir = V3, Xmal = V4)

#read in gtf file
gtf <- read.csv("/Users/kusters/Documents/XiphoStartingOver_Aug2024/00_StartingDataFiles/xiphophorus_birchmanni_10x_12Sep2018_yDAA6.gtf", 
                header = F, quote = "", sep = "\t")
colnames(gtf) <- c("ScaffoldID", "annotSoftware", "feature", "start", "end", "score", "strand", "frame", "attribute", "geneName", "description")

#read in the avgLD data - can be any metric bc they are all filtered the same
setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/08_partialCorrelation/averagePcor/", population))

#read in mtnuc LD result data
dir <- read.table(paste0("mt1_results_forLD_directnmts_par1", population, "data.tsv"), header = TRUE, sep = "\t") %>% 
  mutate(class = "direct_n-mt")
indir <- read.table(paste0("mt1_results_indirectnmt_", population, "_ALL.tsv"), header = TRUE, sep = "\t") %>% 
  mutate(class = "indirect_n-mt")
non <- read.table(paste0("mt1_results_nonnmt_", population, "_ALL.tsv"), header = TRUE, sep = "\t") %>% 
  mutate(class = "non_n-mt")

dat <- rbind(dir, indir, non) %>%
  filter(samples > 50)

#read in AIM info
snp <- read.table(paste0(population, "_SNPinfo.tsv"), header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) 

snp$ID <- gsub("-", ".", snp$ID)

joined_SNPs <- right_join(snp, dat %>% mutate(ID = ScaffoldID) %>% select(-locus), 
                          by = "ID")

# join sub df with SNP df
goodDat <- inner_join(joined_SNPs, subs, by = c("SNP_scaffold", "locus"))

#read in fasta file of CDS transcripts (nucleotide)
cds <- read.fasta(paste0(path, "X_birchmanni_yDAA6_CDS.fasta"),
                  seqtype = "DNA")

#identify non-synonymous and synonymous variants

#iterate through each gene
for (geneID in 1:length(unique(goodDat$attribute))) {
  geneOfInterest <- unique(goodDat$attribute)[geneID]
  gene <- gtf %>% filter(grepl(geneOfInterest, attribute))
  
  gene_snpID <- snp %>% 
    filter(attribute == geneOfInterest) %>%
    select(SNP_scaffold, locus, ID) #3rd AIM in a CDS, others in introns
  
  #iterate through each AIM for a given gene
  for (k in 1:length(gene_snpID$ID)) {
    
    #determine which CDS the AIM is in
    cdsIndex <- NA
    for (i in 1:length(gene$feature)) {
      if (gene_snpID$locus[k] >= gene$start[i] & gene_snpID$locus[k] <= gene$end[i]) {
        print(paste("the AIM is in CDS number", as.character(i)))
        cdsIndex <- i
        if (cdsIndex == 1) {
          cdsIndex <- 2
        }
        break #breaks for loop through CDS rows once a match is found. assumes AIM only falls in one exon, which it should
      } else {
        print(paste("the AIM is not in cds number", as.character(i)))
      }
    }
    
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
    if (gene$strand[2] == "-") { 
      startCodonStart <- startCodon$end 
      CDSBeg <- gene[cdsIndex, "end"]
      prevExonLength <- 0
      if (cdsIndex == (length(gene$feature) - 1) ) { #if the AIM is in the 1st CDS (2nd to last row for a - encoded gene)
        diffFromStart <- (CDSBeg - gene_snpID$locus[k]) + 1 
      } else {
        for (j in (cdsIndex + 1):(length(gene$feature) - 1)) { #need to remove start / stop codon -- make sure that it always goes stop, cds, start if - strand
          prevExonLength <- prevExonLength + gene[j,"end"] - gene[j,"start"] + 1 #+1 bc of counting vs position number
        }
        diffFromStart <- prevExonLength + (CDSBeg - gene_snpID$locus[k]) + 1 
      }
      
    } else if (gene$strand[2] == "+") {
      startCodonStart <- startCodon$start 
      CDSBeg <- gene[cdsIndex, "start"]
      
      prevExonLength <- 0
      if (cdsIndex == 2) { #if the AIM is in the 1st CDS (2nd row for a gene)
        diffFromStart <- (gene_snpID$locus[k] - CDSBeg) + 1
      } else {
        for (j in 2:(cdsIndex - 1)) { 
          prevExonLength <- prevExonLength + gene[j,"end"] - gene[j,"start"] + 1 #+1 bc of the weird subtraction vs position thing
        }
        diffFromStart <- prevExonLength + (gene_snpID$locus[k] - CDSBeg) + 1
      }
      
    }
    
    #figure out what codon position the AIM is
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
    xmal <- ""
    
    indexGoodDat <- which(goodDat$ID == gene_snpID$ID[k])
    xbir <- seqinr::translate(aimFasta, frame = 0)
    
    #check that the seqs are what they should be
    xbirSub <- toupper(aimFasta[aaPos])
    xmalSub <- toupper(aimFasta[aaPos])
    if (gene$strand[2] == "-") {
      checkSub <- toupper(comp(goodDat[which(goodDat$ID == gene_snpID$ID[k]), "Xbir"])) #must r/c if on minus strand
    } else if (gene$strand[2] == "+") {
      checkSub <- toupper(goodDat[which(goodDat$ID == gene_snpID$ID[k]), "Xbir"])
    }
    
    #if the searched for AIM cannot be found in the data
    if (length(checkSub) == 0) {
      goodDat$xbirAA[indexGoodDat] <- "unk"
      goodDat$xmalAA[indexGoodDat] <- "unk"
      goodDat$AAresid[indexGoodDat] <- NA
      goodDat$codonPos[indexGoodDat] <- NA
      goodDat$subType[indexGoodDat] <- "unk"
      next #and then move to the next AIM
    } else if (is.na(checkSub)) {     #if there is not substitution data for the AIM found in a CDS (and is in the allele freq data), it is missing from the substitution data, in which case, you should note that in an easy-to-find-way later
      goodDat$xbirAA[indexGoodDat] <- "unk"
      goodDat$xmalAA[indexGoodDat] <- "unk"
      goodDat$AAresid[indexGoodDat] <- NA
      goodDat$codonPos[indexGoodDat] <- NA
      goodDat$subType[indexGoodDat] <- "unk"
      next #and then move to the next AIM
    }
    
    
    if (xbirSub != checkSub) {
      stop(paste("FATAL ERROR: it appears you do not have the right gene position. AIM", 
                 as.character(gene_snpID$ID[k]), 
                 "in",
                 as.character(geneOfInterest), 
                 "should have", 
                 as.character(goodDat[which(goodDat$ID == gene_snpID$ID[k]), "Xbir"]),
                 "but you have pulled",
                 as.character(xbirSub)))
    }
    
    #substitute xmal seq
    xmalseq <- aimFasta
    if (gene$strand[2] == "-") {
      xmalseq[aaPos] <- comp(goodDat$Xmal[indexGoodDat]) #must r/c if on minus strand
    } else if (gene$strand[2] == "+") {
      xmalseq[aaPos] <- goodDat$Xmal[indexGoodDat]
    }
    
    
    xmal <- seqinr::translate(xmalseq, frame = 0)
    
    if (xbir == xmal) {
      goodDat$subType[indexGoodDat] <- "S"
    } else {
      goodDat$subType[indexGoodDat] <- "NS"
    }
    
    goodDat$xbirAA[indexGoodDat] <- xbir
    goodDat$xmalAA[indexGoodDat] <- xmal
    goodDat$AAresid[indexGoodDat] <- ceiling(diffFromStart / 3)
    goodDat$codonPos[indexGoodDat] <- aaPos
    
  } #end loop through AIMs per gene
  
} #end loop through goodDat genes

# setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")
# population <- "CHAF"
# goodDat <- read.table(paste0(population, "goodDat-w-pcor.tsv"), 
#                        sep = "\t", header = T, quote = "")


#create data frame with # AIMs/gene, # syn subs / gene, # non-syn subs / gene
aimCount <- goodDat %>%
  group_by(attribute) %>%
  summarise(numAIMs.gene = n_distinct(ID))

classInfo <- goodDat %>%
  distinct(attribute, class.x) %>%
  select(attribute, class.x)

avgLD <- goodDat %>%
  group_by(attribute, subType) %>% ##must group by subType or else this data is lost
  summarise(avgP_cor = mean(pcor),
            numFishSamples = sum(samples) / n(),
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
                             levels = c("non-cds", "S", "NS"))

#write data to output

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

write.table(goodDat, paste0(population, "goodDat-w-pcor.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

write.table(joined_avg, paste0(population, "joined-avg-w-pcor.tsv"),
            sep = "\t", row.names = F, col.names = T, quote = F)

#read and join in LD vals for D' and r^2

### CHANGE TO  filter # samples to > 50 

#D' -- now all dat
setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/d-prime/", population, "/"))

Ddir <- read.table(paste0(population, "directnmt_avgd-prime.tsv"), header = TRUE, sep = "\t") %>% 
  mutate(class = "direct_n-mt")
Dindir <- read.table(paste0(population, "indirectnmt_avgd-prime.tsv"), header = TRUE, sep = "\t") %>% 
  mutate(class = "indirect_n-mt")
Dnon <- read.table(paste0(population, "nonnmt_avgd-prime.tsv"), header = TRUE, sep = "\t") %>% 
  mutate(class = "non_n-mt")

Ddat <- rbind(Ddir, Dindir, Dnon)

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
##add a catch if the gene has AIMs w/ diff subTypes 
#aka if 1 AIM is S and another is NS, how will you handle that?
##currently are being avgd separately, which is not what needs to happen!

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
quantDat <- rbind(quantDat, synSubs) %>%
  left_join(Ddat %>% select(attribute, avgD_prime, avgR_squared), by = "attribute")# %>%
  #left_join(rdat %>% select(attribute, avgR_squared), by = "attribute")
quantDat$class.x <- factor(quantDat$class.x,
                           levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt")) ###NEED TO CHANGE TO MAKE NON-N-MT THE REF

#write final output df with all LD values and non-syn calculations
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

write.table(quantDat, paste0(population, "quantDat-all-LD.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)


