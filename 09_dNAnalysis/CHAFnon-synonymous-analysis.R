##dN/dS analysis -- for CHAF

library(ape)
library(seqinr)
library(ggridges)
library(dplyr)
library(ggplot2)
library(car)
library(ggrepel)

###DEPRECATED -- DO NOT USE!@
# SEE:
# NON-SYNONYMOUS-ANALYSIS.R
# NON-SYN-STAT-ANALYSIS.R
# NON-SYN-FIG-MAKER.R
#FOR BOTH CALL, CHAF

#do dN and dS differ in relationship to mean gene class?
#does number of dN per gene affect this relationship as well?


#so probably some kind of linear regression? 
#LD correlation coefficient ~ numNSS
#LD correlation coeff ~ subType

#read in subsitution data
path <- "/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/"
subs <- read.csv(paste0(path,"Xbirchmanni10xgenome_addmito_ancestry_informative_sites_filterF1"), 
                 header = F, quote = "", sep = "\t") %>%
  mutate(subType = "", xbirAA = "", xmalAA = "", AAresid = NA, codonPos = NA) %>%
  rename(SNP_scaffold = V1, locus = V2, Xbir = V3, Xmal = V4)

#read in gtf file
gtf <- read.csv("/Users/kusters/Documents/XiphoStartingOver_Aug2024/00_StartingDataFiles/xiphophorus_birchmanni_10x_12Sep2018_yDAA6.gtf", 
                header = F, quote = "", sep = "\t")
colnames(gtf) <- c("ScaffoldID", "annotSoftware", "feature", "start", "end", "score", "strand", "frame", "attribute", "geneName", "description")

#read in the avgLD data
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/08_partialCorrelation/averagePcor/CHAF")

#read it in
dir <- read.table("mt1_results_forLD_directnmts_par1CHAFdata.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "direct_n-mt")
indir <- read.table("mt1_results_indirectnmt_CHAF_ALL.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "indirect_n-mt")
non <- read.table("mt1_results_nonnmt_CHAF_ALL.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "non_n-mt")

dat <- rbind(dir, indir, non)

snp <- read.table("CHAF_SNPinfo.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) 

snp$ID <- gsub("-", ".", snp$ID)

joined_SNPs <- right_join(snp, dat %>% mutate(ID = ScaffoldID) %>% select(-locus), 
                          by = "ID")

# join sub df with SNP df
goodDat <- inner_join(joined_SNPs, subs, by = c("SNP_scaffold", "locus"))


# dat <- dat %>%
#   filter(numAIMs > 1) #numSamples >180)

#read in fasta file
#use to do gene length, determine syn vs non-syn
cds <- read.fasta(paste0(path, "X_birchmanni_yDAA6_CDS.fasta"),
                  seqtype = "DNA")

#identify non-synonymous and synonymous variants
for (geneID in 1:length(unique(goodDat$attribute))) {
  geneOfInterest <- unique(goodDat$attribute)[geneID]
  gene <- gtf %>% filter(grepl(geneOfInterest, attribute))
  
  gene_snpID <- snp %>% 
    filter(attribute == geneOfInterest) %>%
    select(SNP_scaffold, locus, ID) #3rd AIM in a CDS, others in introns
  
  for (k in 1:length(gene_snpID$ID)) {
    #SAK testing
    #k <- 2
    
    #determine which CDS the SNP is in -- currently only works for 1 AIM / gene
    # WHAT IF THE AIM IS IN THE START CODON?? -- shouldn't ever be the case...
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
    #keep in mind strandedness
    #make note of codon frame??
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
          prevExonLength <- prevExonLength + gene[j,"end"] - gene[j,"start"] + 1 #+1 bc of the weird subtraction vs position thing
        }
        diffFromStart <- prevExonLength + (CDSBeg - gene_snpID$locus[k]) + 1 
      }
      
    } else if (gene$strand[2] == "+") {
      startCodonStart <- startCodon$start 
      CDSBeg <- gene[cdsIndex, "start"]
      #need to do the above but for + strand seq
      
      prevExonLength <- 0
      if (cdsIndex == 2) { #if the AIM is in the 1st CDS (2nd row for a gene)
        diffFromStart <- (gene_snpID$locus[k] - CDSBeg) + 1
      } else {
        for (j in 2:(cdsIndex - 1)) { #need to remove start / stop codon -- make sure that it always goes stop, cds, start if + strand
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
      xmalseq[aaPos] <- comp(goodDat$Xmal[indexGoodDat]) #must r/c if on minus strand, assuming the file molly gave you is for genomic positions, not by gene
    } else if (gene$strand[2] == "+") {
      xmalseq[aaPos] <- goodDat$Xmal[indexGoodDat]
    }
    #i might have the xbir and xmal cols mixed up. continuing as if i do not
    
    
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


# #average LD vals based off of gene
# avgLD <- joined_SNPs %>%
#   group_by(attribute, subType) %>%
#   summarize(avgP_cor = mean(pcor),
#             numFishSamples = sum(samples) / n(),
#             numAIMs = n_distinct(ID)) 

aimCount <- goodDat %>%
  group_by(attribute) %>%
  summarise(numAIMs.gene = n_distinct(ID))

classInfo <- goodDat %>%
  distinct(attribute, class.x) %>%
  select(attribute, class.x)

avgLD <- goodDat %>%
  group_by(attribute, subType) %>%
  summarise(avgP_cor = mean(pcor), #add D', r^2 and do same thing??
            numFishSamples = sum(samples) / n(),
            numAIMs.gene.sub = n_distinct(ID)) %>%
  left_join(aimCount, by = "attribute") %>%
  left_join(classInfo, by = "attribute")

goodDat$subType[is.na(goodDat$subType)] <- "non-cds" 
goodDat %>%
  group_by(class.x, subType) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(pct_within_class = n / sum(n))


#calculate gene length for each gene
geneInfo <- data.frame(attribute = names(cds), 
                       geneLength = NA)
for (i in 1:length(geneInfo$attribute)) {
  geneInfo$geneLength[i] = length(cds[[i]]) 
} 

joined_avg <- left_join(avgLD, geneInfo, by = "attribute")

joined_avg$subType[is.na(joined_avg$subType)] <- "non-cds" 

joined_avg$subType <- factor(joined_avg$subType,
                             levels = c("non-cds", "S", "NS"))

#### ----
#make plots

# # sub type vs gene length
# p <- ggplot(joined_avg, aes(x = subType, y = log(geneLength))) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.03, aes(color = subType)) +
#   theme_bw()
# 
# ggplot_build(p)$data
# 
# # num of AIMs per gene per subType vs geneLength
# ggplot(joined_avg, aes(x = log(geneLength), y = numAIMs.gene.sub, fill = subType)) +
#   geom_point(aes(color = subType)) +
#   scale_color_manual(values = c("non-cds" = "lightgray", "S" = "#00BA38", "NS" = "#619CFF")) +
#   # geom_point(aes(color = subType), data = joined_avg %>%
#   #              filter(!is.na(subType))) +
#   #facet_grid(~subType) +
#   theme_bw() +
#   labs(title = "Are longer genes more likely to have more substitutions?") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# NAmod <- lm(numAIMs.gene.sub ~ geneLength,
#             data = joined_avg %>% filter(is.na(subType)))
# Anova(mod1, type = "III") #v significant
# 
# Smod <- lm(numAIMs.gene.sub ~ geneLength,
#            data = joined_avg %>% filter(subType == "S"))
# Anova(Smod, type = "III") #v significant
# 
# NSmod <- lm(numAIMs.gene.sub ~ geneLength,
#             data = joined_avg %>% filter(subType == "NS"))
# Anova(NSmod, type = "III") #also v signficant
# 
# allMod <- lm(numAIMs.gene.sub ~ geneLength + subType,
#              data = joined_avg)
# Anova(allMod, type = "III")
# 
# pairs(emmeans(allMod, ~ subType))
# 
# # num of AIMs per gene per subType vs geneLength but NO NON-CDS
# ggplot(joined_avg %>% filter(subType != "non-cds"), 
#        aes(x = geneLength, y = numAIMs.gene.sub, fill = subType)) +
#   geom_point(aes(color = subType)) +
#   scale_color_manual(values = c("S" = "#00BA38", "NS" = "#619CFF")) +
#   # geom_point(aes(color = subType), data = joined_avg %>%
#   #              filter(!is.na(subType))) +
#   facet_grid(~subType) +
#   theme_bw() +
#   labs(title = "Are longer genes more likely to have more substitutions?") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # num AIMS per gene vs sub type
# ggplot(joined_avg, aes(x = subType, y = log(numAIMs.gene))) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.03, aes(color = subType)) +
#   theme_bw()
# 
# # across gene classes, what does sub type look like?
# ggplot(joined_avg, aes(x = subType)) +
#   geom_histogram(stat = "count") +
#   facet_wrap(~class.x) +
#   scale_y_log10() +
#   theme_bw() +
#   labs(y = "log(count)")


twoCat <- joined_avg %>%
  group_by(attribute) %>%
  mutate(
    genericSubCat = if (any(subType == "NS")) "atLeast1NS" else "none"
  ) %>%
  distinct(attribute, .keep_all = TRUE) #g.1000.t1 has 5 NA, 1 S, and 1 NS sub. it should be classified as NS

twoCat2 <- joined_avg %>%
  filter(subType != "non-cds") %>%
  group_by(attribute) %>%
  mutate(
    genericSubCat = if (any(subType == "NS")) "atLeast1NS" else "none"
  ) %>%
  distinct(attribute, .keep_all = TRUE)




# do dir n-mt differ in num of ns subs?
ggplot(twoCat, aes(x = log(geneLength), y = numAIMs.gene)) +
  geom_point(aes(color = "no NS substitutions"),
             alpha = 0.5,
             data = twoCat %>% filter(genericSubCat == "none")) +
  geom_point(aes(color = "at least 1 NS sub"),
             alpha = 0.5,
             data = twoCat %>% filter(genericSubCat == "atLeast1NS")) +
  #scale_color_manual(values = c("non-cds" = "lightgray", "non-syn" = "#619CFF")) +
  # geom_point(aes(color = subType), data = joined_avg %>%
  #              filter(!is.na(subType))) +
  facet_grid(~class.x) +
  theme_bw() +
  labs(title = "Are longer genes more likely to have more substitutions?") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#boxplot of NS vs none grouping
ggplot(twoCat, aes(x = genericSubCat, y = log(geneLength))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, 
              aes(color = genericSubCat)) +
  facet_grid(~class.x) +
  theme_bw()


#write data to output

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

write.table(goodDat, "CHAFgoodDat-w-pcor.tsv", 
            sep = "\t", row.names = F, col.names = T, quote = F)

write.table(twoCat, "CHAFtwoCat-w-pcor.tsv",
            sep = "\t", row.names = F, col.names = T, quote = F)

write.table(joined_avg, "CHAFjoined-avg-w-pcor.tsv",
            sep = "\t", row.names = F, col.names = T, quote = F)

joined_avg <- read.table("CHAFjoined-avg-w-pcor.tsv", 
                         sep = "\t", header = T, quote = "")

#read and join in LD vals for D' and r^2

#D'
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/d-prime/CHAF/")

Ddir <- read.table("CHAFdirectnmt_avgd-prime.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "direct_n-mt")
Dindir <- read.table("CHAFindirectnmt_avgd-prime.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "indirect_n-mt")
Dnon <- read.table("CHAFnonnmt_avgd-prime.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "non_n-mt")

Ddat <- rbind(Ddir, Dindir, Dnon)

Djoined_SNPs <- right_join(snp, dat %>% mutate(ID = ScaffoldID) %>% select(-locus), 
                           by = "ID")

#r^2
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/r-squared/CHAF/")

rdir <- read.table("CHAFdirectnmt_avgr-squared.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "direct_n-mt")
rindir <- read.table("CHAFindirectnmt_avgr-squared.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "indirect_n-mt")
rnon <- read.table("CHAFnonnmt_avgr-squared.tsv", header = TRUE, sep = "\t") %>% 
  mutate(class = "non_n-mt")

rdat <- rbind(rdir, rindir, rnon)

rjoined_SNPs <- right_join(snp, dat %>% mutate(ID = ScaffoldID) %>% select(-locus), 
                           by = "ID")

#all together now
twoCat2 <- left_join(twoCat2, Ddat %>% select(attribute, avgD_prime), by = "attribute") %>%
  left_join(rdat %>% select(attribute, avgR_squared), by = "attribute")

# run statistical tests -- D PRIME

#   Mtnuc LD ~ geneClass * [at least one NS sub or only S subs] 
twoCat2$class.x <- factor(twoCat2$class.x,
                          levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt")) ###NEED TO CHANGE TO MAKE NON-N-MT THE REF)

catMod <- lm(avgD_prime ~ class.x * genericSubCat, 
             data = twoCat2)
Anova(catMod, type = "III") #for D'

twoCat2 %>%
  group_by(class.x, genericSubCat) %>%
  summarise(n = n(), 
            meanDprime = mean(avgD_prime, na.rm = T), 
            sdDprime = sd(avgD_prime, na.rm = T))

mypalette <- c("#D589B1", "#A4BDBE", "#457373")
ggplot(twoCat2 %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = genericSubCat, y = avgD_prime)) +
  geom_jitter(width = 0.2, alpha = 0.3, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  geom_violin(trim = F, aes(fill = class.x)) +
  scale_fill_manual(values = mypalette) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Substitution Category of a Gene",
       y = "D'") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Mtnuc LD ~ geneClass * [ # NS subs / CDS length] 
synSubs <- joined_avg %>%
  filter(subType == "S") %>%
  mutate(NS = 0) 
synSubs$numAIMs.gene.sub <- 0

quantDat <- joined_avg %>%
  filter(subType == "NS") %>%
  mutate(NS = numAIMs.gene.sub / geneLength) 
quantDat <- rbind(quantDat, synSubs) %>%
  left_join(Ddat %>% select(attribute, avgD_prime), by = "attribute") %>%
  left_join(rdat %>% select(attribute, avgR_squared), by = "attribute")
quantDat$class.x <- factor(quantDat$class.x,
                           levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt")) ###NEED TO CHANGE TO MAKE NON-N-MT THE REF

quantMod <- lm(avgD_prime ~ class.x * NS, #for D'
               data = quantDat)
Anova(quantMod, type = "III")
summary(quantMod)
# quantMod2 <- lm(avgD_prime ~ numAIMs.gene.sub * class.x + geneLength,
#                 data = quantDat)
# Anova(quantMod2, type = "III")
# summary(quantMod2)

quantDat %>%
  group_by(class.x) %>%
  summarise(n = n(),
            meanDprime = mean(avgD_prime, na.rm = T),
            sdDprime = sd(avgD_prime, na.rm = T),
            meanNS = mean(NS),
            sdNS = sd(NS))

#this fig is PER GENE. each dot is a gene. pcor is gene-averaged pcor
mypalette <- c("#D589B1", "#A4BDBE", "#457373")
ggplot(quantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = NS, y = avgD_prime)) +
  geom_point(alpha = 0.15, aes(color = class.x)) +
  geom_smooth(method = "lm", 
              se = FALSE, alpha = 0.2, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Non-synonymous substitutions count / CDS length (bp)",
       y = "D'") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(quantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = numAIMs.gene.sub, y = avgD_prime)) +
  geom_point(alpha = 0.15, aes(color = geneLength)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2) +
  scale_color_viridis_c() +
  #scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Number of non-synonymous substitutions in a gene",
       y = "D'") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


# R-SQUARED

#   Mtnuc LD ~ geneClass * [at least one NS sub or only S subs] 
catMod <- lm(avgR_squared ~ class.x * genericSubCat, 
             data = twoCat2)
Anova(catMod, type = "III") #for r-sq

twoCat2 %>%
  group_by(class.x, genericSubCat) %>%
  summarise(n = n(), 
            meanRsq = mean(avgR_squared, na.rm = T), 
            sdRsq = sd(avgR_squared, na.rm = T))

mypalette <- c("#D589B1", "#A4BDBE", "#457373")
ggplot(twoCat2 %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = genericSubCat, y = avgR_squared)) +
  geom_jitter(width = 0.2, alpha = 0.3, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  geom_violin(trim = F, aes(fill = class.x)) +
  scale_fill_manual(values = mypalette) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Substitution Category of a Gene",
       y = "R-squared") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Mtnuc LD ~ geneClass * [ # NS subs / CDS length] 

quantMod <- lm(avgR_squared ~ class.x * NS, #for D'
               data = quantDat)
Anova(quantMod, type = "III")
summary(quantMod)
# quantMod2 <- lm(avgR_squared ~ numAIMs.gene.sub * class.x + geneLength,
#                 data = quantDat)
# Anova(quantMod2, type = "III")
# summary(quantMod2)

quantDat %>%
  group_by(class.x) %>%
  summarise(n = n(),
            meanRsq = mean(avgR_squared, na.rm = T),
            sdRsq = sd(avgR_squared, na.rm = T),
            meanNS = mean(NS),
            sdNS = sd(NS))

#this fig is PER GENE. each dot is a gene. pcor is gene-averaged pcor
mypalette <- c("#D589B1", "#A4BDBE", "#457373")
ggplot(quantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = NS, y = avgR_squared)) +
  geom_point(alpha = 0.15, aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Non-synonymous substitutions count / CDS length (bp)",
       y = "R-squared") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(quantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = numAIMs.gene.sub, y = avgR_squared)) +
  geom_point(alpha = 0.15, aes(color = geneLength)) +
  scale_color_viridis_c() +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2) +
  #scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Number of non-synonymous substitutions in a gene",
       y = "r^2") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))






### PARTIAL CORRELATION
#   Mtnuc LD ~ geneClass * [at least one NS sub or only S subs] 
catMod <- lm(avgP_cor ~ class.x * genericSubCat, 
             data = twoCat2)
Anova(catMod, type = "III") #for pcor

twoCat2 %>%
  group_by(class.x, genericSubCat) %>%
  summarise(n = n(), 
            meanPcor = mean(avgP_cor, na.rm = T), 
            sdPcor = sd(avgP_cor, na.rm = T))

mypalette <- c("#D589B1", "#A4BDBE", "#457373")
ggplot(twoCat2 %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = genericSubCat, y = avgP_cor)) +
  geom_jitter(width = 0.2, alpha = 0.3, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  geom_violin(trim = F, aes(fill = class.x)) +
  scale_fill_manual(values = mypalette) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Substitution Category of a Gene",
       y = "Partial Correlation Score") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Mtnuc LD ~ geneClass * [ # NS subs / CDS length] 
quantMod <- lm(avgP_cor ~ class.x * NS, #for pcor
               data = quantDat)
Anova(quantMod, type = "III")
summary(quantMod)
# quantMod2 <- lm(avgP_cor ~ numAIMs.gene.sub * class.x + geneLength,
#                 data = quantDat)
# Anova(quantMod2, type = "III")
# summary(quantMod2)



quantDat %>%
  group_by(class.x) %>%
  summarise(n = n(),
            meanPcor = mean(avgP_cor, na.rm = T),
            sdPcor = sd(avgP_cor, na.rm = T),
            meanNS = mean(NS),
            sdNS = sd(NS))

#this fig is PER GENE. each dot is a gene. pcor is gene-averaged pcor
mypalette <- c("#D589B1", "#A4BDBE", "#457373")
ggplot(quantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = NS, y = avgP_cor)) +
  geom_point(alpha = 0.15, aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Non-synonymous substitutions count / CDS length (bp)",
       y = "Partial correlation score") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(quantDat%>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = numAIMs.gene.sub, y = avgP_cor)) +
  geom_point(alpha = 0.15, aes(color = geneLength)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2) +
  scale_color_viridis_c() +
  #scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Number of non-synonymous substitutions in a gene",
       y = "Partial correlation score") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))




#does adding NS actually improve the model fit?
ogModD <- lm(avgD_prime ~ class.x, data = quantDat)
quantModD <- lm(avgD_prime ~ class.x*NS, data = quantDat)
anova(ogModD, #model w/ just gene class
      quantModD) #model w/ gene class and NS interaction
#for Dprime, model NOT improved by addition of NS

ogModR <- lm(avgR_squared ~ class.x, data = quantDat)
quantModR <- lm(avgR_squared ~ class.x*NS, data = quantDat)
anova(ogModR,
      quantModR) #r^2 NOT improved by addition of NS

ogModP <- lm(avgP_cor ~ class.x, data = quantDat)
quantModP <- lm(avgP_cor ~ class.x*NS, data = quantDat)
anova(ogModP,
      quantModP) #pcor model improved by addition of NS (p = 0.048)







ggplot(twoCat, aes(x = genericSubCat, y = avgP))


# joined_avg %>%
#   filter(subType == "NS") %>%
#   mutate(NS = numAIMs.gene.sub / geneLength) %>%
#   group_by(class.x) %>%
#   summarise(n = n(),
#             meanPcor = mean(avgP_cor, na.rm = T),
#             sdPcor = sd(avgP_cor, na.rm = T),
#             meanNS = mean(NS),
#             sdNS = sd(NS))

# do permutation to confirm positive results??


# longer genes tend to have NS subs:
joined_avg %>% 
  distinct(attribute, .keep_all = T) %>% #all genes
  group_by(class.x) %>% 
  summarise(mean = mean(geneLength))

quantDat %>% 
  distinct(attribute, .keep_all = T) %>% #only genes with NS subs
  group_by(class.x) %>% 
  summarise(mean = mean(geneLength))
