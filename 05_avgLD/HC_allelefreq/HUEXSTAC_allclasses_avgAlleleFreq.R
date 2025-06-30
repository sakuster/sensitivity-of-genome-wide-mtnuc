### Shady Kuster
### 13 February 2024
### Made for the purpose of avg allele frequency across genes

#load libraries
library(dplyr)
library(ggridges)
library(dplyr)
library(ggplot2)


#change me!!
population <- "CHAF"
species <- "X. birchmanni" #remember to change par1/2 for input files depending on what species this is!
speciesShort <- "Xbir"
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/HC_allelefreq/CHAF/CHAF_allelefreq")
#setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/HC_allelefreq/", speciesShort, "AlleleFreq"))


#read in the snp file
snp <- read.table("../CHAF_SNPinfo.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) 
## this only has 738 direct nmts... what's up with that

snp$ID <- gsub("-", ".", snp$ID)

#read in separate allele freq files
direct_fn <- "nucallelefreq_forLD_directnmts_par1CHAFdata.tsv"
indir_fn <- "nucallelefreq_forLD_indirectnmts_par1CHAFdata.tsv"
non_fn <- "nucallelefreq_forLD_allnonnmts_par1CHAFdata.tsv"
output_fn <- paste0(population, "-averageA-freq_prettyPlot.pdf")


dir <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "direct_n-mt", ID = ScaffoldID)
indir <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "indirect_n-mt", ID = ScaffoldID)
non <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "non-n-mt", ID = ScaffoldID)

ld <- rbind(dir, indir, non)

#merge based off of Scaffold ID & locus
joined_SNPs <- right_join(snp, ld, by = c("ID", "class"))

#average vals based off of gene for each class
dirAvgLD <- joined_SNPs %>%
  filter(class == "direct_n-mt") %>%
  group_by(attribute) %>%
  summarize(avgAFreq = mean(allele_freq))

indirAvgLD <- joined_SNPs %>%
  filter(class == "indirect_n-mt") %>%
  group_by(attribute) %>%
  summarize(avgAFreq = mean(allele_freq))

nonAvgLD <- joined_SNPs %>%
  filter(class == "non-n-mt") %>%
  group_by(attribute) %>%
  summarize(avgAFreq = mean(allele_freq))

avgLD <- joined_SNPs %>%
  group_by(attribute, class) %>%
  summarize(avgAFreq = mean(allele_freq))

#export avg LD tsvs
write.table(dirAvgLD, paste0(population, "_avg", speciesShort, "a-freq_directnmt.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)
write.table(indirAvgLD, paste0(population, "_avg", speciesShort, "a-freq_indirectnmt.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)
write.table(nonAvgLD, paste0(population, "_avg", speciesShort, "a-freq_nonnmt.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

write.table(avgLD, paste0(population, "_avg", speciesShort, "a-freq_allclasses.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#create plot??
p1 <- ggplot(dirAvgLD, aes(x = attribute, y = avgAFreq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  labs(title = paste0(species, " Interacting N-mt Allele Frequency"), x = "Allele Frequency")
p1

ggsave(paste0("dirn-mt_avg", speciesShort, "a-freq.pdf"), p1, width = 7, height = 4)

p2 <- ggplot(indirAvgLD, aes(x = attribute, y = avgAFreq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  labs(title = paste0(species, " Non-interacting N-mt Allele Frequency"), x = "Allele Frequency")
p2

ggsave(paste0("indirn-mt_avg", speciesShort, "a-freq.pdf"), p2, width = 7, height = 4)

p3 <- ggplot(nonAvgLD, aes(x = attribute, y = avgAFreq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  labs(title = paste0(species, " Non-N-mt Allele Frequency"), x = "Allele Frequency")
p3

ggsave(paste0("non-n-mt_avg", speciesShort, "a-freq.pdf"), p3, width = 7, height = 4)

p4 <- ggplot(avgLD, aes(x = attribute, y = avgAFreq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  labs(title = paste0(species, " Allele Frequency"), x = "Allele Frequency")
p4

ggsave(paste0("allclasses_avg", speciesShort, "a-freq.pdf"), p4, width = 7, height = 4)


#get top 1% transcripts -- all classes
numTranscripts <- length(avgLD$attribute)*0.01

topVals <- as.data.frame(head(avgLD[order(-avgLD$avgAFreq),], n = numTranscripts))

nomTopVals <- right_join(joined_SNPs, topVals, by = "attribute") #can maintain order with inner_join() but doesn't work

order.scores <- order(-nomTopVals$avgAFreq) 
top1 <- nomTopVals[order.scores,]

write.table(top1, paste0(species, "_", population,"allclasses_topOnePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts -- all classes
numTranscripts5 <- round(length(avgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(avgLD[order(-avgLD$avgAFreq),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgAFreq) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, paste0(species, "_", population,"allclasses_topFivePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)


#get top 1% transcripts -- interacting n-mt
numTranscriptsDir <- ceiling(length(dirAvgLD$attribute)*0.01)

topValsDir <- as.data.frame(head(dirAvgLD[order(-dirAvgLD$avgAFreq),], n = numTranscriptsDir))

nomTopValsDir <- right_join(joined_SNPs, topValsDir, by = "attribute") #can maintain order with inner_join() but doesn't work

order.scoresDir <- order(-nomTopValsDir$avgAFreq) 
top1Dir <- nomTopValsDir[order.scoresDir,]

write.table(top1Dir, paste0("directn-mt_", species, "_", population,"_topOnePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 1% transcripts -- non-interacting n-mt
numTranscriptsIndir <- length(indirAvgLD$attribute)*0.01

topValsIndir <- as.data.frame(head(indirAvgLD[order(-indirAvgLD$avgAFreq),], n = numTranscriptsIndir))

nomTopValsIndir <- right_join(joined_SNPs, topValsIndir, by = "attribute") #can maintain order with inner_join() but doesn't work

order.scoresIndir <- order(-nomTopValsIndir$avgAFreq) 
top1Indir <- nomTopValsIndir[order.scoresIndir,]

write.table(top1Indir, paste0("indirectn-mt_", species, "_", population,"_topOnePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 1% transcripts -- non-n-mt
numTranscriptsNon <- length(nonAvgLD$attribute)*0.01

topValsNon <- as.data.frame(head(nonAvgLD[order(-nonAvgLD$avgAFreq),], n = numTranscriptsNon))

nomTopValsNon <- right_join(joined_SNPs, topValsNon, by = "attribute") #can maintain order with inner_join() but doesn't work

order.scoresNon <- order(-nomTopValsNon$avgAFreq) 
top1Non <- nomTopValsNon[order.scoresNon,]

write.table(top1Non, paste0("non-n-mt_", species, "_", population,"_topOnePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts -- interacting n-mt
numTranscripts5 <- round(length(dirAvgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(dirAvgLD[order(-dirAvgLD$avgAFreq),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgAFreq) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, paste0("directn-mt_", species, "_", population,"_topFivePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts -- non-interacting n-mt
numTranscripts5 <- round(length(indirAvgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(indirAvgLD[order(-indirAvgLD$avgAFreq),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgAFreq) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, paste0("indirectn-mt_", species, "_", population,"_topFivePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts -- non-n-mt
numTranscripts5 <- round(length(nonAvgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(nonAvgLD[order(-nonAvgLD$avgAFreq),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgAFreq) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, paste0("nonn-mt_", species, "_", population,"_topFivePercentAFreq.tsv"), row.names = F, sep = '\t',
            col.names = T, quote = F)

