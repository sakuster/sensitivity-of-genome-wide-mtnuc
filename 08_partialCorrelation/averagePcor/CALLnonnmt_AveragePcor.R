### Shady Kuster
### 9 February 2024
### Made for the purpose of avg LD values (D' or r^2) across genes

#load libraries
library(dplyr)
library(ggridges)
library(dplyr)
library(ggplot2)

#read in the snp file
snp <- read.table("CALL_SNPinfo.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) %>%
  filter(class == "non-n-mt")

snp$ID <- gsub("-", ".", snp$ID)

#read in the LD file
ld <- read.table("mt1_results_nonnmt_CALL_ALL.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = ScaffoldID)

#merge based off of Scaffold ID & locus
joined_SNPs <- right_join(snp, ld, by = "ID")

#average LD vals based off of gene
avgLD <- joined_SNPs %>%
  group_by(attribute) %>%
  summarize(avgP_cor = mean(pcor),
            numSamples = sum(samples) / n(),
            numAIMs = n_distinct(ID))

#export avg LD tsv
write.table(avgLD, "CALLnonnmt_avgpcor.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

#create plot??
p1 <- ggplot(avgLD, aes(x = attribute, y = avgP_cor)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  labs(title = "CALL Gene Average R Squared", x = "Gene")

ggsave("CALLnonnmt_avgpcor.pdf", p1, width = 7, height = 4)

#get top 1% transcripts
numTranscripts <- length(avgLD$attribute)*0.01

topVals <- as.data.frame(head(avgLD[order(-avgLD$avgP_cor),], n = numTranscripts))

nomTopVals <- right_join(joined_SNPs, topVals, by = "attribute") #can maintain order with inner_join() but doesn't work
order.scores <- order(-nomTopVals$avgP_cor) 
top1 <- nomTopVals[order.scores,]

write.table(top1, "CALLnonnmt_topOnePercentPcor.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts
numTranscripts5 <- round(length(avgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(avgLD[order(-avgLD$avgP_cor),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgP_cor) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, "CALLnonnmt_topFivePercentPcor.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

