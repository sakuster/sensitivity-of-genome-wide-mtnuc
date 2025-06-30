### Shady Kuster
### 9 February 2024
### Made for the purpose of avg LD values (D' or r^2) across genes

#load libraries
library(dplyr)

library(ggridges)
library(dplyr)
library(ggplot2)

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/r-squared/CHAF")

#read in the snp file
snp <- read.table("CHAF_SNPinfo.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) %>%
  filter(class == "indirect_n-mt")

snp$ID <- gsub("-", ".", snp$ID)

#read in the LD file
ld <- read.table("mt1_results_indirectnmt_CHAF_ALL.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = ScaffoldID)

#merge based off of Scaffold ID & locus
joined_SNPs <- right_join(snp, ld, by = "ID")

#average LD vals based off of gene
avgLD <- joined_SNPs %>%
  group_by(attribute) %>%
  summarize(avgR_squared = mean(r_squared))

#export avg LD tsv
write.table(avgLD, "CHAFindirectnmt_avgr-squared.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

#create plot??
p1 <- ggplot(avgLD, aes(x = attribute, y = avgR_squared)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) +
  labs(title = "CHAF Gene Average R Squared", x = "Gene")

ggsave("CHAFindirectnmt_avgr-squared.pdf", p1, width = 7, height = 4)

#get top 1% transcripts
numTranscripts <- length(avgLD$attribute)*0.01

topVals <- as.data.frame(head(avgLD[order(-avgLD$avgR_squared),], n = numTranscripts))

nomTopVals <- right_join(joined_SNPs, topVals, by = "attribute") #can maintain order with inner_join() but doesn't work
order.scores <- order(-nomTopVals$avgR_squared) 
top1 <- nomTopVals[order.scores,]

write.table(top1, "CHAFindirectnmt_topOnePercentRSquared.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts
numTranscripts5 <- round(length(avgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(avgLD[order(-avgLD$avgR_squared),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgR_squared) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, "CHAFindirectnmt_topFivePercentRSquared.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)
