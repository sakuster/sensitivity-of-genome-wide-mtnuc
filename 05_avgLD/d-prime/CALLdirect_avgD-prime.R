### Shady Kuster
### 9 February 2024
### Made for the purpose of avg LD values (D' or r^2) across genes

#load libraries
library(dplyr)

library(ggridges)
library(dplyr)
library(ggplot2)

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/d-prime/CALL")

#read in the snp file
snp <- read.table("CALL_SNPinfo.tsv", header = T, sep = '\t', quote = "") %>%
  mutate(ID = paste(SNP_scaffold, locus, sep = '.')) %>%
  filter(class == "direct_n-mt")
## this only has 738 direct nmts... what's up with that
  
snp$ID <- gsub("-", ".", snp$ID)

#read in the LD file
ld <- read.table("mt1_results_forLD_directnmts_par1CALLdata.tsv", header = T, sep = '\t') %>%
  mutate(ID = ScaffoldID)

#merge based off of Scaffold ID & locus
joined_SNPs <- right_join(snp, ld, by = "ID")
## this, per above, also is only 738 (the LD is 991)

#average LD vals based off of gene
avgLD <- joined_SNPs %>%
  group_by(attribute) %>%
  summarize(avgD_prime = mean(D_prime))

#export avg LD tsv
write.table(avgLD, "CALLdirectnmt_avgd-prime.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

#create plot??
p1 <- ggplot(avgLD, aes(x = attribute, y = avgD_prime)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  labs(title = "CALL Gene Average D Prime", x = "D'")

ggsave("CALLdirectnmt_avgd-prime.pdf", p1, width = 7, height = 4)

#get top 1% transcripts
numTranscripts <- length(avgLD$attribute)*0.01

topVals <- as.data.frame(head(avgLD[order(-avgLD$avgD_prime),], n = numTranscripts))

nomTopVals <- right_join(joined_SNPs, topVals, by = "attribute") #can maintain order with inner_join() but doesn't work

order.scores <- order(-nomTopVals$avgD_prime) 
top1 <- nomTopVals[order.scores,]

write.table(top1, "CALLdirectnmt_topOnePercentDPrime.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

#get top 5% transcripts
numTranscripts5 <- round(length(avgLD$attribute)*0.05, 0)

topVals5 <- as.data.frame(head(avgLD[order(-avgLD$avgD_prime),], n = numTranscripts5))

nomTopVals5 <- right_join(joined_SNPs, topVals5, by = "attribute")
order.scores <- order(-nomTopVals5$avgD_prime) 
top5 <- nomTopVals5[order.scores,]

write.table(top5, "CALLdirectnmt_topFivePercentDPrime.tsv", row.names = F, sep = '\t',
            col.names = T, quote = F)

  
  
  
  
  