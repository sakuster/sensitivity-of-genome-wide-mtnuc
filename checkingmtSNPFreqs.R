##checking that mtSNPs all show the same or similar post probab

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/00_StartingDataFiles")

call <- read.table("CALL_mtSNPs.tsv", sep = "\t", header = T, quote = "")
chaf <- read.table("CHAF_mtSNPs.tsv", sep = "\t", header = T, quote = "")
hc <- read.table("HC_mtSNPs.tsv", sep = "\t", header = T, quote = "")

for (i in 1:nrow(call)) {
  snpNum <- length(call[1,]) #how many SNPs are there for each fish? 
  sumProbs <- sum(call[1,]) #what is the sum of post. probabilities for the fish?
  if (snpNum != sumProbs & sumProbs > 1) { #does this match either 0 or the # SNPs?
    print(paste0("row ", i, " in CALL equals ", sumProbs, " not 0 or number of SNPs"))
  } 
}
print("call has finished. if you did not receive another message, this population's mtDNA posterior probabilities match the expectation that mtDNA is non-recombining and linked")

for (i in 1:nrow(chaf)) {
  snpNum <- length(chaf[1,]) #how many SNPs are there for each fish? 
  sumProbs <- sum(chaf[1,]) #what is the sum of post. probabilities for the fish?
  if (snpNum != sumProbs & sumProbs > 1) { #does this match either 0 or the # SNPs?
    print(paste0("row ", i, " in CHAF equals ", sumProbs, " not 0 or number of SNPs"))
  } 
}
print("CHAF has finished. if you did not receive another message, this population's mtDNA posterior probabilities match the expectation that mtDNA is non-recombining and linked")

for (i in 1:nrow(hc)) {
  snpNum <- length(hc[1,]) #how many SNPs are there for each fish? 
  sumProbs <- sum(hc[1,]) #what is the sum of post. probabilities for the fish?
  if (snpNum != sumProbs & sumProbs > 1) { #does this match either 0 or the # SNPs?
    print(paste0("row ", i, " in HUEXSTAC equals ", sumProbs, " not 0 or number of SNPs"))
  } 
}
print("HUEXSTAC has finished. if you did not receive another message, this population's mtDNA posterior probabilities match the expectation that mtDNA is non-recombining and linked")

