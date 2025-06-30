#permutation test to confirm HC significant ANOVA
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(lmPerm)
library(coin)
library(gtools)

#calculate difference in means for HC gene classes
dir <- read.table("HUEXSTAC_avgXcora-freq_directnmt.tsv", header = T, quote = "", sep = "\t") %>%
  mutate(class = "direct_n-mt")
indir <- read.table("HUEXSTAC_avgXcora-freq_indirectnmt.tsv", header = T, quote = "", sep = "\t") %>%
  mutate(class = "non_direct_n-mt")
non <- read.table("HUEXSTAC_avgXcora-freq_nonnmt.tsv", header = T, quote = "", sep = "\t")  %>%
  mutate(class = "non_direct_n-mt")

nonDir <- rbind(indir, non)

dirMean <- mean(dir$avgAFreq)
nonDirMean <- mean(nonDir$avgAFreq)
obsDiffMean <- dirMean - nonDirMean

allDat <- rbind(dir, indir, non)

#code modified from https://mac-theobio.github.io/QMEE/lectures/permutation_examples.notes.html
set.seed(101) ## for reproducibility
nsim <- 9999
res <- numeric(nsim) ## set aside space for results

for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(allDat))
  bdat <- transform(allDat,avgAFreq=avgAFreq[perm])
  ## compute & store difference in means; store the value
  res[i] <- mean(bdat$avgAFreq[bdat$class=="direct_n-mt"]) -
    mean(bdat$avgAFreq[bdat$class=="non_direct_n-mt"])
}

## append the observed value to the list of results
res <- c(res,obsDiffMean)

pdf(file = "permutationTestResults.pdf")
hist(res,col="gray",las=1,
main="Results of HUEXSTAC Allele Frequency Permutation Test",
     xlab = "(Mean Interacting N-mt Allele Frequency) - (Mean All Other Allele Frequencies)")
abline(v=obsDiffMean,col="red")
dev.off()

num <- length(which(res >= obsDiffMean))
pVal <- num / 10000

fileConn<-file("PermutationTestResults.txt")
writeLines(c("After running ", as.character(nsim), " simulations with your permutated data, your p-value is: ", as.character(pVal)), fileConn)
close(fileConn)

dat <- as.data.frame(res)
prettyPlot <- ggplot(dat, aes(x = res)) +
  geom_histogram(fill = "gray", color = "gray4") +
  geom_vline(xintercept = obsDiffMean, color = "red", linewidth = 1) +
  theme_classic() +
  labs(x = "Difference Between Permutated Interacting N-mt Allele Frequency\n and Allele Frequency of All Other Gene Classes",
       y = "Density of Permutated Differences")

prettyPlot

ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "SuppFig-HUEXSTAC-permutationTest.pdf", prettyPlot, width = 8, height = 8)
