###make pretty figures & stats

library(ggridges)
library(dplyr)
library(ggplot2)
library(car)
library(ggstatsplot)
library(ggrepel)
library(tidystats)


#set inputs
population <- "HUEXSTAC"
species <- "Xcor"
speciesLong <- "X. cortezi"

all_fn <- paste0(population, "_avg", species, "a-freq_allclasses.tsv")
output_fn <- paste0(population, "-averagea-freq_prettyPlot.pdf")
output_fn2 <- paste0(population, "-averagea-freq_prettyPlot_withIncompats.pdf")
mosPlot2_fn <- paste0(population, "_2groupsmosaicPlot.pdf")
mosPlot3_fn <- paste0(population, "_3groupsmosaicPlot.pdf")
fischerbarstat_fn <- paste0(population, "fischerBarPlot.pdf")

setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/HC_allelefreq/", species, "AlleleFreq"))

dat <- read.table(all_fn, header = TRUE, sep = "\t", quote = "")

mypalette <- c("#D589B1", "#A4BDBE", "#457373") #rose, light gray, dark teal

geneList <- data.frame(attribute = c("g8054.t1", #Ndufs5
                                     "g3279.t1", #Ndufa13
                                     "g3322.t1", #MTERF4 ; UQCR11 not found
                                     "g6088.t1", #ATP5MG -- isn't in CHAF for some reason??
                                     "g6066.t1", #C1QBP
                                     "g15308.t1", #MMUT
                                     "g15060.t1", #SMIM8
                                     "g15051.t1", #LYRM2
                                     "g15056.t1", #RMDN3
                                     "g11575.t1"), #UQCRC2
                       GN = c("NDUFS5", "NDUFA13","MTERF4", "ATP5MG", "C1QBP",
                              "MMUT", "SMIM8", "LYRM2", "RMDN3", "UQCRC2")) 

incompat_genes <- inner_join(geneList, dat, by = "attribute")

dat$class <- factor(dat$class,
                        levels = c("direct_n-mt", "indirect_n-mt", "non-n-mt"),
                        labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))
incompat_genes$class <- factor(incompat_genes$class,
                        levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                        labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))


all_plot <- ggplot(dat, aes(y = class, x = avgAFreq, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette) +
  theme_classic() +
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  labs(title = paste0(population, " Gene Averaged ", species, " Allele Frequencies"), 
       x = "Allele Frequency", y = "Gene Class") +
  xlim(min(dat$avgAFreq), max(dat$avgAFreq))

all_plot

ggsave(output_fn, all_plot, width = 7, height = 4)

mypalette <- c("#457373", "#A4BDBE", "#D589B1")
#plot with incompatibility genes marked
incompat_plot <- ggplot(dat, aes(y = class, x = avgAFreq, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette,
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  geom_vline(xintercept = mean(dat$avgAFreq), linetype = "dashed") +
  geom_point(data = incompat_genes, 
             aes(y = class, x = avgAFreq),
             show.legend=FALSE) +
  geom_label_repel(data = incompat_genes, 
                   aes(label = GN),
                   alpha = 0.7,
                   show.legend=FALSE) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(axis.text.y=element_blank(), #remove y axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2,0.8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("HUEX-STAC")),
       x = "Average Allele Frequency", 
       y = "Gene Class Density") +
  xlim(min(dat$avgAFreq), max(dat$avgAFreq))

incompat_plot

ggsave(output_fn2, incompat_plot, width = 7, height = 4)
ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "SuppFig5_HUEX-STAC-incompatibilityPlot.pdf", width = 8, height = 5)

#one way anova test
dat2 <- dat
dat2$class <- as.factor(dat2$class)

sumStats <- dat2 %>%
  group_by(class) %>% 
  summarise(count = n(), mean = mean(avgAFreq), sd(avgAFreq))

OneWayFit <- lm(avgAFreq ~ class, data = dat2)
aov <- anova(OneWayFit)

write.table(sumStats, paste(population, species, "summaryStats.txt", sep = "_"),
            row.names = F, sep = '\t', col.names = T, quote = F)

write.table(aov, paste(population, species, "_ANOVA.txt", sep = "_"),
            row.names = T, sep = '\t', col.names = T, quote = F)



#get top 1% transcripts -- done in HUEXSTAC_allclasses_avgAlleleFreq.R
top1 <- read.table(paste0(speciesLong, "_", population, "allclasses_topOnePercentAFreq.tsv"),
                     header = T,
                     sep = "\t", 
                     quote = "") %>%
  group_by(attribute)

# Fischer's exact test -- does the top 1% represent n-mts disporportionately?

# split by n-mt and non-n-mt
statDat <- dat %>%
  mutate(topOne = case_when(attribute %in% top1$attribute ~ "yes",
                            !attribute %in% top1$attribute ~ "no"),
         n_mt =  case_when(class == "direct_n-mt" | class == "indirect_n-mt" ~ "yes",
                           class == "non-n-mt" ~ "no")) %>%
  group_by(topOne, n_mt) %>%
  summarise(n = n())

numsStatTest <- data.frame(
  "topOne_no" = c(as.numeric(statDat[2,3]),as.numeric(statDat[1,3])),
  "topOne_yes" = c(as.numeric(statDat[4,3]), as.numeric(statDat[3,3])),
  row.names = c("n-mt", "non-n-mt"),
  stringsAsFactors = F
)

pdf(file = mosPlot2_fn)
mosaicplot(numsStatTest,
                      main = "Mosaic plot",
                      color = TRUE)
dev.off()


expectedChi <- chisq.test(numsStatTest)$expected #indicates whether a chi-sq or fisher should be used

if (expectedChi[1,1] > 5 & expectedChi[1,2] > 5 & expectedChi[2,1] > 5 & expectedChi[2,2] > 5) {
  g2Test <- chisq.test(numsStatTest)
} else {
  g2Test <- fisher.test(numsStatTest)
}


# split by interacting n-mt, non-interacting n-mt, non-n-mt
statDat2 <- dat %>%
  mutate(topOne = case_when(attribute %in% top1$attribute ~ "yes",
                            !attribute %in% top1$attribute ~ "no")) %>%
  group_by(topOne, class) %>%
  summarise(n = n())

newRow <- statDat2 %>% 
  mutate(category = paste(topOne, class, sep = "_"))

numsStatTest2 <- data.frame(
  "topOne_no" = c(as.numeric(newRow[which(newRow$category == "no_direct_n-mt"),3]), 
                  as.numeric(newRow[which(newRow$category == "no_indirect_n-mt"),3]),
                  as.numeric(newRow[which(newRow$category == "no_non-n-mt"),3])),
  "topOne_yes" = c(as.numeric(newRow[which(newRow$category == "yes_direct_n-mt"),3]), 
                   as.numeric(newRow[which(newRow$category == "yes_indirect_n-mt"),3]),
                   as.numeric(newRow[which(newRow$category == "yes_non-n-mt"),3])),
  row.names = c("interacting-n-mt", "non-interacting-n-mt", "non-n-mt"),
  stringsAsFactors = F
)

for (row in 1:length(numsStatTest2$topOne_no)) {
  for (col in 1:length(numsStatTest2[1,])) {
    if (is.na(numsStatTest2[row,col]) == T) {
      numsStatTest2[row,col] <- 0
    }
  }
}

pdf(file = mosPlot3_fn)
mosPlot <- mosaicplot(numsStatTest2,
                      main = "Mosaic plot",
                      color = TRUE)
dev.off()

expectedChi2 <- chisq.test(numsStatTest2)$expected #indicates whether a chi-sq or fisher should be used

if (expectedChi2[1,1] > 5 & expectedChi2[1,2] > 5 & expectedChi2[2,1] > 5 & expectedChi2[2,2] > 5) {
  test <- chisq.test(numsStatTest2)
} else {
  test <- fisher.test(numsStatTest2, alternative = "greater")
}

x <- c() #starting here, code modified from https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/
for (row in rownames(numsStatTest2)) {
  for (col in colnames(numsStatTest2)) {
    x <- rbind(x, matrix(rep(c(row, col), numsStatTest2[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
colnames(df) <- c("GeneClass", "TopOnePercentAvgLD")
df

#test <- chisq.test(table(df))

# combine plot and statistical test with ggbarstats
#library(ggstatsplot)
pdf(file = fischerbarstat_fn)
ggbarstats(
  df, GeneClass, TopOnePercentAvgLD,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
)
dev.off()

statList <- list()
statList <- statList |>
  add_stats(g2Test) |>
  add_stats(test)

write_stats(statList, paste0(population, "statistics.txt"))


