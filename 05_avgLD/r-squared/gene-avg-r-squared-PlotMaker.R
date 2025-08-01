###make pretty figures & stats

library(ggridges)
library(dplyr)
library(ggplot2)
library(car)

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/r-squared/CHAF")

#set inputs
population <- "CHAF"
direct_fn <- paste0(population, "directnmt_avgr-squared.tsv")
indir_fn <- paste0(population, "indirectnmt_avgr-squared.tsv")
non_fn <- paste0(population, "nonnmt_avgr-squared.tsv")
output_fn <- paste0(population, "-averageR-squared-prettyPlot.pdf")
output_fn2 <- paste0(population, "-averageR-squared_prettyPlot_withIncompats.pdf")
mosPlot2_fn <- paste0(population, "_2groupsmosaicPlot.pdf")
mosPlot3_fn <- paste0(population, "_3groupsmosaicPlot.pdf")
fischerbarstat_fn <- paste0(population, "fischerBarPlot.pdf")

dir <- read.table(direct_fn, header = TRUE, sep = "\t") %>% mutate(class = "direct_n-mt")
indir <- read.table(indir_fn, header = TRUE, sep = "\t") %>% mutate(class = "indirect_n-mt")
non <- read.table(non_fn, header = TRUE, sep = "\t") %>% mutate(class = "non_n-mt")

dat <- rbind(dir, indir, non)

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

incompat_genes$class <- factor(incompat_genes$class,
                               levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                               labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))

dat$class <- factor(dat$class,
                    levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                    labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))


#make plot
all_plot <- ggplot(dat, aes(y = class, x = avgR_squared, fill = class)) + 
  geom_density_ridges() +
  theme_classic() +
  coord_cartesian(clip = "off") + # prevent cutting off of top ridge
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank()) +
  scale_fill_manual(values = mypalette) +
  labs(title = paste0(population, " Gene Averaged LD"), x = "r^2", y = "Gene Class") +
  xlim(min(0), max(dat$avgR_squared))

all_plot

ggsave(output_fn, all_plot, width = 7, height = 4)

mypalette <- c("#457373", "#A4BDBE", "#D589B1")
incompat_plot <- ggplot(dat, aes(y = class, x = avgR_squared, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette,
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  geom_vline(xintercept = mean(dat$avgR_squared), linetype = "dashed") +
  geom_point(data = incompat_genes, 
             aes(y = class, x = avgR_squared),
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
        legend.position = c(0.8,0.8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("CHAF")), 
       x = bquote(italic(r^2)), 
       y = "Gene Class Density") +
  xlim(min(dat$avgR_squared), max(dat$avgR_squared))

incompat_plot

ggsave(output_fn2, incompat_plot, width = 7, height = 4)
ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "SuppFig_CHAF-r-squred-incompatibilityPlot.pdf", plot = incompat_plot, width = 8, height = 5)


#one way anova test
dat2 <- rbind(non, indir, dir)
dat2$class <- as.factor(dat2$class)
dat_df <- as.data.frame(dat2)

sumStats <- group_by(dat_df, class) %>% summarise(count = n(), mean = mean(avgR_squared), sd(avgR_squared))

OneWayFit <- lm(avgR_squared ~ class, data = dat_df)
aov <- anova(OneWayFit)

write.table(sumStats, paste0(population, "_summaryStats.txt"),
            row.names = F, sep = '\t', col.names = T, quote = F)

write.table(aov, paste0(population, "_ANOVA.txt"),
            row.names = T, sep = '\t', col.names = T, quote = F)

#get top 1% transcripts
numTranscripts <- length(dat$attribute)*0.01

topVals <- as.data.frame(head(dat[order(-dat$avgR_squared),], n = numTranscripts))
snp <- read.table(paste0(population, "_SNPinfo.tsv"), header = T, sep = '\t')
nomTopVals <- right_join(snp, topVals, by = "attribute") #can maintain order with inner_join() but doesn't work
order.scores <- order(-nomTopVals$avgR_squared) 
top1 <- nomTopVals[order.scores,]


write.table(top1, paste0(population, "_ALLGENES_topOnePercentRSquared.tsv"), col.names = T,
            row.names = F)

# Fischer's exact test -- does the top 1% represent n-mts disporportionately?

# split by n-mt and non-n-mt
statDat <- dat %>%
  mutate(topOne = case_when(attribute %in% top1$attribute ~ "yes",
                            !attribute %in% top1$attribute ~ "no"),
         n_mt =  case_when(class == "direct_n-mt" | class == "indirect_n-mt" ~ "yes",
                           class == "non_n-mt" ~ "no")) %>%
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
                  as.numeric(newRow[which(newRow$category == "no_non_n-mt"),3])),
  "topOne_yes" = c(as.numeric(newRow[which(newRow$category == "yes_direct_n-mt"),3]), 
                   as.numeric(newRow[which(newRow$category == "yes_indirect_n-mt"),3]),
                   as.numeric(newRow[which(newRow$category == "yes_non_n-mt"),3])),
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
