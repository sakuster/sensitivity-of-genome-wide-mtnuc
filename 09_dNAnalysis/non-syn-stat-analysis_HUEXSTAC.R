### non-synonymous analysis -- stats + individual figs

library(dplyr)
library(car)
library(olsrr)
library(broom)

#does number of non-syn subs per gene affect the relationship between mtnuc LD and gene class?
#do interacting n-mts with more non-syn subs have a higher mtnuc LD, indicating higher selection for 'matched' genotypes?

population <- "HUEXSTAC"

#read in dat
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

quantDat <- read.table(paste0(population, "quantDat-all-LD.tsv"), 
                       sep = "\t", header = T, quote = "")

quantDat$class.x <- factor(quantDat$class.x,
                           levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt"), #make non-n-mt the ref group
                           labels = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) 

sumStats <- quantDat %>% #is already filtered by fish depth
  group_by(class.x) %>%
  summarise(n = n(),
            meanNS = mean(NS),
            sdNS = sd(NS),
            meanAfreq = mean(avgAFreq, na.rm = T))

write.table(sumStats, paste0(population, "summary-stats-NS-analysis.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)


#does adding NS actually improve the model fit?
ogMod <- lm(avgAFreq ~ class.x, data = quantDat)
quantModD <- lm(avgAFreq ~ class.x*NS, data = quantDat)
anova(ogMod, 
      quantModD) #for HC, including NS marginally improves model (p = 0.02957)



#allle frequency

#create linear model
quantMod <- lm(avgAFreq ~ class.x * NS,
               data = quantDat)
anovaT <- data.frame(Anova(quantMod, type = "III")) #run type III Anova and save to df
summary(quantMod) #test interacting n-mt vs non-n-mt

Dsum <- tidy(quantMod) #save and write summary() -- tidy() does same thing but allows for saving as a table
write.table(Dsum, 
            file = paste0(population,"-non-syn-allelefreq-summary.tsv"), 
            sep = "\t", row.names = FALSE)



#write ANOVA stats to file
DAnova2 <- anovaT %>%
  mutate(population = population,
         stat = "allelefreq") %>%
  mutate(across(where(is.numeric), round, 4))


write.table(DAnova2, paste0(population, "-ANOVA-III_results.tsv"), 
            sep = "\t")










