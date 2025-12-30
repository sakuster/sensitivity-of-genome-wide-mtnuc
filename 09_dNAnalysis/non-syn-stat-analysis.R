### non-synonymous analysis -- stats + individual figs

library(dplyr)
library(car)
library(olsrr)
library(broom)

#does number of non-syn subs per gene affect the relationship between mtnuc LD and gene class?
#do interacting n-mts with more non-syn subs have a higher mtnuc LD, indicating higher selection for 'matched' genotypes?

population <- "CALL"

#read in dat
setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")

quantDat <- read.table(paste0(population, "quantDat-all-LD.tsv"), 
                         sep = "\t", header = T, quote = "")

quantDat$class.x <- factor(quantDat$class.x,
                           #levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt"), #make non-n-mt the ref group
                           levels = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) 

sumStats <- quantDat %>%
  group_by(class.x) %>%
  summarise(n = n(),
            meanNS = mean(NS),
            sdNS = sd(NS),
            meanDprime = mean(avgD_prime, na.rm = T),
            sdDprime = sd(avgD_prime, na.rm = T),
            meanRsq = mean(avgR_squared, na.rm = T),
            sdRsq = sd(avgR_squared, na.rm = T),
            meanPcor = mean(avgP_cor, na.rm = T),
            sdPcor = sd(avgP_cor, na.rm = T))

write.table(sumStats, paste0(population, "summary-stats-NS-analysis.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)


#does adding NS actually improve the model fit?
ogMod <- lm(avgD_prime ~ class.x, data = quantDat)
quantModD <- lm(avgD_prime ~ class.x*NS, data = quantDat)
anova(ogMod, 
      quantModD) #for Dprime, including NS improves the model (p = 0.04849)
#not for CHAF though -- p = 0.5374

ogModR <- lm(avgR_squared ~ class.x, data = quantDat)
quantModR <- lm(avgR_squared ~ class.x*NS, data = quantDat)
anova(ogModR,
      quantModR) #r^2 improved by NS (p-val = 0.01019), CHAF p = 0.03519

ogModP <- lm(avgP_cor ~ class.x, data = quantDat)
quantModP <- lm(avgP_cor ~ class.x*NS, data = quantDat)
anova(ogModP,
      quantModP) #pcor improved by NS (p-val = 0.0464), CHAF p = 0.01839


#FOR D'

#create linear model
quantMod <- lm(avgD_prime ~ class.x * NS,
               data = quantDat)
DAnova <- data.frame(Anova(quantMod, type = "III")) #run type III Anova and save to df
summary(quantMod) #test interacting n-mt vs non-n-mt

Dsum <- tidy(quantMod) #save and write summary() -- tidy() does same thing but allows for saving as a table
write.table(Dsum, 
            file = paste0(population,"-non-syn-D-summary.tsv"), 
            sep = "\t", row.names = FALSE)


# R-SQUARED

#create linear model
quantMod <- lm(avgR_squared ~ class.x * NS,
               data = quantDat)
RAnova <- data.frame(Anova(quantMod, type = "III"))
summary(quantMod)

Rsum <- tidy(quantMod)
write.table(Rsum, 
            file = paste0(population, "-non-syn-R-sq-summary.tsv"), 
            sep = "\t", row.names = FALSE)


### PARTIAL CORRELATION

#create linear model
quantMod <- lm(avgP_cor ~ class.x * NS,
               data = quantDat)
pcorAnova <- data.frame(Anova(quantMod, type = "III"))
summary(quantMod)

pcorSum <- tidy(quantMod)
write.table(pcorSum, 
            file = paste0(population, "-non-syn-pcor-summary.tsv"), 
            sep = "\t", row.names = FALSE)


#write ANOVA stats to file
DAnova2 <- DAnova %>%
  mutate(population = population,
         stat = "D'") %>%
  mutate(across(where(is.numeric), round, 4))
RAnova2 <- RAnova %>%
  mutate(population = population,
         stat = "r^2") %>%
  mutate(across(where(is.numeric), round, 4))
pcorAnova2 <- pcorAnova %>%
  mutate(population = population,
         stat = "partial correlation") %>%
  mutate(across(where(is.numeric), round, 4))
allRes <- rbind(DAnova2, RAnova2, pcorAnova2)

write.table(allRes, paste0(population, "-ANOVA-III_results.tsv"), 
            sep = "\t")



#determine influential points for CALL -- this section must be run line by line

interquantDat <- quantDat %>%
  filter(class.x == "interacting n-mt")


#for D' -- no longer need to do bc is not sig post-filtering!
# dMod <- lm(avgD_prime ~ NS,
#            data = interquantDat)
# summary(dMod) #NS p-val = 0.06317, est = 22.04466, r^2 = 0.06317
# 
# ols_plot_dffits(dMod) #points 27, 43, 30 in order of magnitude
# 
# dMod27 <- lm(avgD_prime ~ NS,
#              data = interquantDat[-27,])
# summary(dMod27) #NS p-val = 0.097, est = 16.10471, r^2 = 0.02932
# 
# dMod43 <- lm(avgD_prime ~ NS, 
#              data = interquantDat[-43,])
# summary(dMod43) #NS p-val = 0.0113, est = 24.28277, r^2 = 0.08739
# 
# dMod30 <- lm(avgD_prime ~ NS, 
#              data = interquantDat[-30,])
# summary(dMod30) #NS p-val = 0.0125, est = 25.42222, r^2 = 0.08442
# 
# #see if each affects the OG model
# row27 <- interquantDat[27,]
# row43 <- interquantDat[43,]
# row30 <- interquantDat[30,]
# 
# #ndufs5
# allD28 <- lm(avgD_prime ~ class.x * NS, 
#              data = quantDat %>%
#                anti_join(row27, by = "attribute"))
# Anova(allD28, type = "III") #NS*class p-val = 0.2347
# summary(allD28) #no sig diff for direct n-mts (p = 0.251)
# 
# #mrps30
# allD43 <- lm(avgD_prime ~ class.x * NS, 
#              data = quantDat %>%
#                anti_join(row43, by = "attribute"))
# Anova(allD43, type = "III") #NS*class p-val = 0.0817
# summary(allD43) #no sig diff for direct n-mts (p = 0.0624)
# 
# #mrpl53
# allD30 <- lm(avgD_prime ~ class.x * NS, 
#              data = quantDat %>%
#                anti_join(row30, by = "attribute"))
# Anova(allD30, type = "III") #NS*class p-val = 0.07456
# summary(allD30) #NS*directn-mt p-val = 0.0558



#do the same for r-sq model

rMod <- lm(avgR_squared ~ NS,
           data = interquantDat) 
summary(rMod) #NS p-val = 0.0356, est = 16.70198, r^2 = 0.05511
ols_plot_dffits(rMod) #points 27 in order magnitude

rMod27 <- lm(avgR_squared ~ NS,
             data = interquantDat[-27,])
summary(rMod27) #NS p-val = 0.202, est = 8.959287, r^2 = 0.0108

#see how removing them affects the whole model
row27 <- interquantDat[27,]

allR27 <- lm(avgR_squared ~ class.x * NS, 
             data = quantDat %>%
               anti_join(row27, by = "attribute"))
Anova(allR27, type = "III") #NS*class p-val = 0.07854
summary(allR27) #NS*dirn-mt p-val = 0.2375

##for CHAF r^2:
#do the same for r-sq model

rMod <- lm(avgR_squared ~ NS,
           data = interquantDat) 
summary(rMod) #NS p-val = 0.0356, est = 16.70198, r^2 = 0.05511
ols_plot_dffits(rMod) #points 27, 22, 23, 14 in order magnitude

rMod27 <- lm(avgR_squared ~ NS,
             data = interquantDat[-27,])
summary(rMod27) #NS p-val = 0.154, est = 3.422349, r^2 = 0.01754

rMod22 <- lm(avgR_squared ~ NS,
             data = interquantDat[-22,])
summary(rMod22) #NS p-val = 0.0509, est = 5.179334, r^2 = 0.0464

rMod23 <- lm(avgR_squared ~ NS,
             data = interquantDat[-23,])
summary(rMod23) #NS p-val = 0.0179, est = 8.535311, r^2 = 0.0747

rMod14 <- lm(avgR_squared ~ NS,
             data = interquantDat[-14,])
summary(rMod14) #NS p-val = 0.00651, est = 8.004279, r^2 = 0.1023

#see how removing them affects the whole model
row27 <- interquantDat[27,]
row22 <- interquantDat[22,]
row23 <- interquantDat[23,]
row14 <- interquantDat[14,]

allR27 <- lm(avgR_squared ~ class.x * NS, 
             data = quantDat %>%
               anti_join(row27, by = "attribute"))
Anova(allR27, type = "III") #NS*class p-val = 0.5789
summary(allR27) #NS*dirn-mt p-val = 0.309

allR22 <- lm(avgR_squared ~ class.x * NS, 
             data = quantDat %>%
               anti_join(row22, by = "attribute"))
Anova(allR22, type = "III") #NS*class p-val = 0.2443
summary(allR22) #NS*dirn-mt p-val = 0.0962

allR23 <- lm(avgR_squared ~ class.x * NS, 
             data = quantDat %>%
               anti_join(row23, by = "attribute"))
Anova(allR23, type = "III") #NS*class p-val = 0.0639
summary(allR23) #NS*dirn-mt p-val = 0.0195

allR14 <- lm(avgR_squared ~ class.x * NS, 
             data = quantDat %>%
               anti_join(row14, by = "attribute"))
Anova(allR14, type = "III") #NS*class p-val = 0.03137
summary(allR14) #NS*dirn-mt p-val = 0.00871





#do same for pcor mod

pcorMod <- lm(avgP_cor ~ NS,
              data = interquantDat) 
summary(pcorMod) #NS p-val = 0.0321, est = 21.408839, r^2 = 0.05792
ols_plot_dffits(pcorMod) #only 27

pcorMod27 <- lm(avgP_cor ~ NS,
                data = interquantDat[-27,])
summary(pcorMod27) #NS p-val = 0.145, est = 13.675937, r^2 = 0.01896


#how does removing them affect the whole model?
allp27 <- lm(avgP_cor ~ class.x * NS, 
             data = quantDat %>%
               anti_join(row27, by = "attribute"))
Anova(allp27, type = "III") #NS*class p-val = 0.30248
summary(allp27) #NS*dirn-mt p-val = 0.07115










