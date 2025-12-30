library(ggridges)
library(dplyr)
library(ggplot2)
library(car)
library(ggstatsplot)
library(ggrepel)
library(tidystats)


#set inputs for CALL D'
population <- "CALL"
direct_fn <- paste0(population, "directnmt_avgd-prime.tsv")
indir_fn <- paste0(population, "indirectnmt_avgd-prime.tsv")
non_fn <- paste0(population, "nonnmt_avgd-prime.tsv")
#output_fn <- paste0(population, "-averageD-prime_prettyPlot.pdf")
#output_fn2 <- paste0(population, "-averageD-prime_prettyPlot_withIncompats.pdf")

setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/d-prime/", population))

CALLdir <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% mutate(class = "direct_n-mt")
CALLindir <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% mutate(class = "indirect_n-mt")
CALLnon <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% mutate(class = "non_n-mt")

CALLdat <- rbind(CALLdir, CALLindir, CALLnon)

#mypalette <- c("#", "#A4BDBE", "#457373") #rose, light gray, dark teal
mypalette <- c("#457373", "#A4BDBE", "#D589B1") #rose, light gray, dark teal


CALLdat$class <- factor(CALLdat$class,
                    levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                    labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))

call_plot <- ggplot(CALLdat, aes(y = class, x = avgD_prime, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(axis.text.y=element_blank(), #remove y axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.3, 0.8),
        plot.title = element_text(hjust = 0.5)) +
        #axis.text.x = element_text(size = 9),
        #legend.text = element_text(size = 9), 
        #axis.title = element_text(size = 10)) +
  labs(title = expression(bold("CALL")),
       x = expression(italic("D'")), 
       y = "Gene Class Density") +
  xlim(min(CALLdat$avgD_prime), max(CALLdat$avgD_prime))


call_plot

#ggsave(output_fn, call_plot, width = 7, height = 4)




#set inputs for CHAF D'
population <- "CHAF"
direct_fn <- paste0(population, "directnmt_avgd-prime.tsv")
indir_fn <- paste0(population, "indirectnmt_avgd-prime.tsv")
non_fn <- paste0(population, "nonnmt_avgd-prime.tsv")
#output_fn <- paste0(population, "-averageD-prime_prettyPlot.pdf")
#output_fn2 <- paste0(population, "-averageD-prime_prettyPlot_withIncompats.pdf")

setwd(paste0("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/d-prime/", population))

CHAFdir <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% mutate(class = "direct_n-mt")
CHAFindir <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% mutate(class = "indirect_n-mt")
CHAFnon <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% mutate(class = "non_n-mt")

CHAFdat <- rbind(CHAFdir, CHAFindir, CHAFnon)


CHAFdat$class <- factor(CHAFdat$class,
                        levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                        labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))

chaf_plot <- ggplot(CHAFdat, aes(y = class, x = avgD_prime, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(axis.text.y=element_blank(), #remove y axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  #axis.text.x = element_text(size = 9),
  #legend.text = element_text(size = 9), 
  #axis.title = element_text(size = 10)) +
  labs(title = expression(bold("CHAF")),
       x = expression(italic("D'")), 
       y = "Gene Class Density") +
  xlim(min(CHAFdat$avgD_prime), max(CHAFdat$avgD_prime))


chaf_plot

#joining all plots into one figure
everyone <- plot_grid(call_plot, chaf_plot, 
                      labels = c('A', 'B'), 
                      label_size = 12)
everyone

ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "Figure1-mitonuclearLD_Oct2025.pdf", everyone, height = 4, width = 8, units = "in")



#make supp fig 1
#set inputs for CALL r^2

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/r-squared/CALL")

population <- "CALL"
direct_fn <- paste0(population, "directnmt_avgr-squared.tsv")
indir_fn <- paste0(population, "indirectnmt_avgr-squared.tsv")
non_fn <- paste0(population, "nonnmt_avgr-squared.tsv")
#output_fn <- paste0(population, "-averageR-squared-prettyPlot.pdf")
#output_fn2 <- paste0(population, "-averageR-squared_prettyPlot_withIncompats.pdf")

CALLdirR <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "direct_n-mt")
CALLindirR <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "indirect_n-mt")
CALLnonR <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "non_n-mt")

CALLdatR <- rbind(CALLdirR, CALLindirR, CALLnonR)



CALLdatR$class <- factor(CALLdatR$class,
                         levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                         labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))


CALL_plotR <- ggplot(CALLdatR, aes(y = class, x = avgR_squared, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") + # prevent cutting off of top ridge
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("CALL")),
       x = bquote(italic(r^2)),
       #expression(italic("r %subset% 2" * phantom(0)^"2")), 
       y = "Gene Class Density") +
  xlim(min(0), max(CALLdatR$avgR_squared))


CALL_plotR



#set inputs for CHAF r^2

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/r-squared/CHAF")

population <- "CHAF"
direct_fn <- paste0(population, "directnmt_avgr-squared.tsv")
indir_fn <- paste0(population, "indirectnmt_avgr-squared.tsv")
non_fn <- paste0(population, "nonnmt_avgr-squared.tsv")
#output_fn <- paste0(population, "-averageR-squared-prettyPlot.pdf")
#output_fn2 <- paste0(population, "-averageR-squared_prettyPlot_withIncompats.pdf")

CHAFdirR <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "direct_n-mt")
CHAFindirR <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "indirect_n-mt")
CHAFnonR <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "non_n-mt")

CHAFdatR <- rbind(CHAFdirR, CHAFindirR, CHAFnonR)



CHAFdatR$class <- factor(CHAFdatR$class,
                         levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                         labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))


CHAF_plotR <- ggplot(CHAFdatR, aes(y = class, x = avgR_squared, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") + # prevent cutting off of top ridge
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("CHAF")),
       x = bquote(italic(r^2)),
       #expression(italic("r %subset% 2" * phantom(0)^"2")), 
       y = "Gene Class Density") +
  xlim(min(0), max(CHAFdatR$avgR_squared))


CHAF_plotR


#set inputs for CALL pcor

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/08_partialCorrelation/averagePcor/CALL/")

population <- "CALL"
direct_fn <- paste0(population, "directnmt_avgpcor.tsv")
indir_fn <- paste0(population, "indirectnmt_avgpcor.tsv")
non_fn <- paste0(population, "nonnmt_avgpcor.tsv")
#output_fn <- paste0(population, "-averageR-squared-prettyPlot.pdf")
#output_fn2 <- paste0(population, "-averageR-squared_prettyPlot_withIncompats.pdf")

CALLdirPcor <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "direct_n-mt")
CALLindirPcor <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "indirect_n-mt")
CALLnonPcor <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "non_n-mt")

CALLdatPcor <- rbind(CALLdirPcor, CALLindirPcor, CALLnonPcor)



CALLdatPcor$class <- factor(CALLdatPcor$class,
                         levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                         labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))


CALL_plotPcor <- ggplot(CALLdatPcor, aes(y = class, x = avgP_cor, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") + # prevent cutting off of top ridge
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("CALL")),
       x = "Partial correlation",
       #expression(italic("r %subset% 2" * phantom(0)^"2")), 
       y = "Gene Class Density") +
  xlim(min(CALLdatPcor$avgP_cor, na.rm = T), max(CALLdatPcor$avgP_cor, na.rm = T) + 0.07)


CALL_plotPcor



#set inputs for CHAF pcor

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/08_partialCorrelation/averagePcor/CHAF/")

population <- "CHAF"
direct_fn <- paste0(population, "directnmt_avgpcor.tsv")
indir_fn <- paste0(population, "indirectnmt_avgpcor.tsv")
non_fn <- paste0(population, "nonnmt_avgpcor.tsv")
#output_fn <- paste0(population, "-averageR-squared-prettyPlot.pdf")
#output_fn2 <- paste0(population, "-averageR-squared_prettyPlot_withIncompats.pdf")

CHAFdirPcor <- read.table(direct_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "direct_n-mt")
CHAFindirPcor <- read.table(indir_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "indirect_n-mt")
CHAFnonPcor <- read.table(non_fn, header = TRUE, sep = "\t", quote = "") %>% 
  mutate(class = "non_n-mt")

CHAFdatPcor <- rbind(CHAFdirPcor, CHAFindirPcor, CHAFnonPcor)



CHAFdatPcor$class <- factor(CHAFdatPcor$class,
                         levels = c("direct_n-mt", "indirect_n-mt", "non_n-mt"),
                         labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))


CHAF_plotPcor <- ggplot(CHAFdatPcor, aes(y = class, x = avgP_cor, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") + # prevent cutting off of top ridge
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("CHAF")),
       x = "Partial correlation",
       #expression(italic("r %subset% 2" * phantom(0)^"2")), 
       y = "Gene Class Density") +
  xlim(min(CHAFdatPcor$avgP_cor, na.rm = T), max(CHAFdatPcor$avgP_cor, na.rm = T) + 0.07)


CHAF_plotPcor



#joining all plots into one figure
suppfig1 <- plot_grid(CALL_plotR, CHAF_plotR, CALL_plotPcor, CHAF_plotPcor,
                      labels = c('A', 'B', 'C', 'D'), 
                      label_size = 12)

suppfig1

ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "SuppFigure1-mitonuclearLDmetrics_Oct2025.pdf", suppfig1, height = 8, width = 8, units = "in")







