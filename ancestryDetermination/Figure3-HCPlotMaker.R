#makes Figure 3 for sensitivity of genome-wide tests for mitonuclear coevolution MS
#which is huex-stac allele freq and nuclear ancestry

library(cowplot)

#read in HC ancestry data
setwd("~/Documents/XiphoStartingOver_Aug2024/07_ancestryDetermination")

HCNucAnc <- read.table("HC_ancestryResults_v2.tsv", header = T, sep = "\t", nrow = 255)

#HC merge mt with nuc
HCMtAnc <- read.table("HC_MtAncestryCall.tsv", header = T, sep = "\t")
HC_merged <- merge(HCNucAnc, HCMtAnc, by = "FishID")

HC_merged <- HC_merged %>%
  filter(Conclusion!= "Heterozygous",
         Conclusion != "No assignment", 
         Conclusion != "No mtSNPs",
         FishID != "HUEX-XI-19-80.R1.fastq") #something is up w/ ancestry calls


hcPlot <- ggplot(HC_merged, aes(Percent.Par2.Ancestry, fill = Conclusion)) +
  geom_histogram(position = "stack", bins = 12) +
  scale_fill_manual(values = c("#FE6100")) + 
  scale_x_continuous(limits = c(0,1), oob = scales::oob_keep) + #kept saying missing vals, is an issue with geom_histogram, this is the workaround
  labs(fill = "Mitochondrial Haplotype", 
       x = expression(paste("Proportion of ", italic("X. cortezi "), "Nuclear Ancestry")), 
       y = "Number of Individual Fish") +
  theme_light() +
  theme(#axis.text = element_text(size = 15),
    #axis.title = element_text(size = 17),
    #legend.text = element_text(size = 12),
    #legend.title = element_text(size = 15),
    legend.position = "none")


#read in HC allele freq (xcor) data

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

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/05_avgLD/HC_allelefreq/10-21-25_HCallelefreq/v4_HCalleleFreq/")

HCdat <- read.table(all_fn, header = TRUE, sep = "\t", quote = "")

HCdat$class <- factor(HCdat$class,
                         levels = c("direct_n-mt", "indirect_n-mt", "non-n-mt"),
                         labels = c("interacting n-mt","non-interacting n-mt", "non-n-mt"))



mypalette <- c("#457373", "#A4BDBE", "#D589B1" ) #rose, light gray, dark teal

HC_plot <- ggplot(HCdat, aes(y = class, x = avgAFreq, fill = class)) + 
  geom_density_ridges() +
  scale_fill_manual(values = mypalette, 
                    limits = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) +
  theme_classic() +
  coord_cartesian(clip = "off") + # prevent cutting off of top ridge
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.4, 0.8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = expression(bold("HUEX-STAC")),
       x = "Average Allele Frequency",
       #expression(italic("r %subset% 2" * phantom(0)^"2")), 
       y = "Gene Class Density") +
  xlim(min(0), max(HCdat$avgAFreq))


HC_plot

#ggsave(output_fn, all_plot, width = 7, height = 4)



#joining all plots into one figure
everyone <- plot_grid(HC_plot, hcPlot,
                      labels = c('A', 'B'), 
                      label_size = 12)

ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "Figure3-HUEX-STAC_Oct2025.pdf", everyone, height = 4, width = 8, units = "in")





