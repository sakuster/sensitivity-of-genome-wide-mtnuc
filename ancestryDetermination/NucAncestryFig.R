library(ggplot2)
library(ggtext)
library(cowplot)
library(viridis)
library(dplyr)

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/07_ancestryDetermination")
CALLNucAnc <- read.table("CALL_ancestryResults_v2.tsv", header = T, sep = "\t", nrow = 281) #282 is avg ancestry number

hist(as.numeric(CALLNucAnc$Percent.Par1.Ancestry), 
     xlim = c(0,1), 
     xlab = "Percent Xbir Ancestry", 
     main = "CALL Nuclear Genome Ancestry")

CHAFNucAnc <- read.table("CHAF_ancestryResults_v2.tsv", header = T, sep = "\t", nrow = 250)
hist(CHAFNucAnc$Percent.Par1.Ancestry, xlim = c(0,1), xlab = "Percent Xbir Ancestry", main = "CHAF Nuclear Genome Ancestry")
hist(CHAFNucAnc$Percent.Par2.Ancestry, xlim = c(0,1), xlab = "Percent Xmal Ancestry", main = "CHAF Nuclear Genome Ancestry")


HCNucAnc <- read.table("HC_ancestryResults_v2.tsv", header = T, sep = "\t", nrow = 255)
hist(HCNucAnc$Percent.Par2.Ancestry, xlim = c(0,1), xlab = "Percent Xcor Ancestry", main = "HUEX-STAC Nuclear Genome Ancestry")

max(HCNucAnc$Percent.Xcor.Ancestry)


#merge mt with nuc CALL
CALLMtAnc <- read.table("CALL_MtAncestryCall.tsv", header = T, sep = "\t")
CALL_merged <- merge(CALLNucAnc, CALLMtAnc, by = "FishID")

CALL_merged <- CALL_merged %>%
  filter(Conclusion!= "Heterozygous",
         Conclusion != "No assignment", 
         Conclusion != "No mtSNPs")

callPlot <- ggplot(CALL_merged, aes(Percent.Par1.Ancestry, fill = Conclusion)) +
  geom_histogram(position = "stack", bins = 12) +#, binwidth = 0.1) +
  scale_fill_manual(values = c("#648FFF", "#DC267F"), 
                    labels=c(expression(italic("X. birchmanni")), 
                             expression(italic("X. malinche")))) + 
  scale_x_continuous(limits = c(0,1), oob = scales::oob_keep) + #kept saying missing vals, is an issue with geom_histogram, this is the workaround
  labs(fill = "Mitochondrial Haplotype", 
       x = expression(paste("Proportion of ", italic("X. birchmanni "), "Nuclear Ancestry")), 
       y = "Number of Individual Fish") +
  theme_light() +
  theme(#axis.text = element_text(size = 15), 
        #axis.title = element_text(size = 17),
        #legend.text = element_text(size = 10),
        #legend.title = element_text(size = 12),
        legend.position = c(0.3, 0.8))

ggsave("CALL_NucByMtAncestry_v2.pdf", callPlot, height = 8, width = 8, units = "in")

#CHAF merge mt with nuc
CHAFMtAnc <- read.table("CHAF_MtAncestryCall.tsv", header = T, sep = "\t")
CHAF_merged <- left_join(CHAFNucAnc, CHAFMtAnc, by = "FishID")

#CHAF_merged <- CHAF_merged[CHAF_merged$Conclusion!="Heterozygous",]
CHAF_merged <- CHAF_merged %>%
  filter(Conclusion!= "Heterozygous",
         Conclusion != "No assignment", 
         Conclusion != "No mtSNPs")

chafPlot <- ggplot(CHAF_merged, aes(Percent.Par2.Ancestry, fill = Conclusion)) +
  geom_histogram(position = "stack", bins = 12) +#, binwidth = 0.1) +
  scale_fill_manual(values = c("#648FFF", "#DC267F"), 
                    labels=c(expression(italic("X. birchmanni")), 
                             expression(italic("X. malinche")))) + 
  #guides(fill = guide_legend(postion = "inside")) +
  scale_x_continuous(limits = c(0,1), oob = scales::oob_keep) + #kept saying missing vals, is an issue with geom_histogram, this is the workaround
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200)) +
  labs(fill = "Mitochondrial Haplotype", 
       x = expression(paste("Proportion of ", italic("X. malinche "), "Nuclear Ancestry")), 
       y = "Number of Individual Fish") +
  theme_light() +
  theme(#axis.text = element_text(size = 15), 
        #axis.title = element_text(size = 17),
        #legend.text = element_text(size = 10),
        #legend.title = element_text(size = 12),
        legend.position = c(0.3, 0.8))

ggsave("CHAF_NucByMtAncestry_v2.pdf", chafPlot, height = 8, width = 8, units = "in")


#joining CALL and CHAF into one figure
callchafPlot <- plot_grid(callPlot, chafPlot, labels = c('A', 'B'), label_size = 12)

ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "Figure2-CALL-CHAF_NucAncestry_v2_Oct2025.pdf", callchafPlot, height = 4, width = 8, units = "in")

