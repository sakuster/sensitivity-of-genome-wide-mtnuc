library(ggplot2)
library(ggtext)
library(cowplot)
library(viridis)
library(dplyr)

#to visualize nuclear & mitochondrial ancestry per fish
#CALL par1 = Xbir (major parent), par2 = Xmal
#CHAF par1 = Xbir, par2 = Xmal (major parent)
#HUEX-STAC par1 = Xbir, par2 = Xcor (major parent)

setwd("")
CALLNucAnc <- read.table("CALL_ancestryResults.tsv", header = T, sep = "\t", nrow = 281) #282 is avg ancestry number

hist(as.numeric(CALLNucAnc$Percent.Par1.Ancestry), 
     xlim = c(0,1), 
     xlab = "Percent Xbir Ancestry", 
     main = "CALL Nuclear Genome Ancestry")

CHAFNucAnc <- read.table("CHAF_ancestryResults.tsv", header = T, sep = "\t", nrow = 250)
hist(CHAFNucAnc$Percent.Par1.Ancestry, xlim = c(0,1), xlab = "Percent Xbir Ancestry", main = "CHAF Nuclear Genome Ancestry")
hist(CHAFNucAnc$Percent.Par2.Ancestry, xlim = c(0,1), xlab = "Percent Xmal Ancestry", main = "CHAF Nuclear Genome Ancestry")


HCNucAnc <- read.table("HC_ancestryResults.tsv", header = T, sep = "\t", nrow = 255)
hist(HCNucAnc$Percent.Par2.Ancestry, xlim = c(0,1), xlab = "Percent Xcor Ancestry", main = "HUEX-STAC Nuclear Genome Ancestry")



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

ggsave("CALL_NucByMtAncestry.pdf", callPlot, height = 8, width = 8, units = "in")

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
  labs(fill = "Mitochondrial Haplotype", 
       x = expression(paste("Proportion of ", italic("X. malinche "), "Nuclear Ancestry")), 
       y = "Number of Individual Fish") +
  theme_light() +
  theme(#axis.text = element_text(size = 15), 
        #axis.title = element_text(size = 17),
        #legend.text = element_text(size = 10),
        #legend.title = element_text(size = 12),
        legend.position = c(0.3, 0.8))

ggsave("CHAF_NucByMtAncestry.pdf", chafPlot, height = 8, width = 8, units = "in")


#joining CALL and CHAF into one figure
callchafPlot <- plot_grid(callPlot, chafPlot, labels = c('A', 'B'), label_size = 12)

ggsave(path = "/Users/kusters/Library/CloudStorage/OneDrive-Colostate/Research/Xiphophorus/Writing",
       "Figure2-CALL-CHAF_NucAncestry.pdf", callchafPlot, height = 4, width = 8, units = "in")


#HC merge mt with nuc
HCMtAnc <- read.table("HC_MtAncestryCall.tsv", header = T, sep = "\t")
HC_merged <- merge(HCNucAnc, HCMtAnc, by = "FishID")

HC_merged <- HC_merged %>%
  filter(Conclusion!= "Heterozygous",
         Conclusion != "No assignment", 
         Conclusion != "No mtSNPs",
         FishID != "HUEX-XI-19-80.R1.fastq") #is F1 hybrid somehow?

hcPlot <- ggplot(HC_merged, aes(Percent.Par2.Ancestry, fill = Conclusion)) +
  geom_histogram(position = "stack", bins = 12) +
  scale_fill_manual(values = c("#FE6100"), 
                    labels=c(expression(italic("X. cortezi")))) + 
  scale_x_continuous(limits = c(0,1), oob = scales::oob_keep) + #kept saying missing vals, is an issue with geom_histogram, this is the workaround
  labs(fill = "Mitochondrial Haplotype", 
       x = expression(paste("Proportion of ", italic("X. cortezi "), "Nuclear Ancestry")), 
       y = "Number of Individual Fish") +
  theme_light() +
  theme(#axis.text = element_text(size = 15),
        #axis.title = element_text(size = 17),
        #legend.text = element_text(size = 12),
        #legend.title = element_text(size = 15),
        legend.position = c(0.3, 0.8))

ggsave("HC_NucByMtAncestry.pdf", hcPlot, height = 8, width = 8, units = "in")


sd(CALLNucAnc$Percent.Par1.Ancestry)
sd(CHAFNucAnc$Percent.Par2.Ancestry)
sd(HCNucAnc$Percent.Par2.Ancestry)

mean(CALLNucAnc$Percent.Par1.Ancestry)
mean(CHAFNucAnc$Percent.Par2.Ancestry)
mean(HCNucAnc$Percent.Par2.Ancestry)

x - 2*a
x+ 2*a

y - 2*b
y + 2*b

z - 2*c
z + 2*c

