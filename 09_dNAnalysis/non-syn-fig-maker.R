#non-synonymous analysis
#figure maker
library(dplyr)
library(ggrepel)
library(ggplot2)
library(cowplot)

setwd("/Users/kusters/Documents/XiphoStartingOver_Aug2024/09_dNdSAnalysis/")


#read in CALL data:
CALLquantDat <- read.table("CALLquantDat-all-LD.tsv", 
                       sep = "\t", header = T, quote = "")

CALLquantDat$class.x <- factor(CALLquantDat$class.x,
                           #levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt"), #make non-n-mt the ref group
                           levels = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) 

#read in CHAF data:
CHAFquantDat <- read.table("CHAFquantDat-all-LD.tsv", 
                           sep = "\t", header = T, quote = "")

CHAFquantDat$class.x <- factor(CHAFquantDat$class.x,
                               levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt"), #make non-n-mt the ref group
                               labels = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) 

#read in HUEXSTAC data:
HCquantDat <- read.table("HUEXSTACquantDat-all-LD.tsv", 
                         sep = "\t", header = T, quote = "")

HCquantDat$class.x <- factor(HCquantDat$class.x,
                               levels = c("non-n-mt", "indirect_n-mt", "direct_n-mt"), #make non-n-mt the ref group
                               labels = c("non-n-mt", "non-interacting n-mt", "interacting n-mt")) 


#make CALL plots:
#make figs
mypalette <- c("#D589B1", "#A4BDBE", "#457373")

#D'
CALLD <- ggplot(CALLquantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = NS, y = avgD_prime)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", 
              se = FALSE, alpha = 0.2, show.legend = FALSE,
              aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = element_blank(), #"Non-synonymous substitution count / CDS length (bp)",
       y = expression(italic("D'"))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("CALLD-non-syn-fig.pdf", CALLD, width = 7, height = 4)


#r-sq
CALLR <- ggplot(CALLquantDat %>%
                  mutate(class.x = fct_rev(class.x)), 
                         aes(x = NS, y = avgR_squared)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, show.legend = FALSE,
              aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Non-synonymous substitution count / CDS length (bp)",
       y = bquote(italic(r^2))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("CALLRsq-non-syn-fig.pdf", CALLR, width = 7, height = 4)


#pcor
CALLPcor <- ggplot(CALLquantDat %>%
         mutate(class.x = fct_rev(class.x)), 
       aes(x = NS, y = avgP_cor)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, show.legend = FALSE,
              aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = element_blank(), #"Non-synonymous substitution count / CDS length (bp)",
       y = "Partial correlation score") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("CALLPcor-non-syn-fig.pdf", CALLPcor, width = 7, height = 4)


#make CHAF individual plots:
#D'
CHAFD <- ggplot(CHAFquantDat %>%
                  mutate(class.x = fct_rev(class.x)), 
                aes(x = NS, y = avgD_prime)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", show.legend = FALSE,
              se = FALSE, alpha = 0.2, aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = element_blank(), #"Non-synonymous substitution count / CDS length (bp)",
       y = expression(italic("D'"))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("CHAFD-non-syn-fig.pdf", CHAFD, width = 7, height = 4)


#r-sq
CHAFR <- ggplot(CHAFquantDat %>%
                  mutate(class.x = fct_rev(class.x)), 
                aes(x = NS, y = avgR_squared)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, show.legend = FALSE,
              aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Non-synonymous substitution count / CDS length (bp)",
       y = bquote(italic(r^2))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("CHAFR-sq-non-syn-fig.pdf", CHAFR, width = 7, height = 4)


#pcor
CHAFPcor <- ggplot(CHAFquantDat %>%
                     mutate(class.x = fct_rev(class.x)), 
                   aes(x = NS, y = avgP_cor)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, show.legend = FALSE,
              aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = element_blank(), #"Non-synonymous substitution count / CDS length (bp)",
       y = "Partial correlation score") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("CHAFPcor-non-syn-fig.pdf", CHAFPcor, width = 7, height = 4)


#HUEXSTAC allelefreq
HCaf <- ggplot(HCquantDat %>%
                     mutate(class.x = fct_rev(class.x)), 
                   aes(x = NS, y = avgAFreq)) +
  geom_point(alpha = 0.15, show.legend = FALSE,
             aes(color = class.x)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, show.legend = FALSE,
              aes(color = class.x)) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~class.x) +
  theme_bw() +
  labs(x = "Non-synonymous substitution count / CDS length (bp)",
       y = "Allele Frequency") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("HC-allelefreq-non-syn-fig.pdf", HCaf, width = 7, height = 4)



#make fig 4 -- D' only for both CALL, CHAF
Dplot <- plot_grid(CALLD, CHAFD, HCaf,
          labels = c('A', 'B', 'C'), 
          nrow = 3, 
          ncol = 1,
          label_size = 12)

ggsave("Fig4-non-syn-analysis.pdf", Dplot, 
       width = 8, height = 8, units = "in")

Pcorplot <- plot_grid(CALLPcor, CHAFPcor, HCaf,
                      labels = c('CALL', 'CHAF', 'HUEXSTAC'), 
                      nrow = 3, 
                      ncol = 1,
                      label_size = 12,
                      hjust = 0, label_x = 0.01)

ggsave("Fig4-non-syn-analysis.pdf", Pcorplot, 
       width = 8, height = 8, units = "in")


#make supp figs
#for CALL
CALLsupp <- plot_grid(CALLD, CALLR,
                      labels = c('A', 'B'), 
                      nrow = 2, 
                      ncol = 1,
                      label_size = 12)

ggsave("CALL_r-sq-d-prime-non-syn.pdf", CALLsupp, 
       width = 8, height = 8, units = "in")

#for CHAF
CHAFsupp <- plot_grid(CHAFD, CHAFR,
                      labels = c('A', 'B'), 
                      nrow = 2, 
                      ncol = 1,
                      label_size = 12)

ggsave("CHAF_r-sq-pcor-non-syn.pdf", CHAFsupp, 
       width = 8, height = 8, units = "in")

















