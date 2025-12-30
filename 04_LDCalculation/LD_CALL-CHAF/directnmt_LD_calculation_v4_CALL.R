###Calculate LD on Xiphophorus data
###Shady Kuster
###21 July 2022
library(tidyr)
library(ggplot2)

#this has been modified for ovis

#set output filenames (modify each time)
out_fn_ALL <- "mt1_results_directnmt_CALL_ALL.tsv"
plot1_fn <- "CALL_mt1_LD_plot_DPRIME_directnmts.pdf"
plot2_fn <- "CALL_mt1_LD_plot_D_directnmts.pdf"

#set input filenames (modify each time)
filename_pattern <- "forLD_directnmts_par1CALL"
filename_pattern2 <- "forLD_directnmts_par2CALL"
mtpar1_filename <- "mtSNP1_par1_CALL.tsv"
mtpar2_filename <- "mtSNP1_par2_CALL.tsv"

# #SAK testing
# mtpar1_filename <- "/Users/kusters/Documents/XiphoStartingOver_Aug2024/04_LDCalculation/fromOvis/LD_d-prime/CALL_d-prime/mtSNP1_par1_CALL.tsv"
# mtpar2_filename <- "/Users/kusters/Documents/XiphoStartingOver_Aug2024/04_LDCalculation/fromOvis/LD_d-prime/CALL_d-prime/mtSNP1_par2_CALL.tsv"

path_to_files <- "CALL_forovis/"
#SAK testing
# path_to_files <- "/Users/kusters/Documents/XiphoStartingOver_Aug2024/04_LDCalculation/CALL_subsetted/"

basep1fn <- list.files(pattern = filename_pattern, path = path_to_files)
basep2fn <- list.files(pattern = filename_pattern2, path = path_to_files)

p1_fn_l <- paste0(path_to_files, basep1fn)
p2_fn_l <- paste0(path_to_files, basep2fn)

out_fn_l <- paste0("mt1_CALL/mt1_results_", basep1fn)

#SAK testing
#setwd("~/Downloads")

for (filename in 1:length(p1_fn_l)) {
  #SAK testing
  #filename <- 1
  p1_fn <- p1_fn_l[filename]
  p2_fn <- p2_fn_l[filename]
  out_fn <- out_fn_l[filename]
  
  par1 <- read.csv(p1_fn, header = TRUE, sep = '\t') 
  par2 <- read.csv(p2_fn, header = TRUE, sep = '\t')
  
  mt1 <- read.csv(mtpar1_filename, header = TRUE, sep = '\t')
  mt2 <- read.csv(mtpar2_filename, header = TRUE, sep = '\t')
  
  #add mt column at end of par1 and par2
  par1 <- cbind(par1, mt1)
  par2 <- cbind(par2, mt2)
  
  par1Anc <- par1
  
  #filter bad probabilities & convert to absolute allele calls
  for (j in 1:ncol(par1)) {
    for (k in 1:length(par1[,1])) {
      if (par1[k,j] >= 0.9 & par2[k,j] <= 0.1) {
        par1Anc[k,j] <- 1.0
      } else if (par1[k,j] <= 0.1 & par2[k,j] >= 0.9) {
        par1Anc[k,j] <- 0.0
      } else if (par1[k,j] + par2[k,j] <= 0.2) {
        par1Anc[k,j] <- 0.5
      } else {
        par1Anc[k,j] <- NA
      }
    }
  }
 
  n <- length(par1[,1])
  
  if(n != length(par2[,1])) {
    stop("Input columns are not the same length")
  }
  
  #calculate LD for each AIM
  for (locus1 in 1:(ncol(par1))) { 
    
    print(paste("starting on ", locus1, sep = ""))
    
    locus2 <- ncol(par1) #is the mt locus
    
    #drop rows with missing values
    filtDat <- par1Anc[,c(locus1, locus2)] %>%
      drop_na()
    
    if (length(filtDat[,1]) == 0) {
      LD <- NA #differentiate that LD could not be calculated bc there was no data to calculate it with
      d_prime <- NA #will also have samples = 0 though
      r_sq <- NA
      
    } else {
      haplotype_freq <- sum(filtDat[,1] * filtDat[,2]) / length(filtDat[,1]) #1 is nuc locus, 2 is mt locus
      
      allele_freq_locus1 <- sum(filtDat[,1]) / length(filtDat[,1])
      allele_freq_locus2 <- sum(filtDat[,2]) / length(filtDat[,1]) #should have same # rows in it due to drop_na()
      
      if (allele_freq_locus1 == 0 | allele_freq_locus2 == 0 | allele_freq_locus1 == 1 | allele_freq_locus2 == 1) {
        LD <- NaN #differentiate that these had data but that locus was fixed
        d_prime <- NaN
        r_sq <- NaN
      } else {
        #LD calculation
        LD <- haplotype_freq - allele_freq_locus1*allele_freq_locus2
        
        #D' calculation
        if (LD > 0) {
          d_max <- min(allele_freq_locus1*(1 - allele_freq_locus2), 
                       allele_freq_locus2*(1 - allele_freq_locus1))
        } else if (LD < 0) {
          d_max <- min(allele_freq_locus1*allele_freq_locus2, 
                       (1 - allele_freq_locus1)*(1 - allele_freq_locus2))
        }
        
        d_prime <- LD/d_max
        
        #R squared calculation
        r_sq <- (LD*LD) / (allele_freq_locus1*(1 - allele_freq_locus1)*allele_freq_locus2*(1 - allele_freq_locus2))
        
      }
    }

    nom <- names(par1)[locus1]
    
    fin <- data.frame(locus = locus1, D = LD, D_prime = d_prime, r_squared = r_sq, ScaffoldID = nom, par1 = allele_freq_locus1, samples = length(filtDat[,1]))
        
    if (locus1 == 1) {
      results <- data.frame(locus = locus1, D = LD, D_prime = d_prime, r_squared = r_sq, ScaffoldID = nom, par1 = allele_freq_locus1, samples = length(filtDat[,1]))
    } else if (locus1 == ncol(par1)) {
      next
    } else{
      results <- rbind(results, fin)
    }
    
  } #end huge for
  
  write.table(results, out_fn, sep = "\t", quote = FALSE, row.names = FALSE)
} #end of giant for going through filenames



d <- results
 
p1 <- ggplot(d, aes(x = row.names(d), y = D_prime)) + geom_point() + labs(title = out_fn_ALL)
# p1
 
p2 <- ggplot(d, aes(x = row.names(d), y = D)) + geom_point() + labs(title = out_fn_ALL)
# p2

ggsave(plot1_fn, p1, width = 9, height = 5)
ggsave(plot2_fn, p2, width = 9, height = 5)

