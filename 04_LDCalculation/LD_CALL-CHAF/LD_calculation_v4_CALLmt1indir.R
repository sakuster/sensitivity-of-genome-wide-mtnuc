###Calculate LD on Xiphophorus data
###Shady Kuster
###21 July 2022

###this script currently is for use of smaller sets of data
###eventually, will be translated into a bash script or something
###else bc R can't handle the data I have

#this has been modified for ovis

library(tidyr)

#set output filenames (modify each time)
out_fn_ALL <- "mt1_CALL/mt1_results_indirectnmt_CALL_ALL.tsv"
plot1_fn <- "CALLmt1_LD_plot_DPRIME_indirectnmts.pdf"
plot2_fn <- "CALLmt1_LD_plot_D_indirectnmts.pdf"

#set input filenames (modify each time)
filename_pattern <- "forLD_indirectnmts_par1CALL"
filename_pattern2 <- "forLD_indirectnmts_par2CALL"
mtpar1_filename <- "mtSNP1_par1_CALL.tsv"
mtpar2_filename <- "mtSNP1_par2_CALL.tsv"

path_to_files <- "CALL_forovis/"

basep1fn <- list.files(pattern = filename_pattern, path = path_to_files)
basep2fn <- list.files(pattern = filename_pattern2, path = path_to_files)

p1_fn_l <- paste0(path_to_files, basep1fn)
p2_fn_l <- paste0(path_to_files, basep2fn)

out_fn_l <- paste0("mt1_CALL/mt1_results_", basep1fn)

for (filename in 1:length(p1_fn_l)) {
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
  
  #calculate LD
  for (locus1 in 1:(ncol(par1))) { #- 1)) { #992 is the mt locus
    
    print(paste("starting on ", locus1, sep = ""))
  
    #locus1 <- 2 #remember that col1 is the fish IDs so this cannot be 1
    locus2 <- ncol(par1) #is 992, the mt locus
    
    #remove any rows w/ missing allele calls
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
  #can also just append to OG file here
  #YOU ARE INCLUDING THE MT LOCI 
} #end of giant for going through filenames

###then need to read in all 3 and combine?
base_fn <- "mt1_results_forLD_indirectnmt"

result_files <- list.files(path = "mt1_CALL", pattern = base_fn)
result_files <- paste0("mt1_CALL/", result_files)

for (f in 1:length(result_files)) {
  #whole_fn <- paste0(base_fn, f, ".tsv")
  print(paste("adding ", f))
  result_data <- read.table(result_files[f], sep = "\t", header = TRUE)
  
  if (f ==1) {
    d <- data.frame(result_data) #need to add which subset it's from so the locus means something
  } else {
    d <- rbind(d, result_data)
  }
}

write.table(d, out_fn_ALL,sep = "\t", quote = FALSE, row.names = FALSE)

library(ggplot2)
 
# d <- results
 
p1 <- ggplot(d, aes(x = row.names(d), y = D_prime)) + geom_point() + labs(title = out_fn_ALL)
# 
p2 <- ggplot(d, aes(x = row.names(d), y = D)) + geom_point() + labs(title = out_fn_ALL)

# 
ggsave(plot1_fn, p1, width = 9, height = 5)
ggsave(plot2_fn, p2, width = 9, height = 5)
# 
# library(dplyr)
# 
# #dir <- d[1:length(d$locus) -1,] %>% mutate(class = "direct_n-mt")
# dir <- read.table("results_directnmt_CALLredo.tsv", header = TRUE, sep = "\t", nrows = 991) %>% mutate(class = "direct_n-mt")
# indir <- read.table("results_indirectnmt_CALL_ALL.tsv", header = TRUE, sep = "\t") %>% mutate(class = "indirect_n-mt")
# non <- read.table("results_nonnmt_CALL_ALL.tsv", header = TRUE, sep = "\t") %>% mutate(class = "non_n-mt")
# 
# #accidentally included the mt locus in the final output, removing it
# indir <- filter(indir, D_prime < 0.8)
# non <- filter(non, D_prime < 0.8) #the ones above 0.8 are the mt locus
# 
# 
# dat <- rbind(dir, indir, non)
# 
# all_plot <- ggplot(dat, aes(x = class, y = D_prime, color = class)) + 
#   geom_boxplot() + theme_classic() + labs(title = "LD for one mt locus")
# all_plot
# 
# 
# ggsave("One_mt_all_LD.pdf",all_plot, width = 7, height = 4)
