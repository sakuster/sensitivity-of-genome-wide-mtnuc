###Partial correlation approach -- non-n-mts
###Shady Kuster
###16 July 2025

library(ggplot2)
library(dplyr)
library(tidyr)
library(ppcor)

#set output filenames (modify each time)
out_fn_ALL <- "mt1_CALL/mt1_results_nonnmt_CALL_ALL.tsv"
plot1_fn <- "CALLmt1_LD_plot_rSquared_nonnmts.pdf"
#plot2_fn <- "CALLmt1_LD_plot_D_nonnmts.pdf"

#set input filenames (modify each time)
filename_pattern <- "forLD_nonnmts_par1CALL"
filename_pattern2 <- "forLD_nonnmts_par2CALL"
mtpar1_filename <- "mtSNP1_par1_CALL.tsv"
mtpar2_filename <- "mtSNP1_par2_CALL.tsv"
nuc_filename <- "CALL_ancestryResults_v2.tsv"

path_to_files <- "CALL_forovis/"

basep1fn <- list.files(pattern = filename_pattern, path = path_to_files)
basep2fn <- list.files(pattern = filename_pattern2, path = path_to_files)

p1_fn_l <- paste0(path_to_files, basep1fn)
p2_fn_l <- paste0(path_to_files, basep2fn)

#ALSO MUST MODIFY EACH TIME
out_fn_l <- paste0("mt1_CALL/mt1_results_", basep1fn)

for (filename in 1:length(p1_fn_l)) {
  #sak testing
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
  
  #read in nuclear ancestry values
  nucAnc <- read.csv(nuc_filename, 
                     header = TRUE, sep = "\t") %>%
    head(-1) %>% #remove text about population mean
    dplyr::select(Percent.Par1.Ancestry) 
  par1Anc <- cbind(par1Anc, nucAnc)
 
  n <- length(par1[,1])
  
  if(n != length(par2[,1])) {
    stop("Input columns are not the same length")
  }
  
  #Dan's suggestion
  for (locus1 in 1:(ncol(par1) - 1)) { #- 1)) { #992 is the mt locus
    #SAK TESTING
    #locus1 <- 1
    
    print(paste("starting on ", locus1, sep = ""))
  
    #locus1 <- 2 #remember that col1 is the fish IDs so this cannot be 1
    locus2 <- ncol(par1) #is 992, the mt locus
    
    locus2 <- which(colnames(par1Anc) == "xbir.mito.45")
    
    nucAncCol <- which(colnames(par1Anc) == "Percent.Par1.Ancestry")
    
    testDat <- par1Anc[,c(locus1, locus2, nucAncCol)] %>%
      drop_na()
    
    if (nrow(testDat) > 3) {
      parCor <- pcor.test(testDat[,1], 
                          testDat[,2], 
                          testDat[,3])$estimate #will give you estimate of partial correlation between 1 nuclear locus for all fish and the mtDNA locus 
      
    } else {
      message(paste(locus1, "from subset", filename, 
                    "thrown out. You need more than 3 rows to perform partial correlation test."))
      parCor <- NA
    }
    
    
    nom <- names(par1)[locus1]
    
    #fin <- data.frame(locus = locus1, D = LD, D_prime = d_prime, ScaffoldID = nom)
    fin <- data.frame(locus = locus1, pcor = parCor, ScaffoldID = nom, samples = length(testDat[,1]))
    
    if (locus1 == 1) {
      results <- data.frame(locus = locus1, pcor = parCor, ScaffoldID = nom, samples = length(testDat[,1]))
    } else if (locus1 == ncol(par1)) {
      next
    } else{
      results <- rbind(results, fin)
    }
    
  } #end huge for
  
  write.table(results, out_fn, sep = "\t", quote = FALSE, row.names = FALSE)

} #end of giant for going through filenames

###then need to read in all 3 and combine?
base_fn <- "mt1_results_forLD_nonnmts"

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

