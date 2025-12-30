### HUEX-STAC direct n-mt allele frequency calculation

library(dplyr)
library(tidyr)

#set input filenames (modify each time)
filename_pattern <- "forLD_directnmts_par1HUEX"
filename_pattern2 <- "forLD_directnmts_par2HUEX"

path_to_files <- "./HUEXSTAC_forovis"

basep1fn <- list.files(pattern = filename_pattern, path = path_to_files)
basep2fn <- list.files(pattern = filename_pattern2, path = path_to_files)

p1_fn_l <- paste0(path_to_files, "/", basep1fn)
p2_fn_l <- paste0(path_to_files, "/", basep2fn)

##modified name here
out_fn_l <- paste0("nucallelefreq_", basep1fn)

#loop through all input files
for (filename in 1:length(p1_fn_l)) {
  p1_fn <- p1_fn_l[filename]
  p2_fn <- p2_fn_l[filename]
  out_fn <- out_fn_l[filename]

  par1 <- read.csv(p1_fn, header = TRUE, sep = '\t') 
  par2 <- read.csv(p2_fn, header = TRUE, sep = '\t')
  
  par1Anc <- par1
  par2Anc <- par2
  
  #filter bad probabilities & convert to absolute allele calls
  for (j in 1:ncol(par1)) {
    for (k in 1:length(par1[,1])) {
      if (par1[k,j] >= 0.9 & par2[k,j] <= 0.1) {
        par1Anc[k,j] <- 1.0
        par2Anc[k,j] <- 0.0
      } else if (par1[k,j] <= 0.1 & par2[k,j] >= 0.9) {
        par1Anc[k,j] <- 0.0
        par2Anc[k,j] <- 1.0
      } else if (par1[k,j] + par2[k,j] <= 0.2) {
        par1Anc[k,j] <- 0.5
        par2Anc[k,j] <- 0.5
      } else {
        par1Anc[k,j] <- NA
        par2Anc[k,j] <- NA
      }
    }
  }

  n <- length(par1[,1])
  
  if(n != length(par2[,1])) {
    stop("Input columns are not the same length")
  }
  
  #calculate allele frequency
  for (locus1 in 1:(ncol(par1))) { 
    
    print(paste("starting on ", locus1, sep = ""))
    
    #set 0 each time
    # l1sum = 0
    # l2sum = 0
    # 
    #rem = 0
    
    #remove any rows w/ missing allele calls
    filtDat <- par1Anc[,c(locus1)] %>%
      cbind(par2Anc[,c(locus1)]) %>%
      as.data.frame() %>%
      drop_na()
    probabDat <- par1[,c(locus1)] %>% 
      cbind(par2[,c(locus1)]) %>%
      as.data.frame() %>%
      drop_na()
    
    #allele freqs for allele calls
    filtAlleFreq1 <- sum(filtDat[,1]) / length(filtDat[,1])
    filtAlleFreq2 <- sum(filtDat[,2]) / length(filtDat[,2])
    
    #allele freqs using probabilities only 
    probAlleFreq1 <- sum(probabDat[,1]) / length(probabDat[,1])
    probAlleFreq2 <- sum(probabDat[,2]) / length(probabDat[,2])
    
    #allele freqs using probabilities AND old calculation 
    probAlleFreqCalc1 <- sum(probabDat[,1] + 0.5*(1-probabDat[,1]-probabDat[,2])) / length(probabDat[,1])
    probAlleFreqCalc2 <- sum(probabDat[,2] + 0.5*(1-probabDat[,2]-probabDat[,1])) / length(probabDat[,2])
    
    nom <- names(par1)[locus1]
    
    fin <- data.frame(ScaffoldID = nom, 
                      par1allele_freq = filtAlleFreq1, 
                      par1samples = length(filtDat[,1]),
                      par2allele_freq = filtAlleFreq2,
                      par2samples = length(filtDat[,2]),
                      par1ProbAF = probAlleFreqCalc1,
                      par2ProbAF = probAlleFreqCalc2,
                      par1Wrong = probAlleFreq1,
                      par2Wrong = probAlleFreq2,
                      nProbSamples = length(probabDat[,1]),
                      nProbSamplespar2 = length(probabDat[,2]))
    
    #add results to dataframe
    if (locus1 == 1) {
      results <- data.frame(ScaffoldID = nom, 
                            par1allele_freq = filtAlleFreq1, 
                            par1samples = length(filtDat[,1]),
                            par2allele_freq = filtAlleFreq2,
                            par2samples = length(filtDat[,2]),
                            par1ProbAF = probAlleFreqCalc1,
                            par2ProbAF = probAlleFreqCalc2,
                            par1Wrong = probAlleFreq1,
                            par2Wrong = probAlleFreq2,
                            nProbSamples = length(probabDat[,1]),
                            nProbSamplespar2 = length(probabDat[,2]))
    } else {
      results <- rbind(results, fin)
    }
    
    print(paste("finished", locus1))
    
  } #end allele frequency for loop
  
  write.table(results, out_fn, sep = "\t", quote = FALSE, row.names = FALSE)
  
  #write results per file to a table
  if (!dir.exists("alleleCalls")) {
    dir.create("alleleCalls")
  }
  
  fnShort1 <- gsub(".tsv", "", basename(p1_fn_l[filename]))
  fnShort2 <- gsub(".tsv", "", basename(p2_fn_l[filename]))
  write.table(par1Anc, paste0("alleleCalls/", fnShort1, "alleleCalls_par1.tsv"),
              quote = F, row.names = F, sep = "\t")
  write.table(par2Anc, paste0("alleleCalls/", fnShort2, "alleleCalls_par1.tsv"),
              quote = F, row.names = F, sep = "\t")
  
} #end loop through file names

