### HUEX-STAC direct n-mt allele frequency calculation

#set input filenames (modify each time)
filename_pattern <- "forLD_directnmts_par2HUEX"
filename_pattern2 <- "forLD_directnmts_par1HUEX"

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

  n <- length(par1[,1])
  
  if(n != length(par2[,1])) {
    stop("Input columns are not the same length")
  }
  
  #calculate allele frequency
  for (locus1 in 1:(ncol(par1))) { 
    
    print(paste("starting on ", locus1, sep = ""))
    
    #set 0 each time
    l1sum = 0
    l2sum = 0
    
    for (i in 1:n) {
      l1sum <- l1sum + par1[i,locus1] + 0.5*(1-par1[i,locus1]-par2[i,locus1])
    }
    
    allele_freq_locus1 <- l1sum / n
    
    nom <- names(par1)[locus1]
    
    fin <- data.frame(ScaffoldID = nom, allele_freq = allele_freq_locus1)
    
    #add results to dataframe
    if (locus1 == 1) {
      results <- data.frame(ScaffoldID = nom, allele_freq = allele_freq_locus1)
    } else {
      results <- rbind(results, fin)
    }
    
    print(paste("finished", locus1))
    
  } #end allele frequency for loop
  
  write.table(results, out_fn, sep = "\t", quote = FALSE, row.names = FALSE)
  
} #end loop through file names

hist(results$allele_freq) 
mean(results$allele_freq)
