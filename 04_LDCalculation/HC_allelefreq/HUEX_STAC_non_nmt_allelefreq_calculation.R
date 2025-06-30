#Calculate Allele Frequency for Non-n-mts
#on ancestry files, par1 is X. birchmanni
#and X. cortezi is par2

#set input filenames
filename_pattern <- "forLD_nonnmts_par1HUEX" 
filename_pattern2 <- "forLD_nonnmts_par2HUEX"

path_to_files <- "HUEXSTAC_forovis/"

basep1fn <- list.files(pattern = filename_pattern, path = path_to_files)
basep2fn <- list.files(pattern = filename_pattern2, path = path_to_files)

p1_fn_l <- paste0(path_to_files, basep1fn)
p2_fn_l <- paste0(path_to_files, basep2fn)

##modified name here
out_fn_l <- paste0("Xcor_nucallelefreq_", basep1fn)
all_fn <- "nucallelefreq_forLD_allnonnmts_par1HUEXSTACdata.tsv"

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
  
  #Allele Frequency Calculation
  for (locus1 in 1:(ncol(par1))) { 
    
    print(paste("starting on ", locus1, sep = ""))
    
    #calc allele freq for the desired alleles
    l1sum = 0
    l2sum = 0
    
    for (i in 1:n) { #this is adding the header for some reason???
      l1sum <- l1sum + par1[i,locus1] + 0.5*(1-par1[i,locus1]-par2[i,locus1])
    }
    
    allele_freq_locus1 <- l1sum / n
    
    nom <- names(par1)[locus1]
    
    fin <- data.frame(ScaffoldID = nom, allele_freq = allele_freq_locus1)
    
    if (locus1 == 1) {
      results <- data.frame(ScaffoldID = nom, allele_freq = allele_freq_locus1)
    } else{
      results <- rbind(results, fin)
    }
    
  } #end allele frequency for loop
  
  write.table(results, out_fn, sep = "\t", quote = FALSE, row.names = FALSE)

} #end filename for loop

base_fn <- "nucallelefreq_forLD_nonnmts"

result_files <- list.files(path = ".", pattern = base_fn)
result_files <- paste0("./", result_files)

for (f in 1:length(result_files)) {
  #whole_fn <- paste0(base_fn, f, ".tsv")
  print(paste("adding ", f))
  result_data <- read.table(result_files[f], sep = "\t", header = TRUE)
  
  if (f ==1) {
    d <- data.frame(result_data) 
  } else {
    d <- rbind(d, result_data)
  }
}

write.table(d, all_fn, sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file = "HC_non-n-mtAlleleFreqHistogram.pdf")
histogram <- hist(result_data$allele_freq)
dev.off()
