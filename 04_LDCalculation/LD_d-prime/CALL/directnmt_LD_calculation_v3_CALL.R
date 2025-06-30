###Calculate LD on Xiphophorus data
###Shady Kuster
###21 July 2022
library(ggplot2)

#set output filenames (modify each time)
out_fn_ALL <- "mt1_results_directnmt_CALL_ALL.tsv"
plot1_fn <- "CALL_mt1_LD_plot_DPRIME_directnmts.pdf"
plot2_fn <- "CALL_mt1_LD_plot_D_directnmts.pdf"

#set input filenames (modify each time)
filename_pattern <- "forLD_directnmts_par1CALL"
filename_pattern2 <- "forLD_directnmts_par2CALL"
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
 
  n <- length(par1[,1])
  
  if(n != length(par2[,1])) {
    stop("Input columns are not the same length")
  }
  
  #Dan's suggestion
  for (locus1 in 1:(ncol(par1))) { #- 1)) { #992 is the mt locus
    
    print(paste("starting on ", locus1, sep = ""))
  
    #locus1 <- 2 #remember that col1 is the fish IDs so this cannot be 1
    locus2 <- ncol(par1) #is 992, the mt locus
    hap_sum <- 0
    
    for (i in 1:n) { 
      homo1prob <- par1[i,locus1]
      homo2prob <- par2[i,locus1]
      heteroprob <- 1 - homo1prob - homo2prob
      
      homo1prob_locus2 <- par1[i,locus2]
      homo2prob_locus2 <- par2[i,locus2]
      heteroprob_locus2 <- 1 - homo1prob_locus2 - homo2prob_locus2
      
      hap_sum <- hap_sum + ((homo1prob + (0.5 * heteroprob)) * 
        (homo1prob_locus2 + (0.5 * heteroprob_locus2)))
    }
    
    haplotype_freq <- hap_sum/n
    
    #calc allele freq for the desired alleles
    l1sum = 0
    l2sum = 0
    
    for (i in 1:n) { #this is adding the header for some reason???
      l1sum <- l1sum + par1[i,locus1] + 0.5*(1-par1[i,locus1]-par2[i,locus1])
      l2sum <- l2sum + par1[i,locus2] + 0.5*(1-par1[i,locus2]-par2[i,locus2])
    }
    
    allele_freq_locus1 <- l1sum / n
    allele_freq_locus2 <- l2sum / n
                       
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

    nom <- names(par1)[locus1]
    
    fin <- data.frame(locus = locus1, D = LD, D_prime = d_prime, ScaffoldID = nom)
        
    if (locus1 == 1) {
      results <- data.frame(locus = locus1, D = LD, D_prime = d_prime, ScaffoldID = nom)
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
 
p2 <- ggplot(d, aes(x = row.names(d), y = D)) + geom_point() + labs(title = out_fn_ALL)

ggsave(plot1_fn, p1, width = 9, height = 5)
ggsave(plot2_fn, p2, width = 9, height = 5)

