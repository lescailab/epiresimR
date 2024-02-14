#!/usr/bin/env R

library(readr)
library(dplyr)
library(data.table)


vcf_to_mat <- function(vcf_file, ref_file) {

#use readVcf()
  vcf <- read_delim(vcf_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)
  
  # replace genotypes: '0|0'-> 1 , '1|0'/'0|1'-> 2, '1|1'-> 3
  sample_cols <- 10:ncol(vcf)
  for (col in sample_cols) {
    vcf[,col] <- as.numeric(factor(vcf[,col], levels = c('0|0', '0|1', '1|0', '1|1'), labels = c(1, 2, 2, 3)))
  }
  
  # convert to matrix
  pts <- as.matrix(vcf[, sample_cols])
  pts <- as.integer(pts)
  
  # Sample ids
  iid <- colnames(vcf)[sample_cols]
  

  
  # here we create a PED like file, take into consideration the use of a real or already available ped  
  fid <- paste0("FAM_", 1:numRows)
  numRows <- length(iid)
  PHENOTYPE <- c(rep(1, numRows/2), rep(2, numRows/2))
  pedtable <- data.frame(IID = iid, FID = fid, col3 = rep(0, numRows), col4 = rep(0, numRows), PHENOTYPE = PHENOTYPE, col6 = rep(0, numRows))
  
  # create the Structure of TD1.mat input from authors. 

  #  %  struct with fields:
  #  %       SNPInfo: {18268x8 cell}
  #  %    SampleInfo: {4062x6 cell}
  #  %       Counter: [5x18268 double]
  #  %           pts: [4062x18267 int8]

  ref <- list()
  ref$SNPInfo <- matrix(nrow = 18268, ncol = 8) # Adjust dimensions as necessary
  ref$Counter <- matrix(0, nrow = 5, ncol = 18268)
  ref$SampleInfo <- pedtable
  ref$pts <- t(pts) # Transpose to match MATLAB structure
  
  # this is the input to the main function epiReSim.R
  saveRDS(ref, file = ref_file)
}



