EpiReSIM <- function(Case_Num, Control_Num, SNP_num, MAF, PD, H2, n, repeat, SNP_name, data_mat, data_txt, realdata_name) {
  # SNP_num: the number of SNPs
  # Case_Num: the number of case samples
  # Control_Num: the number of control samples
  # MAF: the Minor Allele Frequency
  # PD: the value of prevalence of model
  # H2: the value of heritability of model. If user calculates model with heritability,
  # then user should give the value of heritability, otherwise the values should be 0.
  # n: the order of the model
  # repeat: the number of simulated dataset
  
  a <- Data(realdata_name, SNP_num)
  pt$order <- n
  
  ## Select the control dataset and calculate the MAF of them
  idx <- which(all(a$class[1,] == 2, na.rm = TRUE))
  # idx <- which(all(a$class[1,] == 1, na.rm = TRUE))
  d <- a$pts[idx, ]
  maf <- rep(0, ncol(d))
  
  for (j in 1:ncol(d)) {
    Aa <- 0
    aa <- 0
    
    for (i in 1:nrow(d)) {
      if (d[i, j] == 2) {
        Aa <- Aa + 1
      }
      if (d[i, j] == 3) {
        aa <- aa + 1
      }
    }
    
    maf[j] <- (Aa + aa * 2) / (nrow(d) * 2)
  }
  
  ## Select the site of the model
  site <- integer(n)
  
  for (j in 1:n) {
    s <- 0.01 # Search in steps of 0.01
    candidates <- integer(0)
    
    while (TRUE) {
      for (i in 1:ncol(d)) {
        if (maf[i] >= (MAF[j] - s) && maf[i] <= (MAF[j] + s)) {
          candidates <- c(candidates, i)
        }
      }
      
      if (length(candidates) > 0) {
        break
      } else {
        s <- s + 0.01
      }
    }
    
    site[j] <- sample(candidates, 1)
  }
  
  pt$MAF <- maf[site]
  pt$loci <- site[1]
  
  ## Calculate model using prevalence or heritability
  pt$penetrance <- calculation(pt$MAF, PD, H2, n)$penetrance
  pt$pd <- calculation(pt$MAF, PD, H2, n)$pd
  pt$h <- calculation(pt$MAF, PD, H2, n)$h
  
  ## Generate simulated data and calculate sample labels
  for (k in 1:repeat) {
    SNP <- simulation(d, Control_Num, Case_Num, pt, SNP_num)
    
    if (data_mat) {
      filename1 <- paste0(SNP_name, "_", k, ".mat")
      save(SNP, file = filename1)
    }
    
        if (data_txt) {
      filename2 <- paste0(SNP_name, "_", k, ".txt")
      write.table(SNP, file = filename2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
  
  ## Save model information
  fid <- file("log.txt", "w")
  cat("order:", pt$order, "\n", file = fid)
  cat("MAF:", paste(pt$MAF, collapse = "\t"), "\n\n", file = fid)
  cat("loci:", paste(pt$loci, collapse = "\t"), "\n\n", file = fid)
  cat("Penetrance:\n", file = fid)
  
  for (i in 1:(3^n)) {
    cat(pt$penetrance[i], "\t", file = fid)
    
    if (i %% 3 == 0) {
      cat("\n", file = fid)
    }
    
    if (i %% 9 == 0) {
      cat("\n", file = fid)
    }
  }
  
  cat("\n", file = fid)
  cat("Prevalence:", pt$pd, "\n", file = fid)
  cat("Heritability:", pt$h, "\n", file = fid)
  close(fid)
}

