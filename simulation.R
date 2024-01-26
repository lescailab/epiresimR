# The translation assumes that the functions get_breaks, Splicing, and StatusDecision are defined elsewhere in your code. 
# The R version of the simulation function initializes an empty matrix SNP with dimensions inds by SNP_num + 1. 
# It then loops through the required number of iterations, generates breakpoints, performs splicing, and assigns 
# the resulting SNP values to the SNP matrix. The StatusDecision function is used to determine the status (case or control) 
# based on the calculated values. The status value is stored in the last column of the SNP matrix. 
# Finally, the function returns the SNP matrix.


simulation <- function(D, control_num, case_num, pt, SNP_num) {
  SNP_num <- SNP_num
  Case_Num <- case_num
  Control_Num <- control_num
  inds <- Case_Num + Control_Num
  Current_Case_Num <- 0
  Current_Control_Num <- 0
  SNP <- matrix(0, nrow = inds, ncol = SNP_num + 1)
  
  i <- 0
  while (i < inds) {
    brk <- get_breaks(D)
    i <- i + 1
    snp <- Splicing(brk, D)
    
    for (j in 1:SNP_num) {
      SNP[i, j] <- snp[j]
    }
    
    Status <- StatusDecision(SNP[i, ], pt)
    
    if (Status == 1) {  # case
      if (Current_Case_Num < Case_Num) {
        SNP[i, SNP_num + 1] <- Status
        Current_Case_Num <- Current_Case_Num + 1
      } else {
        i <- i - 2
      }
    } else {  # control
      if (Current_Control_Num < Control_Num) {
        SNP[i, SNP_num + 1] <- Status
        Current_Control_Num <- Current_Control_Num + 1
      } else {
        i <- i - 1
      }
    }
  }
  
  return(SNP)
}
