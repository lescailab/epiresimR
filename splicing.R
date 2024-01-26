# The R version of the Splicing function creates an empty numeric vector SNP 
# with the same length as the number of columns in the input matrix d. 
# It initializes variables brk_num and column to keep track of the breakpoints 
# and the current column index. 
# A random index j is generated between 1 and the number of rows in d. 
# The function then iterates through the columns of d and checks if 
# the current column matches the breakpoint. 
# If it does not match, the value at index j in the column is assigned to 
# the corresponding index in SNP, and the column index is incremented. 
# If the column matches the breakpoint, a new random index j is generated (excluding 0), 
# and the value at index j in the column is assigned to the corresponding index in SNP. 
# Finally, the function returns the resulting vector SNP.


Splicing <- function(brk, d) {
  # Splice data
  # Returns a vector: a vector of spliced data based on breakpoints
  
  brk_num <- 1
  column <- 1
  j <- round(runif(1) * (nrow(d) - 2)) + 1
  
  SNP <- numeric(length = ncol(d))
  
  while (column <= ncol(d)) {
    if (column != brk[brk_num]) {
      SNP[column] <- d[j, column]
      column <- column + 1
    } else {
      j <- round(runif(1) * (nrow(d) - 1))
      
      while (j == 0) {
        j <- round(runif(1) * (nrow(d) - 1))
      }
      
      SNP[column] <- d[j, column]
      column <- column + 1
      brk_num <- brk_num + 1
    }
  }
  
  return(SNP)
}
