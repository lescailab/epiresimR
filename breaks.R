# The R version of the get_breaks function generates a random number of breakpoints (num) 
# based on the number of columns in the input matrix d. It then creates a vector brk 
# of random values between 1 and the number of columns in d. The vector is sorted in ascending order, 
# and duplicate values are removed using the unique function. 
# Any value equal to 0 is excluded from the vector. 
# Finally, the function appends Inf to the end of the vector 
# as the last breakpoint and returns the resulting vector brk.

get_breaks <- function(d) {
  # Randomly generate breakpoints
  # Returns a vector: a vector of integers showing the location of data breakpoints for each sample
  num <- round(runif(1) * ncol(d))
  brk <- sort(round(runif(num) * ncol(d)))
  brk <- unique(brk)
  brk <- brk[brk != 0]
  brk <- c(brk, Inf)
  return(brk)
}
