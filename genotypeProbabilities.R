# The R version of the genotype_probabilities function initializes variables m and M with the given MAF values. 
# It then calculates the values of s using mapply, which applies the specified function 
# (creating a vector of genotype probabilities) to corresponding elements of m and M. 
# The resulting matrix s contains the genotype probabilities for each MAF value. 
# Finally, the function uses expand.grid to generate all possible combinations of genotype probabilities, 
# and apply with prod to calculate the product of probabilities for each combination. 
# The resulting vector p contains the computed genotype probabilities, which is then returned.


genotype_probabilities <- function(mafs) {
  # GENOTYPE_PROBABILITIES Compute the probabilities associated with all genotype combinations given each MAF.
  m <- mafs
  M <- 1 - m
  s <- mapply(function(m, M) c(M^2, 2 * M * m, m^2), m, M)
  p <- apply(do.call(expand.grid, s), 1, prod)
  return(p)
}
