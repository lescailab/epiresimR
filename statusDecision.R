# The R version of the StatusDecision function calculates the value of num 
# by iterating through the order of the model and using it to access 
# the corresponding index in SNPs and ModelInformation. 
# It then calculates R based on the penetrance values in ModelInformation. 
# Using R, it calculates UR (1 - R). 
# The function updates ProCase by multiplying it with UR and then subtracting it from 1. 
# Next, a random value a is generated using runif(1). If a is less than or equal to ProCase, 
# the function assigns Status as 1; otherwise, it assigns Status as 0. 
# Finally, the function returns Status.


StatusDecision <- function(SNPs, ModelInformation) {
  num <- 0
  for (j in 1:ModelInformation$order) {
    num <- num + (SNPs[ModelInformation$loci[j]] - 1) * 3^(ModelInformation$order - j)
  }
  R <- ModelInformation$penetrance[num + 1, 1]
  UR <- 1 - R
  ProCase <- 1
  ProCase <- ProCase * UR
  ProCase <- 1 - ProCase
  a <- runif(1)
  if (a <= ProCase) {
    Status <- 1
  } else {
    Status <- 0
  }
  return(Status)
}
