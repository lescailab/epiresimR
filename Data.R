
# The translation assumes that the MATLAB load function loads a .mat file 
# containing the variables pts and SampleInfo. 
# The resulting R function creates a list object with the properties name, class, 
# and pts and returns it. The class property is extracted from the SampleInfo variable 
# in the MATLAB code.

Data <- function(path, snpnum) {
  name <- basename(path)
  fid <- load(path)
  NUM <- round(runif(1) * (ncol(fid$pts) - snpnum - 1))
  pts <- fid$pts[, NUM:(NUM + snpnum - 1)]
  class <- rep(0, nrow(pts))
  
  for (i in 1:nrow(pts)) {
    class[i] <- fid$SampleInfo[i, 5]
  }
  
  data <- list(name = name, class = class, pts = pts)
  class(data) <- "Data"
  return(data)
}
