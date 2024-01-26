
# Please note that the code snippet above is a direct translation of the provided MATLAB code to R. 
# However, please note that the lsqminnorm function used in the code snippet is not a standard R function. 
# You would need to replace it with an appropriate function for solving least squares problems in R, 
# such as lm or other relevant functions depending on your specific requirements.

# Function to perform calculations
calculation <- function(MAF, PD, H2, n) {
  k <- n - 1
  X <- numeric(3^n)
  X <- matrix(X, nrow = 3^n, ncol = 1)
  D <- matrix(genotype_probabilities(MAF), nrow = 1)


# Please note that some MATLAB-specific functions and operations may not have direct equivalents in R. 
# You will need to replace them with appropriate R functions or implement the equivalent functionality yourself. 
# Additionally, certain functions like lsqminnorm, jacobian, and pre_seidel are not standard R functions, 
# so you will need to replace them with the appropriate R functions or custom implementations.
# Make sure to carefully review and test the translated code to ensure correctness and compatibility with your specific requirements and data.


if (n == 2) {
  # the marginal penetrance values of A
  PD_AA <- 0
  PD_Aa <- 0
  PD_aa <- 0
  PD_BB <- 0
  PD_Bb <- 0
  PD_bb <- 0
  for (i in 1:(3^k)) {
    PD_AA <- PD_AA + X[i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
    PD_Aa <- PD_Aa + X[3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
    PD_aa <- PD_aa + X[2 * 3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
  }
  PD_AA <- PD_AA - PD
  PD_Aa <- PD_Aa - PD
  PD_aa <- PD_aa - PD
  # the marginal penetrance values of A
  for (j in 1:(3^k)) {
    i <- 3 * (j - 1) + 1
    PD_BB <- PD_BB + X[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    PD_Bb <- PD_Bb + X[1 + i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    PD_bb <- PD_bb + X[2 + i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
  }
  PD_BB <- PD_BB - PD
  PD_Bb <- PD_Bb - PD
  PD_bb <- PD_bb - PD
  f10 <- D %*% X - PD
  if (H2 == 0) {
    # using the prevalence to calculate the eNME model
    E <- matrix(0, nrow = 7, ncol = 3^n)
    f <- matrix(0, nrow = 7, ncol = 1)
    E[1, ] <- PD_AA
    E[2, ] <- PD_Aa
    E[3, ] <- PD_aa
    E[4, ] <- PD_BB
    E[5, ] <- PD_Bb
    E[6, ] <- PD_bb
    E[7, ] <- f10
    f <- rbind(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, f10)
    x_result <- lsqminnorm(E, f)
  } else {
    # When using both prevalence and heritability to calculate the eNME model
    f11 <- 0
    for (i in 1:(3^n)) {
      f11 <- f11 + D[i] * ((X[i] - PD)^2)
    }
    f11 <- f11 - H2 * PD * (1 - PD)
    f <- matrix(0, nrow = 8, ncol = 1)
    f[1, ] <- PD_AA
    f[2, ] <- PD_Aa
    f[3, ] <- PD_aa
    f[4, ] <- PD_BB
    f[5, ] <- PD_Bb
    f[6, ] <- PD_bb
    f[7, ] <- f10
    f[8, ] <- f11
    x0 <- matrix(0, nrow = 3^n, ncol = 1)
    error_dxk <- 0.01
    error_fkk <- 0.05
    num <- 10
    jacobi <- jacobian(f, X)  # Replace with appropriate R function
    for (k in 1:num) {
      Ak <- jacobi  # Replace with appropriate R code to evaluate Jacobian
      bk <- f  # Replace with appropriate R code to evaluate f
      dxk <- pre_seidel(Ak, -bk, k)  # Replace with appropriate R function
      x0 <- x0 + dxk
      fkk <- f  # Replace with appropriate R code to evaluate f
      if (norm(dxk) < error_dxk || norm(fkk) < error_fkk) {
        break
      }
    }
    
    if (k < num) {
      x_result <- x0
    } else {
      x_result <- x0
    }
  }
  for (i in 1:(3^n)) {
    if (x_result[i] < 0) {
      x_result[i] <- 0
    }
  }
}

# Please note that some variables and functions in the code are not defined, 
# such as D, pre_seidel, jacobian, and norm. You will need to define or replace them 
# with the appropriate values or functions in your implementation. 
# Additionally, make sure to review and test the translated code to ensure correctness 
# and compatibility with your specific requirements and data.

if (n == 3) {
    # Calculate the marginal penetrance of the A locus alleles
    PD_AA <- 0
    PD_Aa <- 0
    PD_aa <- 0
    for (i in 1:(3^k)) {
        PD_AA <- PD_AA + X[i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
        PD_Aa <- PD_Aa + X[3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
        PD_aa <- PD_aa + X[2 * 3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
    }
    PD_AA <- PD_AA - PD
    PD_Aa <- PD_Aa - PD
    PD_aa <- PD_aa - PD
    
    # Calculate the marginal penetrance of the B locus alleles
    PD_BB <- 0
    PD_Bb <- 0
    PD_bb <- 0
    for (o in 0:2) {
        for (j in 1:3^(k-1)) {
            i <- j + o * (3^k)
            PD_BB <- PD_BB + X[i, 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2 * 3^(k-1) + i])
            PD_Bb <- PD_Bb + X[i + 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2 * 3^(k-1) + i])
            PD_bb <- PD_bb + X[i + 2 * 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2 * 3^(k-1) + i])
        }
    }
    PD_BB <- PD_BB - PD
    PD_Bb <- PD_Bb - PD
    PD_bb <- PD_bb - PD
    
    # Calculate the allele frequencies of the C locus
    PD_CC <- 0
    PD_Cc <- 0
    PD_cc <- 0
    for (o in 0:2) {
        for (j in 0:k) {
            i <- 1 + o * 3^k + 3 * j
            PD_CC <- PD_CC + X[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
            PD_Cc <- PD_Cc + X[i + 1, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
            PD_cc <- PD_cc + X[i + 2, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
        }
    }
    PD_CC <- PD_CC - PD
    PD_Cc <- PD_Cc - PD
    PD_cc <- PD_cc - PD
    
    f10 <- D %*% X - PD
    if (H2 == 0) {
        A <- X
        b <- c(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, f10)
        x_result <- solve(t(A) %*% A, t(A) %*% b)
    } else {
        f11 <- 0
        for (i in 1:(3^n)) {
            f11 <- f11 + (X[i] - PD)^2 * D[i]
        }
        f11 <- f11 - H2 * PD * (1 - PD)
        f <- c(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, f10, f11)
        x0 <- rep(0, 3^n)
        error_dxk <- 0.05
        error_fkk <- 0.05
        num <- 10
        jacobi <- jacobian(f, X)
        for (k in 1:num) {
            Ak <- jacobi
            bk <- f
            dxk <- pre_seidel(Ak, -bk, k)
            x0 <- x0 + dxk
            fkk <- f
            if (norm(dxk) < error_dxk || norm(fkk) < error_fkk) {
                break
            }
        }
        if (k < num) {
            x_result <- x0
        } else {
            x_result <- x0
        }
    }
    for (i in 1:(3^n)) {
        if (x_result[i] < 0) {
            x_result[i] <- 0
        }
    }
}


# Please note that some parts of the code are missing or incomplete,
# such as the calculation of the Jacobian matrix and the pre_seidel function,
# which need to be defined or replaced with appropriate code.


if (n == 4) {
  # Calculation of marginal penetrance of allele A
  PD_AA <- 0
  PD_Aa <- 0
  PD_aa <- 0
  for (i in 1:(3^k)) {
    PD_AA <- PD_AA + X[i, 1] * (D[1, i] + D[1, (3^k) + i] + D[1, 2 * (3^k) + i])
    PD_Aa <- PD_Aa + X[(3^k) + i, 1] * (D[1, i] + D[1, (3^k) + i] + D[1, 2 * (3^k) + i])
    PD_aa <- PD_aa + X[2 * (3^k) + i, 1] * (D[1, i] + D[1, (3^k) + i] + D[1, 2 * (3^k) + i])
  }
  PD_AA <- PD_AA - PD
  PD_Aa <- PD_Aa - PD
  PD_aa <- PD_aa - PD
  
  # Calculation of marginal penetrance of allele B
  PD_BB <- 0
  PD_Bb <- 0
  PD_bb <- 0
  for (o in 0:2) {
    for (j in 1:(3^(k-1))) {
      i <- j + o * (3^k)
      PD_BB <- PD_BB + X[i, 1] * (D[1, i] + D[1, (3^(k-1)) + i] + D[1, 2 * (3^(k-1)) + i])
      PD_Bb <- PD_Bb + X[i + (3^(k-1)), 1] * (D[1, i] + D[1, (3^(k-1)) + i] + D[1, 2 * (3^(k-1)) + i])
      PD_bb <- PD_bb + X[i + 2 * (3^(k-1)), 1] * (D[1, i] + D[1, (3^(k-1)) + i] + D[1, 2 * (3^(k-1)) + i])
    }
  }
  PD_BB <- PD_BB - PD
  PD_Bb <- PD_Bb - PD
  PD_bb <- PD_bb - PD
  
  # Calculation of marginal penetrance of allele C
  PD_CC <- 0
  PD_Cc <- 0
  PD_cc <- 0
  for (o in 0:(3^(k-1)-1)) {
    for (j in 1:3) {
      i <- j + (3^(n-2)) * o
      PD_CC <- PD_CC + X[i, 1] * (D[1, i] + D[1, (3^1) + i] + D[1, 2 * (3^1) + i])
      PD_Cc <- PD_Cc + X[i + 3 * 1, 1] * (D[1, i] + D[1, (3^1) + i] + D[1, 2 * (3^1) + i])
      PD_cc <- PD_cc + X[i + 3 * 2, 1] * (D[1, i] + D[1, (3^1) + i] + D[1, 2 * (3^1) + i])
    }
  }
  PD_CC <- PD_CC - PD
  PD_Cc <- PD_Cc - PD
  PD_cc <- PD_cc - PD
  
  # Calculation of marginal penetrance of allele D
  PD_DD <- 0
  PD_Dd <- 0
  PD_dd <- 0
  for (o in 0:2) {
    for (j in 0:k) {
      i <- 1 + o * (3^k) + 3 * j
      PD_DD <- PD_DD + X[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
      PD_Dd <- PD_Dd + X[i + 1, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
      PD_dd <- PD_dd + X[i + 2, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    }
  }
  PD_DD <- PD_DD - PD
  PD_Dd <- PD_Dd - PD
  PD_dd <- PD_dd - PD
  
  f10 <- D %*% X - PD
  
  if (H2 == 0) {
    # Not using genetic force calculation, converting to coefficient matrix
    E <- cbind(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, PD_DD, PD_Dd, PD_dd, f10)
    f <- c(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, PD_DD, PD_Dd, PD_dd, f10)
    x_result <- solve(t(E) %*% E) %*% t(E) %*% f
  } else {
    # Using disease prevalence calculation with Newton's method
    f11 <- 0
    for (i in 1:(3^n)) {
      f11 <- f11 + (X[i] - PD)^2 * D[i]
    }
    f11 <- f11 - H2 * PD * (1 - PD)
    f <- c(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, PD_DD, PD_Dd, PD_dd, f10, f11)
    
    # Initial values: use column vector
    x0 <- matrix(0, nrow = 3^n, ncol = 1)
    error_dxk <- 0.03
    error_fkk <- 0.03
    num <- 10
    
    for (k in 1:num) {
      Ak <- jacobi # Replace 'jacobi' with the corresponding Jacobian matrix calculation
      bk <- f
      dxk <- pre_seidel(Ak, -bk, k)  # Step length
      x0 <- x0 + dxk
      fkk <- f  # fk+1 only used for comparison
      if (norm(dxk) < error_dxk || norm(fkk) < error_fkk) {
        break
      }
    }
    
    if (k < num) {
      x_result <- x0
    } else {
      x_result <- x0
    }
  }
  
  for (i in 1:(3^n)) {
    if (x_result[i] < 0) {
      x_result[i] <- 0
    }
  }
}


if (n == 5) {
  # Calculation of marginal penetrance of allele A
  PD_AA <- 0
  PD_Aa <- 0
  PD_aa <- 0
  for (i in 1:(3^k)) {
    PD_AA <- PD_AA + X[i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
    PD_Aa <- PD_Aa + X[3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
    PD_aa <- PD_aa + X[2 * 3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2 * 3^k + i])
  }
  PD_AA <- PD_AA - PD
  PD_Aa <- PD_Aa - PD
  PD_aa <- PD_aa - PD
  
  # Calculation of marginal penetrance of allele B
  PD_BB <- 0
  PD_Bb <- 0
  PD_bb <- 0
  for (o in 0:2) {
    for (j in 1:(3^(k-1))) {
      i <- j + o * (3^k)
      PD_BB <- PD_BB + X[i, 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2 * 3^(k-1) + i])
      PD_Bb <- PD_Bb + X[i + 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2 * 3^(k-1) + i])
      PD_bb <- PD_bb + X[i + 2 * 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2 * 3^(k-1) + i])
    }
  }
  PD_BB <- PD_BB - PD
  PD_Bb <- PD_Bb - PD
  PD_bb <- PD_bb - PD
  
  # Calculation of marginal penetrance of allele C
  PD_CC <- 0
  PD_Cc <- 0
  PD_cc <- 0
  for (o in 0:(3^(n-3)-1)) {
    for (j in 1:(3^2)) {
      i <- j + (3^(n-2)) * o
      PD_CC <- PD_CC + X[i, 1] * (D[1, i] + D[1, 3^(k-2) + i] + D[1, 2 * 3^(k-2) + i])
      PD_Cc <- PD_Cc + X[i + 3^(k-2), 1] * (D[1, i] + D[1, 3^(k-2) + i] + D[1, 2 * 3^(k-2) + i])
      PD_cc <- PD_cc + X[i + 2 * 3^(k-2), 1] * (D[1, i] + D[1, 3^(k-2) + i] + D[1, 2 * 3^(k-2) + i])
    }
  }
  PD_CC <- PD_CC - PD
  PD_Cc <- PD_Cc - PD
  PD_cc <- PD_cc - PD
  
  # Calculation of marginal penetrance of allele D
  PD_DD <- 0
  PD_Dd <- 0
  PD_dd <- 0
  for (o in 0:(3^(k-1)-1)) {
    for (j in 1:3) {
      i <- j + (3^(n-3)) * o
      PD_DD <- PD_DD + X[i, 1] * (D[1, i] + D[1, 3^(k-3) + i] + D[1, 2 * 3^(k-3) + i])
      PD_Dd <- PD_Dd + X[i + 3^(k-3), 1] * (D[1, i] + D[1, 3^(k-3) + i] + D[1, 2 * 3^(k-3) + i])
      PD_dd <- PD_dd + X[i + 2 * 3^(k-3), 1] * (D[1, i] + D[1, 3^(k-3) + i] + D[1, 2 * 3^(k-3) + i])
    }
  }
  PD_DD <- PD_DD - PD
  PD_Dd <- PD_Dd - PD
  PD_dd <- PD_dd - PD
  
  # Calculation of marginal penetrance of allele E
  PD_EE <- 0
  PD_Ee <- 0
  PD_ee <- 0
  for (o in 0:(3^k-1)) {
    j <- 1
    i <- j + (3^(n-4)) * o
    PD_EE <- PD_EE + X[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    PD_Ee <- PD_Ee + X[i + 1, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    PD_ee <- PD_ee + X[i + 2, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
  }
  PD_EE <- PD_EE - PD
  PD_Ee <- PD_Ee - PD
  PD_ee <- PD_ee - PD
  
  f10 <- D %*% X - PD
  if (H2 == 0) {
    # Calculation without using genetic force, convert to coefficient matrix
    E <- equationsToMatrix(c(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, PD_DD, PD_Dd, PD_dd, PD_EE, PD_Ee, PD_ee, f10), X)
    E <- as.matrix(E)
    f <- as.matrix(E)
    x_result <- lsqminnorm(E, f)
  } else {
    # Calculation using disease incidence, Newton's method
    f11 <- 0
    for (i in 1:(3^n)) {
      f11 <- f11 + (X[i] - PD)^2 * D[i]
    }
    f11 <- f11 - H2 * PD * (1 - PD)
    f <- c(PD_AA, PD_Aa, PD_aa, PD_BB, PD_Bb, PD_bb, PD_CC, PD_Cc, PD_cc, PD_DD, PD_Dd, PD_dd, PD_EE, PD_Ee, PD_ee, f10, f11)
    x0 <- rep(0, 3^n)
    error_dxk <- 0.05
    error_fkk <- 0.05
    num <- 10
    jacobi <- jacobian(f, X)
    for (k in 1:num) {
      Ak <- jacobi %>% 
        substitute(., X = x0) %>% 
        eval()
      bk <- f %>% 
        substitute(., X = x0) %>% 
        eval()
      dxk <- pre_seidel(Ak, -bk, k)
      x0 <- x0 + dxk
      fkk <- f %>% 
        substitute(., X = x0) %>% 
        eval()
      if (norm(dxk) < error_dxk || norm(fkk) < error_fkk) {
        break
      }
    }
    if (k < num) {
      x_result <- x0
      x_result
    } else {
      x_result <- x0
      x_result
    }
  }
  
  for (i in 1:(3^n)) {
    if (x_result[i] < 0) {
      x_result[i] <- 0
    }
  }
}



result <- verification(x_result, n, D)
v <- result[[1]]
pd <- result[[2]]
h <- result[[3]]

if (v == 0) {
  x_result <- matrix(0, nrow = 3^n, ncol = 1)
  stop("MYFUN:Badcalculation", "There is no solution to the problem defined.")
  # cat("There is no solution to the problem defined.\n")
  # x_result <- matrix(0, nrow = 3^n, ncol = 1)
}


# In the translated R code, the function verification is called, 
# and the returned values are assigned to variables v, pd, and h. 
# The subsequent code checks if v is equal to 0. If true, it assigns 
# a matrix of zeros to x_result and throws an error with a specific message using stop(). 
# The commented lines provide alternative ways to display a message and assign zeros to x_result if desired.


pre_seidel <- function(A, b, n) {
  # Citation: https://github.com/GaoBoYu599/Num_Func/tree/master/Nonlinear_Equations
  # Pre-processing: It's just that simple
  b <- t(A) %*% b
  A <- t(A) %*% A
  
  # Gauss-Seidel operations
  D <- diag(diag(A))
  L <- lower.tri(A, diag = FALSE)
  U <- upper.tri(A, diag = FALSE)
  B2 <- -solve(D + L) %*% U
  g2 <- solve(D + L) %*% b
  
  radius <- max(abs(eigen(B2)$values))
  x <- matrix(0, nrow = length(b), ncol = 1)
  error <- 0.0001
  count <- 0
  
  while (TRUE) {
    tmp <- B2 %*% x + g2
    if (max(abs(tmp - x)) < error) {
      break
    }
    x <- tmp
    count <- count + 1
  }
  
  return(x)
}

verification <- function(x_result, n, D) {
  PD_ <- matrix(0, n, 3)
  
  if (n == 2) {
    k <- n - 1  # order-1
    
    for (i in 1:(3^k)) {
      PD_[1, 1] <- PD_[1, 1] + x_result[i, 1] * (D[i] + D[3^k + i] + D[2*3^k + i])
      PD_[1, 2] <- PD_[1, 2] + x_result[3^k + i, 1] * (D[i] + D[3^k + i] + D[2*3^k + i])
      PD_[1, 3] <- PD_[1, 3] + x_result[2*3^k + i, 1] * (D[i] + D[3^k + i] + D[2*3^k + i])
    }
    
    for (j in 1:(3^k)) {
      i <- 3 * (j - 1) + 1
      
      PD_[2, 1] <- PD_[2, 1] + x_result[i, 1] * (D[i] + D[1 + i] + D[2 + i])
      PD_[2, 2] <- PD_[2, 2] + x_result[1 + i, 1] * (D[i] + D[1 + i] + D[2 + i])
      PD_[2, 3] <- PD_[2, 3] + x_result[2 + i, 1] * (D[i] + D[1 + i] + D[2 + i])
    }
  }
  
  
  if (n == 3) {
  # Calculation of marginal expressivity of alleles at site A
  k <- n - 1  # order-1
  PD_ <- matrix(0, n, 3)
  
  for (i in 1:(3^k)) {
    PD_[1, 1] <- PD_[1, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
    PD_[1, 2] <- PD_[1, 2] + x_result[3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
    PD_[1, 3] <- PD_[1, 3] + x_result[2*3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
  }
  
  # Calculation of marginal expressivity of alleles at site B
  for (o in 0:2) {
    for (j in 1:3^(k-1)) {
      i <- j + o * (3^k)
      
      PD_[2, 1] <- PD_[2, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
      PD_[2, 2] <- PD_[2, 2] + x_result[i + 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
      PD_[2, 3] <- PD_[2, 3] + x_result[i + 2*3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
    }
  }
  
  # Calculation of marginal expressivity of alleles at site C
  for (o in 0:2) {
    for (j in 0:k) {
      i <- 1 + o * 3^k + 3 * j
      
      PD_[3, 1] <- PD_[3, 1] + x_result[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
      PD_[3, 2] <- PD_[3, 2] + x_result[i + 1, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
      PD_[3, 3] <- PD_[3, 3] + x_result[i + 2, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    }
  }
}

# Note: In the original MATLAB code, D(1, i), D(1, 3^k + i), and D(1, 2*3^k + i) 
# are accessed as D[1, i], D[1, 3^k + i], and D[1, 2*3^k + i] in R, respectively.

if (n == 4) {
  # Calculation of marginal expressivity of alleles at site A
  k <- 3  # order-1
  PD_ <- matrix(0, n, 3)
  
  for (i in 1:(3^k)) {
    PD_[1, 1] <- PD_[1, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
    PD_[1, 2] <- PD_[1, 2] + x_result[3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
    PD_[1, 3] <- PD_[1, 3] + x_result[2*3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
  }
  
  # Calculation of marginal expressivity of alleles at site B
  for (o in 0:2) {
    for (j in 1:3^(k-1)) {
      i <- j + o * (3^k)
      
      PD_[2, 1] <- PD_[2, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
      PD_[2, 2] <- PD_[2, 2] + x_result[i + 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
      PD_[2, 3] <- PD_[2, 3] + x_result[i + 2*3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
    }
  }
  
  # Calculation of marginal expressivity of alleles at site C
  for (o in 0:(3^(k-1)-1)) {
    for (j in 1:3) {
      i <- j + (3^(n-2)) * o
      
      PD_[3, 1] <- PD_[3, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^1 + i] + D[1, 2*3^1 + i])
      PD_[3, 2] <- PD_[3, 2] + x_result[i + 3*1, 1] * (D[1, i] + D[1, 3^1 + i] + D[1, 2*3^1 + i])
      PD_[3, 3] <- PD_[3, 3] + x_result[i + 3*2, 1] * (D[1, i] + D[1, 3^1 + i] + D[1, 2*3^1 + i])
    }
  }
  
  # Calculation of marginal expressivity of alleles at site D
  for (o in 0:2) {
    for (j in 0:k) {
      i <- 1 + o*3^k + 3*j
      
      PD_[4, 1] <- PD_[4, 1] + x_result[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
      PD_[4, 2] <- PD_[4, 2] + x_result[i + 1, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
      PD_[4, 3] <- PD_[4, 3] + x_result[i + 2, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    }
  }
}


# Note: In the original MATLAB code, D(1, i), D(1, 3^k + i), and D(1, 2*3^k + i) 
# are accessed as D[1, i], D[1, 3^k + i], and D[1, 2*3^k + i] in R, respectively.

if (n == 5) {
  # Calculation of marginal expressivity of alleles at site A
  k <- 4  # order-1
  PD_ <- matrix(0, n, 3)
  
  for (i in 1:(3^k)) {
    PD_[1, 1] <- PD_[1, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
    PD_[1, 2] <- PD_[1, 2] + x_result[3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
    PD_[1, 3] <- PD_[1, 3] + x_result[2*3^k + i, 1] * (D[1, i] + D[1, 3^k + i] + D[1, 2*3^k + i])
  }
  
  # Calculation of marginal expressivity of alleles at site B
  for (o in 0:2) {
    for (j in 1:3^(k-1)) {
      i <- j + o * (3^k)
      
      PD_[2, 1] <- PD_[2, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
      PD_[2, 2] <- PD_[2, 2] + x_result[i + 3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
      PD_[2, 3] <- PD_[2, 3] + x_result[i + 2*3^(k-1), 1] * (D[1, i] + D[1, 3^(k-1) + i] + D[1, 2*3^(k-1) + i])
    }
  }
  
  # Calculation of marginal expressivity of alleles at site C
  for (o in 0:(3^(n-3)-1)) {
    for (j in 1:(3^2)) {
      i <- j + (3^(n-2)) * o
      
      PD_[3, 1] <- PD_[3, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^(k-2) + i] + D[1, 2*3^(k-2) + i])
      PD_[3, 2] <- PD_[3, 2] + x_result[i + 3^(k-2), 1] * (D[1, i] + D[1, 3^(k-2) + i] + D[1, 2*3^(k-2) + i])
      PD_[3, 3] <- PD_[3, 3] + x_result[i + 2*3^(k-2), 1] * (D[1, i] + D[1, 3^(k-2) + i] + D[1, 2*3^(k-2) + i])
    }
  }
  
  # Calculation of marginal expressivity of alleles at site D
  for (o in 0:(3^(k-1)-1)) {
    for (j in 1:3) {
      i <- j + (3^(n-3)) * o
      
      PD_[4, 1] <- PD_[4, 1] + x_result[i, 1] * (D[1, i] + D[1, 3^(k-3) + i] + D[1, 2*3^(k-3) + i])
      PD_[4, 2] <- PD_[4, 2] + x_result[i + 3^(k-3), 1] * (D[1, i] + D[1, 3^(k-3) + i] + D[1, 2*3^(k-3) + i])
      PD_[4, 3] <- PD_[4, 3] + x_result[i + 2*3^(k-3), 1] * (D[1, i] + D[1, 3^(k-3) + i] + D[1, 2*3^(k-3) + i])
    }
  }
  
  # Calculation of marginal expressivity of alleles at site E
  for (o in 0:(3^k-1)) {
    j <- 1
    i <- j + (3^(n-4)) * o
    
    PD_[5, 1] <- PD_[5, 1] + x_result[i, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    PD_[5, 2] <- PD_[5, 2] + x_result[i + 1, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
    PD_[5, 3] <- PD_[5, 3] + x_result[i + 2, 1] * (D[1, i] + D[1, 1 + i] + D[1, 2 + i])
  }
}

# Note: In the original MATLAB code, D(i), x_result(i), pd, PD_, and h are 
# accessed in the same way in R. 
# The * operator for matrix multiplication in MATLAB is replaced with %*% in R.

pd <- D %*% x_result
h <- 0

for (i in 1:(3^n)) {
  h <- h + D[i] * ((x_result[i] - pd)^2)
}

h <- h / (pd * (1 - pd))

if (max(abs(pd - min(PD_))) < 0.05 & max(abs(pd - max(PD_))) < 0.05) {
  v <- 1
} else {
  v <- 0
}


## not sure where the following was coming from ???
# result <- verification(x_result, n, D)
# print(result)


}