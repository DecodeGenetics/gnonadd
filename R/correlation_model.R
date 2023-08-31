#' Genetic correlation effects
#'
#' @description
#' This function estimates the correlation effect of a genetic variant on a pair of quantitative traits.
#' The effect (kappa) is measured on Fisher transformed correlation values.
#' The genotype is coded as 0 (non-carrier), 1 (single copy of effect allele) and 2 (two copies of effect allele).
#' We use a likelihood ratio test with 1 degree of freedom:
#' * H0: z_0 = z_1 = z_2,
#' * H1: z_j = z_null +kappa*j,
#' where z_j is the Fisher-transformed correlation between the traits amongst subjects of genotype j.
#'
#' If the traits follow a joint normal distribution, then (thoretically) the Fisher-transformed sample correlation (z) is approximately normally distributed with mean 0 and
#' standard error 1/sqrt(n-3).
#' Otherwise, we have to correct for possible bias in the estimators and scale the weights we assign to them appropriately
#'
#' Note, that even though each trait follows a normal distribution individually, that does not necessarily imply that the pair of them follow a joint normal distirbution.
#'
#' @param qt1 A numeric vector.
#' @param qt2 A numeric vector.
#' @param g An integer vector.
#' @param weight_scale Used to appropriately scale the weight assigned to the correlation estimators
#' @param bias_scale Used to appropriately scale the bias of the correlation estimators
#' @returns
#' A list with the values:
#'
#' * kappa, the estimated correlation effect
#' * pval, the p-value of the likelihood ratio test
#' @examples
#' Sigma0 <- matrix(c(1,0,0,1),nrow=2,ncol=2)
#' Sigma1 <- matrix(c(1,0.3,0.3,1),nrow=2,ncol=2)
#' Sigma2 <- matrix(c(1,0.6,0.6,1),nrow=2,ncol=2)
#' geno_vec <- c(rep(0,10000),rep(1,1000),rep(2,100))
#' Q0 <- MASS::mvrnorm(n = 10000, mu = c(0,0), Sigma = Sigma0)
#' Q1 <- MASS::mvrnorm(n = 1000, mu = c(0,0), Sigma = Sigma1)
#' Q2 <- MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = Sigma2)
#' Q <- rbind(Q0,Q1)
#' Q <- rbind(Q,Q2)
#' qt1_val <- Q[,1]
#' qt2_val <- Q[,2]
#' res <- kappa_calc(qt1_val, qt2_val, geno_vec)
#' @export kappa_calc
kappa_calc <- function(qt1, qt2, g, weight_scale=1, bias_scale=0) {
  g <- round(g)
  n0 <- sum(g==0)
  n1 <- sum(g==1)
  n2 <- sum(g==2)

  # We calculate the correlation (r) for each genotype.
  # We use a Fisher-transformation, corrected for possible bias, to obtain a parameter z that has an approximate normal distribution.
  # We work under the assumption that the variance of the random parameter z is approximately porportional to 1/(n-3).
  # Therefore the corresponding weight (w) is porportional to n-3.
  # If qt1 and qt2 are jointly normal then (theoretically) there is no bias and the weight simply equals n-3.
  if(n0 > 3) {
    r0 <- stats::cor(qt1[g==0], qt2[g==0])
    z0 <- log((1 + r0) / (1 - r0)) / 2 - bias_scale / (n0 - 3)
    w0 <- (n0 - 3) / weight_scale
  }else {
    r0 <- 0
    z0 <- 0
    w0 <- 0
  }
  if(n1 > 3) {
    r1 <- stats::cor(qt1[g==1], qt2[g==1])
    z1 <- log((1 + r1) / (1 - r1)) / 2 - bias_scale / (n1 - 3)
    w1 <- (n1 - 3) / weight_scale
  }else {
    r1 <- 0
    z1 <- 0
    w1 <- 0
  }
  if(n2 > 3) {
    r2 <- stats::cor(qt1[g==2], qt2[g==2])
    z2 <- log((1 + r2) / (1 - r2)) / 2 - bias_scale / (n2 - 3)
    w2 <- (n2 - 3) / weight_scale
  }else {
    r2 <- 0
    z2 <- 0
    w2 <- 0
  }

  #Estimators for null model and log-likelihood function
  z_null <- (w0 * z0 + w1 * z1 + w2 * z2) / (w0 + w1 + w2)
  l_null <- (w0 * (z0 - z_null)^2 + w1 * (z1 - z_null)^2 + w2 * (z2 - z_null)^2) / 2

  #Estimators for alt model and log-likelyhood function
  D <- (w1 + 4 * w2) * (w0 + w1 + w2) - (w1 + 2 * w2)^2
  z_alt <- ((w1 + 4 * w2) * (w0 * z0 + w1 * z1 + w2 * z2) - (w1 + 2 * w2) * (z1 * w1 + 2 * z2 * w2)) / D
  kappa <- ((z1 * w1 + 2 * z2 * w2) * (w0 + w1 + w2) - (w0 * z0 + w1 * z1 + w2 * z2) * (w1 + 2 * w2)) / D
  l_alt <- (w0 * (z0 - z_alt)^2 + w1 * (z1 - z_alt - kappa)^2 + w2 * (z2 - z_alt - 2 * kappa)^2) / 2

  #Calculation of significance
  p <- stats::pchisq(2 * (l_null - l_alt), 1,lower.tail = FALSE)
  return(list(kappa = kappa, pval = p))
}
