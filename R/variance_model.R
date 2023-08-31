#' Genetic variance effects
#'
#' @description
#' This function estimates the variance effect of a genetic variant on a quantitative trait.
#' The genotype is coded as 0 (non-carrier), 1 (single copy of effect allele) and 2 (two copies of effect allele).
#' The variance effect (alpha) is modelled to be multiplicitive.
#' We use a likelihood ratio test with 1 degree of freedom,
#' * H0: y~N(mu_g,sigma^2)
#' * H1: y~N(mu_g,sigma^2\*alpha^g)
#' Under the alternative model, the variance of the trait changes multiplicatively
#' with genotype:
#'     Var(y|g=j)=alpha^j*sigma^2, j in 0,1,2
#'
#' @param qt A numeric vector.
#' @param g An integer vector.
#' @returns
#' A list with the values:
#'
#' * alpha, the estimated variance effect
#' * sigma2_alt, The variance estimated for non-carriers
#' * pval, the p-value of the likelihood ratio test
#' * chi2, the chi squared statistics (with one degree of freedom) corrisponding to the p-value
#' @examples
#' n_val <- 50000L
#' geno_vec <- sample(c(0, 1, 2), size = n_val, replace = TRUE)
#' qt_vec <- rnorm(n_val) * (1.1^geno_vec)
#' res <- alpha.calc(qt_vec, geno_vec)
#' @export
alpha.calc <- function(qt, g) {
  g <- round(g)  # We round imputed values to the nearest integer
  n <- length(qt)
  f <- mean(g) / 2.0

  flip <- 0
  if(f > 0.5) {
    # If there are no non-carriers of the effect allele, the programs attempts division by 0.
    # This means that rss0 will be zero, or we would divide by a small number.
    # Therefore we always do the calculations in terms of the minor allele. We flip the alleles back in the end
    g <- 2 - g
    flip <- 1
    f <- mean(g) / 2
  }

  m0 <- mean(qt[g == 0])
  m1 <- mean(qt[g == 1])
  m2 <- mean(qt[g == 2])
  rss0 <- sum((qt[g == 0] - m0)^2)
  rss1 <- sum((qt[g == 1] - m1)^2)
  rss2 <- sum((qt[g == 2] - m2)^2)

  # Estimators for null-model and log-likelihood function
  sigma2_null <- (rss0 + rss1 + rss2) / n
  l0 <- -n * log(sigma2_null) / 2.0

  # Estimators for alt-model and log-likelihood function
  a <- 2 * f * rss0
  b <- (2 * f - 1) * rss1
  c <- (2 * f - 2) * rss2

  alpha <- (- b + sqrt(b^2 - 4 * a * c)) / (2 * a)
  sigma2_alt <- (rss2 / (alpha^2) + rss1 / alpha + rss0) / n

  l1 <- - n * log(sigma2_alt) / 2 - log(alpha) * n * f

  # Calculation of significance
  X2 <- 2 * (l1 - l0)
  pval <- stats::pchisq(X2, 1, lower.tail=F)

  if(flip >0.5){
    # We flip the parameters back, so they are in accordence with the input effect allele
    sigma2_alt <- sigma2_alt * (alpha^2)
    alpha <- 1 / alpha
  }
  return(list(alpha = alpha, sigma2_alt = sigma2_alt, pval = pval, chi2 = X2))
}
