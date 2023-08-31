#' Calibration for the correlation test
#'
#' @description
#' This function finds appropriate scaling_parameters for the kappa_calc function by means of bootstrapping.
#'
#' @param qt1 A numeric vector.
#' @param qt2 A numeric vector.
#' @param reps The number of repeates we want to perform for each sample size
#' @returns
#' A list with the values:
#'
#' * bias_scale, a value to determine the bias for the correlation estimators
#' * weight_scale, a value to determine how much weight we assign to each correlation estimator
#' * safe_weight_scale, same as weight scale, but larger. Using this weight decreases the chance of type 1 error, with the cost of statistical power.
#'
#' @examples
#' n_val <- 10000
#' Q <- MASS::mvrnorm(n = n_val, mu = c(0,0), Sigma = matrix(c(1,0.3,0.3,1),
#'                    nrow = 2, ncol = 2))
#' qt1_val <- Q[,1]
#' qt2_val <- Q[,2]
#' res <- corr.calibration(qt1_val, qt2_val, 10)
#' @export
corr.calibration <- function(qt1, qt2, reps=10000) {

  # The 'true' fisher translated correlation value is obtained from the entire data
  n <- length(qt1)
  r_true <- stats::cor(qt1, qt2)
  z_true <- log((1 + r_true) / (1 - r_true)) / 2

  # k is a vector of sample sizes we are going to extract from the data, s is there to make sure we do not exceed the number of data-points we have
  s <- log(n / 10, 2)
  k <- 10 * 2^(0:s)
  m <- rep(0, length(k)) # A vector containing the mean of the z-values obtained from each sample size
  v <- rep(0, length(k)) # A vector containing the variance of the z-values obtained from each sample size
  for(i in 1:length(k)) {
    z_vector <- rep(0, reps) #A vector containing all the z-values obtained with the current sample size
    for(j in 1:reps) {
      s <- sample(n, k[i], replace=T)
      r_sample <- stats::cor(qt1[s], qt2[s])
      z_vector[j] <- log((1 + r_sample) / (1 - r_sample)) / 2
    }
    m[i] <- mean(z_vector) - z_true
    v[i] <- stats::var(z_vector)
  }
  weight_scale <- exp(mean(log(k - 3) + log(v)))
  safe_weight_scale <- exp(max(log(k - 3) + log(v)))
  bias_scale <- summary(stats::lm(m ~ I(1 / (k - 3)) + 0))$coeff[1, 1]
  return(list(bias_scale=bias_scale, weight_scale=weight_scale, safe_weight_scale=safe_weight_scale))
}
