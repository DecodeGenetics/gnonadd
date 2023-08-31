#' Uncertanty association
#'
#' @description
#' This function finds the association between the predicted uncertanty of some estimates of a trait to the "actual uncertanty" of the estimates.
#' Suppose we have estimates of some trait (this might be a polygenic risk score).
#' Moreover, suppose we have assigned a variance value to each estimate (this might be a variance risk score) to reflect
#' how certain we believe we are about each estimate. Given the true trait values, this function evaluates how well the assigned
#' variance values reflect reality.
#'
#' We use a likelihood ratio test with 1 degree of freedom
#' * H0: y~N(mu+b*m_score,sigma_sq),
#' * H1: y~N(mu+b*m_score,sigma_sq*(v_score)^a),
#' where y is the trait m_score is the estimate of the trait and v_score is the variance assigned to the estimate.
#' Thus H0 has three degrees of freedom (mu,b,sigma_sq), whereas H1 has four (mu,b,sigma_sq,a)
#'
#'
#' @param qt A numeric vector.
#' @param m_score A numeric vector
#' @param v_score A numeric vector with positive values
#' @param iter A number of iterations for the Gauss-Newton algorithm
#' @returns
#' A list with the values:
#'
#' * a, the association between v_score and the actual variance.
#' * pval, the p-value of the likelihood ratio test
#' @examples
#'
#' n_val <- 50000L
#' trait_vec <- rnorm(n_val,0,1)
#' var_vec <- exp(rnorm(n_val,0,0.1))
#' est_vec <- trait_vec+rnorm(n_val,0,var_vec)
#' res <- Var.assoc(trait_vec,est_vec,var_vec, iter = 20)
#' @export
Var.assoc <- function(qt, m_score, v_score, iter=50) {
  n <- length(qt)
  v_score <- log(v_score) # We transform the v_score to an additive scale since it is easier to work with.
  l <- stats::lm(qt~m_score) # We use the values from the regression as initial guess parameters in the Gauss-Newton algorithm
  sigma_sq <- stats::var(l$res) * (n - 1) /n # inital guess for the sigma^2 parameter
  a <- 0 # inital guess for the association with the v_score
  mu <- summary(l)$coeff[1, 1] # Initial guess for the constant term in the m_score association
  b <- summary(l)$coeff[2, 1] # Inital guess for the m_score association
  x <- c(sigma_sq,a,mu,b) # Initial guesses
  H <- matrix(0, 4, 4) # A matrix, that will represent the Hessian of the log-likelyhood function at the point x in each iteration of the algorithm
  for(i in 1:iter){

    # The following values represent the derivative of the log-likelyhood function evaluated at x
    d_sigma_sq <- - n * sigma_sq + sum((qt - mu - b * m_score)^2 * exp(- a * v_score))
    d_a <- -sigma_sq * (sum(v_score)) + sum((qt - mu - b * m_score)^2 * v_score * exp(-a * v_score))
    d_mu <- sum((qt - mu - b * m_score) * exp(- a * v_score))
    d_b <- sum((qt - mu - b * m_score) * m_score * exp(- a * v_score))
    d_x <- c(d_sigma_sq, d_a, d_mu, d_b)

    # H is the Hessian of the log-likelyhood function evaluated at x
    H[1, 1] <- - n
    H[1, 2] <- - sum((qt - mu - b * m_score)^2 * v_score * exp( - a * v_score))
    H[1, 3] <- - 2 * sum((qt - mu - b * m_score) * exp(- a * v_score))
    H[1, 4] <- - 2 * sum((qt - mu - b * m_score) * m_score * exp(- a * v_score))
    H[2, 1] <- - sum(v_score)
    H[2, 2] <- - sum((qt - mu - b * m_score)^2 * (v_score^2) * exp(- a * v_score))
    H[2, 3] <- - 2 * sum((qt - mu - b * m_score) * v_score * exp( - a * v_score))
    H[2, 4] <- - 2 * sum((qt - mu - b * m_score) * v_score * m_score * exp(- a * v_score))
    H[3, 1] <- 0
    H[3, 2] <- - sum((qt - mu - b * m_score) * v_score * exp(- a * v_score))
    H[3, 3] <- - sum(exp(- a * v_score))
    H[3, 4] <- - sum(m_score * exp(- a * v_score))
    H[4, 1] <- 0
    H[4, 2] <- - sum((qt - mu - b * m_score) * v_score * m_score * exp(- a * v_score))
    H[4, 3] <- - sum(m_score * exp(- a * v_score))
    H[4, 4] <- - sum((m_score^2) * exp(- a * v_score))
    x <- x - solve(H) %*% d_x
    sigma_sq <- x[1]
    a <- x[2]
    mu <- x[3]
    b <- x[4]
  }

  # Evaluated likelyhood-functions
  l_null <- - n * log((n-1) * stats::var(l$residuals) / n)
  l_alt <- - n * log(sigma_sq) - a * sum(v_score)
  xi2 <- l_alt - l_null
  p <- stats::pchisq(xi2, 1, lower.tail=F)
  return(list(association=a, pval=p))
}

