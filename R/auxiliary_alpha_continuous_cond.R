#' An auxilliary function for the alpha.continuous.cond function
#'
#' @description
#' Returns two objects. The gradiant and the hessian of the log-likelyhood function corresponding to the alpha.continuous.cond analysis.
#'
#' @param a A numeric vector consisting of variance effect estimates of each variable being tested.
#' @param y A numeric vector, representing the trait being considered
#' @param x A matrix, each column represents a covariate
#' @noRd
Aux_cont_likely <- function(a, y, x) {
  n_subjects <- length(y)
  n_covar <- ncol(x)
  
  #The following vectors/numbers appear several times in the calulations. It is most efficient to calculate them just once
  v <- exp(- (x %*% a)) #A vector of length n_subjects
  colsum_x <- rep(0,n_covar) # A vector of length n_covar. Element k is the sum of column k in x.
  for(k in 1:n_covar){
    colsum_x[k] <- sum(x[,k])
  }
  rss <- sum(y^2 * v) #This number corresponds in some sense to the rss in traditional regression analysis. Hence the name.
  rss_diff <- rep(0, n_covar) #This is the derivative of rss with respect to each a_k. Thus a vector with length n_covar.
  for(k in 1:n_covar){
    rss_diff[k] <- sum(x[,k] * y^2 * v)
  }
  
  rss_hess <- matrix(0, nrow = n_covar, ncol = n_covar) #This is the Hessian of rss with respect to the a_k variables
  for(k in 1:n_covar){
    for(r in 1:n_covar){
      rss_hess[k, r] <- sum(x[, k] * x[, r] * y^2 * v)
    }
  }
  
  #This is not really the gradiant, but rather a function with the same zeros
  Grad <- rep(0, n_covar)
  for(k in 1:n_covar) {
    Grad[k] <- n_subjects * rss_diff[k] - colsum_x[k] * rss
  }
  
  
  #This is not really the Hessian, but rather the Jacobian of Grad
  Hess <- matrix(0, nrow = n_covar, ncol = n_covar)
  for(k in 1:n_covar){
    for(r in 1:n_covar){
      Hess[k, r] <- -n_subjects * rss_hess[k, r] + colsum_x[k] * rss_diff[r]
    }
  }
  return(list(gradiant = Grad, hessian = Hess))
}
