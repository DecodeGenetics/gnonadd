#' Variance parameters
#'
#' @description
#' This function jointly estimates the variance effect of a set of (continuous) variables on a qt trait.
#' More precisely. It finds the maximum likelyhood estimators.
#'
#' @param qt A numeric vector.
#' @param x A numeric matrix, each column represents a covariate.
#' @param iter_num An integer. Represents the number of iterations performed in the Gauss-Newton algorithm
#' @param eps_param A number. The Gauss-Newton algorithm terminates if the incriment change of all variance estimates is smaller than this number.
#' @param initial_guess A vector of length ncol(x). Represents the initial guess of parameters for the Gauss-Newton algorithm.
#'
#' @returns
#' A vector with a variance estimate for each variable.

#' @examples
#' n_val <- 50000
#' x <- matrix(0,nrow = n_val, ncol = 4)
#' for(i in 1:4) {
#' x[, i] <- rnorm(n_val)
#' }
#' var_vec <- exp(0.2 * x[, 1] - 0.3 * x[, 4])
#' qt_vec <- rnorm(n_val, 0, sqrt(var_vec))
#' res <- alpha.multi.est(qt_vec, x)
#' @export alpha.multi.est
alpha.multi.est <- function(qt, x, iter_num = 50, eps_param = 1e-10, initial_guess = rep(0,ncol(x))) {
  n_subjects <- length(qt)
  n_covars <- ncol(x)
  a <- initial_guess  #This is the original guess of variance parameters for the covariates.
  for(i in 1:iter_num){
    A <- Aux_cont_likely(a, qt, x)
    increment <- solve(A$hessian) %*% A$gradiant
    a <- a - increment
    if(max(abs(increment)) < eps_param){
      break
    }
  }
  return(a)
}


