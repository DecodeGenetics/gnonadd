#' variance effect conditioned on continuous variables
#'
#' @description
#' We estimate the variance effect of a variant conditioned on one or more continous variable
#' We apply a likelyhood ratio test with one degree of freedom
#'
#' H0: All covariates have a variance effect, but not the variant
#' H1: The variant has a variance effect, and the covariates as well
#'
#' @param qt A numeric vector.
#' @param g An integer vector.
#' @param x A numeric matrix, each column represents a covariate.
#' @param iter_num An integer. Represents the number of iterations performed in the Gauss-Newton algorithm
#' @param eps_param A number. The Gauss-Newton algorithm terminates if the incriment change of all variance estimates is smaller than this number.
#'
#' @returns
#' A list with the values:
#' * alpha, the estimated variance effect, conditioned on the covariates
#' * pval, the p-value corresponding to alpha
#' @examples
#' n_val <- 50000
#' x <- matrix(0,nrow = n_val, ncol = 4)
#' for(i in 1:4) {
#' x[, i] <- rnorm(n_val)
#' }
#' g_vec <- rbinom(n_val,2,0.3)
#' var_vec <- exp(0.2 * x[, 1] - 0.3 * x[, 4] + 0.3 * g_vec)
#' qt_vec <- rnorm(n_val, 0, sqrt(var_vec))
#' res <- alpha.continuous.cond(qt_vec, g_vec, x)
#' @export alpha.continuous.cond
alpha.continuous.cond <- function(qt, g, x, iter_num = 50, eps_param = 1e-10) {
  n_subjects <- length(qt)
  n_covars <- ncol(x)
  g_factor <- as.factor(g)

  #We fit a linear model to the data and obtain residuals. We then performe the variance test on the residuals.
  #This is the only place where we treat genotype as a discrete variable.
  formula_string <- 'qt ~ g_factor'
  for(i in 1:n_covars){
    formula_string <- paste(formula_string, ' + x[, ', as.character(i), ']', sep="")
  }
  formula_string <- stats::as.formula(formula_string)
  y <- stats::lm(formula_string)$residuals

  #Null model.
  a <- alpha.multi.est(qt, x, iter_num, eps_param)
  v <- exp(- (x %*% a))
  SigmaSq_null <- sum(y^2 * v) / n_subjects
  l_null <- - n_subjects * log(SigmaSq_null)
  for(k in 1:n_covars){
    l_null <- l_null - sum(a[k] * x[,k])
  }

  #Alt model
  x_alt <- cbind(g, x)
  a <- alpha.multi.est(qt, x_alt, iter_num, eps_param)
  v <- exp(- (x_alt %*% a))
  SigmaSq_alt <- sum(y^2 * v) / n_subjects
  l_alt <- - n_subjects * log(SigmaSq_alt)
  for(k in 1:(n_covars + 1)){
    l_alt <- l_alt - sum(a[k] * x_alt[,k])
  }

  #Computation of significance
  X2 <- l_alt - l_null  #Note that I removed a factor of 1/2 when calculating l_null and l_alt. Hence there is no factor of 2 here.
  p <- stats::pchisq(X2, 1, lower.tail = FALSE)

  return(list(alpha = a[1], pval = p, X2 = X2))
}


