#' Mean and variance effect adjustments.
#'
#' @description
#' Given is a set of (continuous) variables and a qt trait. First, this function adjusts the trait for the mean effects of the variables with a linear model.
#' Next, the variance effect of the variables are estimated and the trait is adjusted further by scaling it in accordance with the results.
#'
#' @param qt A numeric vector.
#' @param x A numeric matrix, each column represents a covariate.
#' @param iter_num An integer. Represents the number of iterations performed in the Gauss-Newton algorithm
#' @param eps_param A number. The Gauss-Newton algorithm terminates if the incriment change of all variance estimates is smaller than this number.
#'
#' @returns
#' A vector, representing the adjusted trait.
#' @examples
#' n_val <- 50000
#' x <- matrix(0,nrow = n_val, ncol = 4)
#' for(i in 1:4) {
#' x[, i] <- rnorm(n_val)
#' }
#' var_vec <- exp(0.2 * x[, 1] - 0.3 * x[, 4])
#' qt_vec <- rnorm(n_val, 0, sqrt(var_vec))
#' res <- var.adj(qt_vec, x)
#' @export var.adj
var.adj <- function(qt, x, iter_num = 50, eps_param = 1e-10){
  x <- as.matrix(x)
  n_subjects <- length(qt)
  n_covars <- ncol(x)
  A <- as.data.frame(cbind(qt, x))
  colnames(A)[1] <- "qt"
  qt_mean_adj <- stats::lm(qt ~ ., data = A)$residuals
  a <- alpha.multi.est(qt_mean_adj, x, iter_num = 50, eps_param = 1e-10)
  adjustment_scalar <- 1 / exp((x %*% a) / 2)
  qt_adj <- base::scale(qt_mean_adj * adjustment_scalar)
  return(qt_adj)
}
