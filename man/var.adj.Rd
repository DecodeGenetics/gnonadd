% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Variance_adjustments.R
\name{var.adj}
\alias{var.adj}
\title{Mean and variance effect adjustments.}
\usage{
var.adj(qt, x, iter_num = 50, eps_param = 1e-10)
}
\arguments{
\item{qt}{A numeric vector.}

\item{x}{A numeric matrix, each column represents a covariate.}

\item{iter_num}{An integer. Represents the number of iterations performed in the Gauss-Newton algorithm}

\item{eps_param}{A number. The Gauss-Newton algorithm terminates if the incriment change of all variance estimates is smaller than this number.}
}
\value{
A vector, representing the adjusted trait.
}
\description{
Given is a set of (continuous) variables and a qt trait. First, this function adjusts the trait for the mean effects of the variables with a linear model.
Next, the variance effect of the variables are estimated and the trait is adjusted further by scaling it in accordance with the results.
}
\examples{
n_val <- 50000
x <- matrix(0,nrow = n_val, ncol = 4)
for(i in 1:4) {
x[, i] <- rnorm(n_val)
}
var_vec <- exp(0.2 * x[, 1] - 0.3 * x[, 4])
qt_vec <- rnorm(n_val, 0, sqrt(var_vec))
res <- var.adj(qt_vec, x)
}
