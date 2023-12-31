% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Violin_by_genotype.R
\name{Viol.by.gen}
\alias{Viol.by.gen}
\title{Violin plot by genotype}
\usage{
Viol.by.gen(qt, g, trait_name = "qt trait", title = "")
}
\arguments{
\item{qt}{A numeric vector.}

\item{g}{An integer vector.}

\item{trait_name}{A string.}

\item{title}{A string.}
}
\value{
A violin plot
}
\description{
This tool creates violin plots corresponding to each genotype.
}
\examples{
n_val <- 50000L
geno_vec <- sample(c(0, 1, 2), size = n_val, replace = TRUE)
qt_vec <- rnorm(n_val) * (1.3^geno_vec) + 1 * geno_vec
Viol.by.gen(qt_vec, geno_vec)
}
