% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interaction_model.R
\name{interaction.calc}
\alias{interaction.calc}
\title{Variant-Variant interaction effects}
\usage{
interaction.calc(
  qt,
  g1,
  g2,
  round_imputed = FALSE,
  dominance_terms = FALSE,
  covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))
)
}
\arguments{
\item{qt}{A numeric vector}

\item{g1}{A vector with (possibly imputed) genotype values. All entries should be larger than or equal to 0 and smaller than or equal to 2.}

\item{g2}{A vector with (possibly imputed) genotype values. All entries should be larger than or equal to 0 and smaller than or equal to 2.}

\item{round_imputed}{A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.}

\item{dominance_terms}{A boolian variable determining whether dominance terms for the variants should be included as covariates in the analysis}

\item{covariates}{A dataframe containing any other covariates that should be used; one column per covariate}
}
\value{
A list with the interaction effect and corresponding standard error, t statistic and p-value
}
\description{
This function estimates the interaction effect of a pair of genetic variant on a quantitatvie trait
We apply a simple linear regression model to estimate interaction effects.
We include a linear term for each variant seperately, coded as 0,1 and 2 for non-carriers, heterozygotes and homozygous carriers of the effect allele.
We also include an interaction term, coded as the product of the two genotype values.
Effect size and significance is based on the interaction term.
}
\examples{
g1_vec <- rbinom(100000, 2, 0.9)
g2_vec <- rbinom(100000, 2, 0.1)
qt_vec <- rnorm(100000) + 0.1 * g1_vec + 0.2 * g2_vec +0.4 * g1_vec * g2_vec
res <- interaction.calc(qt_vec, g1_vec, g2_vec)
}
