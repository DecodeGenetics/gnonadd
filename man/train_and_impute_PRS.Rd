% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Train_and_impute_PRS.R
\name{train_and_impute_PRS}
\alias{train_and_impute_PRS}
\title{Trains and imputes a poligenic risk score (PRS)}
\usage{
train_and_impute_PRS(
  qt_training,
  g_training,
  g_impute,
  dominance_effects = rep(FALSE, ncol(g_training)),
  interaction_effects = matrix(0, nrow = 0, ncol = 0)
)
}
\arguments{
\item{qt_training}{A numeric vector. Represents the qt values of the data we train the model on.}

\item{g_training}{A matrix, where each colomn represents a variant and each line represents a subject in the training data}

\item{g_impute}{A matrix, where each column represents a variant and each line represents a subject.}

\item{dominance_effects}{A Boolian vector, each term determines whether a dominance term for the corresponding variant is used in the model.}

\item{interaction_effects}{An integer matrix with two columns. Each line represents a pair of interacting variants that should be included in the model.}
}
\value{
Returns a list with the following objects
* PRS_imputed, the imputed PRS values
* PRS_training, the PRS values for the training data
* Residuals_training, the residuals from the model in the training data
}
\description{
This function trains a poligenic risk score model on a dataset, and then imputes the risk score into another dataset
}
\examples{
g_train_vec <- matrix(0, nrow = 100000, ncol = 5)
freqs <- runif(ncol(g_train_vec), min = 0, max = 1)
for(i in 1:ncol(g_train_vec)){
 g_train_vec[,i] <- rbinom(100000, 2, freqs[i])
}
g_impute_vec <- matrix(0, nrow = 50000, ncol = 5)
for(i in 1:ncol(g_impute_vec)){
 g_impute_vec[,i] <- rbinom(50000, 2, freqs[i])
}
dom_vec <- c(TRUE, FALSE, FALSE, TRUE, FALSE)

int_vec <- matrix(c(1, 2, 4, 5), nrow = 2 , ncol = 2)

qt_vec <- rnorm(100000) + 0.2 * g_train_vec[, 1] + 0.3 * g_train_vec[, 1] * g_train_vec[, 4]

res <- train_and_impute_PRS(qt_vec, g_train_vec, g_impute_vec,
       dominance_effects = dom_vec, interaction_effects = int_vec)
}
