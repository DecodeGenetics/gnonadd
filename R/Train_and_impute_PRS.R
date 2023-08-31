#' Trains and imputes a poligenic risk score (PRS)
#'
#' @description
#' This function trains a poligenic risk score model on a dataset, and then imputes the risk score into another dataset
#'
#'
#' @param qt_training A numeric vector. Represents the qt values of the data we train the model on.
#' @param g_training A matrix, where each colomn represents a variant and each line represents a subject in the training data
#' @param g_impute A matrix, where each column represents a variant and each line represents a subject.
#' @param dominance_effects A Boolian vector, each term determines whether a dominance term for the corresponding variant is used in the model.
#' @param interaction_effects An integer matrix with two columns. Each line represents a pair of interacting variants that should be included in the model.

#' @returns
#' Returns a list with the following objects
#' * PRS_imputed, the imputed PRS values
#' * PRS_training, the PRS values for the training data
#' * Residuals_training, the residuals from the model in the training data
#'
#'
#' @examples
#' g_train_vec <- matrix(0, nrow = 100000, ncol = 5)
#' freqs <- runif(ncol(g_train_vec), min = 0, max = 1)
#' for(i in 1:ncol(g_train_vec)){
#'  g_train_vec[,i] <- rbinom(100000, 2, freqs[i])
#' }
#' g_impute_vec <- matrix(0, nrow = 50000, ncol = 5)
#' for(i in 1:ncol(g_impute_vec)){
#'  g_impute_vec[,i] <- rbinom(50000, 2, freqs[i])
#' }
#' dom_vec <- c(TRUE, FALSE, FALSE, TRUE, FALSE)
#'
#' int_vec <- matrix(c(1, 2, 4, 5), nrow = 2 , ncol = 2)
#'
#' qt_vec <- rnorm(100000) + 0.2 * g_train_vec[, 1] + 0.3 * g_train_vec[, 1] * g_train_vec[, 4]
#'
#' res <- train_and_impute_PRS(qt_vec, g_train_vec, g_impute_vec,
#'        dominance_effects = dom_vec, interaction_effects = int_vec)
#' @export
train_and_impute_PRS <- function(qt_training, g_training, g_impute, dominance_effects = rep(FALSE, ncol(g_training)), interaction_effects = matrix(0, nrow = 0, ncol = 0)){
  Data_training <- as.data.frame(cbind(qt_training, g_training))
  colnames(Data_training) <- c('qt',paste(rep('g_',ncol(g_training)),as.character(1:ncol(g_training)),sep=""))
  Data_impute <- as.data.frame(g_impute)
  colnames(Data_impute) <- paste(rep('g_',ncol(g_training)),as.character(1:ncol(g_training)),sep="")
  dom_names <- paste(rep('dom_',ncol(g_training)),as.character(1:ncol(g_training)),sep="")
  dom_names <- dom_names[dominance_effects]
  dom_data_training <- as.data.frame(matrix(0, nrow = nrow(Data_training), ncol = sum(dominance_effects)))
  colnames(dom_data_training) <- dom_names
  dom_data_impute <- as.data.frame(matrix(0, nrow = nrow(Data_impute), ncol = sum(dominance_effects)))
  colnames(dom_data_impute) <- dom_names
  counter <- 0
  for(i in 1:length(dominance_effects)){
    if(dominance_effects[i]){
      counter <- counter + 1
      dom_data_training[, counter] <- as.numeric(round(g_training[, i]) == 2)
      dom_data_impute[, counter] <- as.numeric(round(g_impute[, i]) == 2)

    }
  }
  Data_training <- cbind(Data_training, dom_data_training)
  Data_impute <- cbind(Data_impute, dom_data_impute)
  int_data_training <- as.data.frame(matrix(0, nrow = nrow(Data_training), ncol = nrow(interaction_effects)))
  int_data_impute <- as.data.frame(matrix(0, nrow = nrow(Data_impute), ncol = nrow(interaction_effects)))
  int_names <- rep('',nrow(interaction_effects))
  for(i in 1:length(int_names)){
    int_names[i] <- paste('int', as.character(interaction_effects[i,1]), as.character(interaction_effects[i,2]),sep="_")
  }
  colnames(int_data_training) <- int_names
  colnames(int_data_impute) <- int_names
  for(i in 1:nrow(interaction_effects)){
    num1 <- interaction_effects[i, 1]
    num2 <- interaction_effects[i, 2]
    int_data_impute[, i] <- g_impute[, num1] * g_impute[, num2]
    int_data_training[, i] <- g_training[, num1] * g_training[, num2]
  }
  Data_training <- cbind(Data_training, int_data_training)
  Data_impute <- cbind(Data_impute, int_data_impute)
  l <- stats::lm(qt~., data = Data_training)
  PRS_impute <- stats::predict(l,Data_impute)
  PRS_training <- stats::predict(l,Data_training)
  r <- l$residuals
  return(list(PRS_imputed = PRS_impute, PRS_training = PRS_training, Residuals_training = r))
}
