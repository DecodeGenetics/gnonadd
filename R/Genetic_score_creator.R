#' Creates poligenic risk scores
#'
#' @description
#' This function creates genetic risk scores, with the option of including interactions or dominance effects.
#' The score is automatically shifted to have mean 0.
#'
#'
#' @param g A matrix, where each colomn represents a variant and each line represents a subject
#' @param betas A numeric vector, representing the (additive) effects of the variants
#' @param dominance_effects A numeric vector, representing dominance effects of the variants
#' @param interaction_effects A matrix with three columns. First two columns are integers that correspond to the variants that are interacting. The third column is the effect size.
#' @param log_scale A Boolian variabe. If true all analysis is done on log-transformed effect values and the resulting score is transformed back to an exponential scale in the end.
#' @returns
#' A numeric vector with a poligenic risk score for each subject
#'
#' @examples
#' g_vec <- matrix(0, nrow = 100000, ncol = 5)
#' freqs <- runif(ncol(g_vec), min = 0, max = 1)
#' for(i in 1:ncol(g_vec)){
#'  g_vec[,i] <- rbinom(100000, 2, freqs[i])
#' }
#' beta_vec <- runif(5, min = -0.5, max = 0.5)
#' dom_vec <- runif(5, min = -0.5, max = 0.5)
#' int_vec <- matrix(0,nrow = 2, ncol = 3)
#' int_vec[, 1] <- c(1, 3)
#' int_vec[, 2] <- c(2, 5)
#' int_vec[, 3] <- runif(2, min = -0.5, max = 0.5)
#' res <- PRS_creator(g_vec, beta_vec, dominance_effects = dom_vec, interaction_effects = int_vec)
#' @export
PRS_creator <- function(g, betas, dominance_effects = rep(0,length(betas)), interaction_effects = matrix(0, nrow = 0, ncol = 0), log_scale = FALSE){
  if(log_scale == TRUE ){
    betas <- log(betas)
    dominance_effects <- log(dominance_effects)
    interaction_effects <- log(interaction_effects)
  }
  g_dom <- 1 * (round(g)==2)
  qt <- g %*% betas + g_dom %*% dominance_effects
  for(i in 1:nrow(interaction_effects)){
    g1 <- g[, interaction_effects[i, 1]]
    g2 <- g[, interaction_effects[i, 2]]
    effect <- interaction_effects[i, 3]
    qt <- qt + g1 * g2 * effect
  }
  qt <- qt - mean(qt)
  if(log_scale == TRUE){
    qt <- exp(qt)
  }
  return(qt)
}
