#' Variant-Variant interaction effects
#'
#' @description
#' This function estimates the interaction effect of a pair of genetic variant on a quantitatvie trait
#' We apply a simple linear regression model to estimate interaction effects.
#' We include a linear term for each variant seperately, coded as 0,1 and 2 for non-carriers, heterozygotes and homozygous carriers of the effect allele.
#' We also include an interaction term, coded as the product of the two genotype values.
#' Effect size and significance is based on the interaction term.
#'
#' @param qt A numeric vector
#' @param g1 A vector with (possibly imputed) genotype values. All entries should be larger than or equal to 0 and smaller than or equal to 2.
#' @param g2 A vector with (possibly imputed) genotype values. All entries should be larger than or equal to 0 and smaller than or equal to 2.
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.
#' @param dominance_terms A boolian variable determining whether dominance terms for the variants should be included as covariates in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate
#'
#' @returns
#' A list with the interaction effect and corresponding standard error, t statistic and p-value
#' @examples
#' g1_vec <- rbinom(100000, 2, 0.9)
#' g2_vec <- rbinom(100000, 2, 0.1)
#' qt_vec <- rnorm(100000) + 0.1 * g1_vec + 0.2 * g2_vec +0.4 * g1_vec * g2_vec
#' res <- interaction.calc(qt_vec, g1_vec, g2_vec)
#' @export
interaction.calc <- function(qt, g1, g2, round_imputed = FALSE, dominance_terms = FALSE,
                             covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))){
  NA_vector <- is.na(g1) | is.na(g2) | is.na(qt)
  qt <- qt[!NA_vector]
  g1 <- g1[!NA_vector]
  g2 <- g2[!NA_vector]
  covariates <- covariates[!NA_vector,]
  flip <- 1   #We flip to the minor allele of each variant. This variable is used to flip the interaction effect in the end so it has the right sign.
  if(mean(g1) > 1){
    g1 <- 2 - g1
    flip <- - flip
  }
  if(mean(g2) > 1){
    flip <- - flip
    g2 <- 2 - g2
  }
  if(round_imputed == TRUE){
    g1 <- round(g1)
    g2 <- round(g2)
  }
  int <- g1 * g2
  if(stats::sd(int) == 0) {
    warning("Interaction undefined. All interaction values are the same.")
    gamma <- NA
    se <- NA
    t <- NA
    p <- NA
  } else{

    #We define a dataframe containing all variables that should be considered
    Int_data <- as.data.frame(cbind(qt, int))
    Int_data <- cbind(Int_data, g1)
    Int_data <- cbind(Int_data, g2)
    if(dominance_terms == TRUE){
      Int_data$g1_dom <- as.numeric(round(g1) == 2)
      Int_data$g2_dom <- as.numeric(round(g2) == 2)
    }
    if(nrow(covariates) > 0){
      Int_data <-cbind(Int_data, covariates)
    }

    #We use linear regression to estimate the interaction effect
    l_interaction <- stats::lm(qt ~ ., data = Int_data)
    param <- "int"
    if(param %in% rownames(stats::coef(summary(l_interaction)))){
      gamma <- summary(l_interaction)$coeff[param, 1]
      se <- summary(l_interaction)$coeff[param, 2]
      t <- summary(l_interaction)$coeff[param, 3]
      p <- summary(l_interaction)$coeff[param, 4]
    }else{
      warning("Singular model matrix")
      gamma <- NA
      se <- NA
      t <- NA
      p <- NA
    }
  }
  gamma <- flip * gamma
  return(list(interaction_effect = gamma, standard_error = se, t = t, pval = p))
}
