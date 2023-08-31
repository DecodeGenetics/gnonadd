#' Variant-Variant interaction effects on a case control variable
#'
#' @description
#' This function estimates the interaction effect of a pair of genetic variant on a case-control variable
#' We apply a logistic regression model to estimate interaction effects.
#' We include a linear term for each variant seperately, coded as 0,1 and 2 for non-carriers, heterozygotes and homozygous carriers of the effect allele.
#' We also include an interaction term, coded as the product of the two genotype values.
#' Effect size and significance is based on the interaction term.
#'
#' @param cc A numeric vector
#' @param g1 A vector with (possibly imputed) genotype values. All entries should be larger than 0 and smaller than 2.
#' @param g2 A vector with (possibly imputed) genotype values. All entries should be larger than 0 and smaller than 2.
#' @param yob A numerical vector containing year of birth. If some are unknown they should be marked as -1
#' @param sex A numerical vector containing sex, coded 0 for males, 1 for females and -1 for unknown
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.
#' @param dominance_terms A boolian variable determining whether dominance terms for the variants should be included as covariates in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate
#'
#' @returns
#' A list with the interaction effect (on log-scale) and corresponding standard error, z statistic and p-value
#' @examples
#' g1_vec <- rbinom(100000, 2, 0.9)
#' g2_vec <- rbinom(100000, 2, 0.1)
#' cc_vec <- rbinom(100000,1,0.1 * (1.05^g1_vec) *
#'           (1.05^g2_vec) * (1.3 ^ (g1_vec * g2_vec)))
#' res <- interaction_CC.calc(cc_vec, g1_vec, g2_vec)
#' @export
interaction_CC.calc <- function(cc, g1, g2, yob = rep(-1,length(cc)),
                                sex = rep(-1,length(cc)), round_imputed = FALSE, dominance_terms = FALSE, covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))){
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
  if(stats::sd(int) == 0){
    warning("Interaction undefined. All interaction values are the same.")
    gamma <- NA
    se <- NA
    z <- NA
    p <- NA
  } else{
    no_date <- yob < 0
    yob[!no_date] <- yob[!no_date] - mean(yob[!no_date])
    yob[no_date] <- 0
    sex <- as.factor(sex)

    #We define a dataframe containing all variables that should be considered
    Int_data <- as.data.frame(cbind(cc, int))
    Int_data <- cbind(Int_data, g1)
    Int_data <- cbind(Int_data, g2)
    if(dominance_terms == TRUE){
      Int_data$g1_dom <- as.numeric(round(Int_data$g1) == 2)
      Int_data$g2_dom <- as.numeric(round(Int_data$g2) == 2)
    }
    if(stats::sd(yob) > 0) {
      Int_data <- cbind(Int_data, yob)
    }
    if(length(unique(no_date)) > 1) {
      Int_data <- cbind(Int_data, no_date)
    }
    if(length(unique(sex)) > 1) {
      Int_data <- cbind(Int_data, sex)
    }
    if(nrow(covariates) > 0) {
      Int_data <- cbind(Int_data, covariates)
    }

    #We use logistic regression to estimate the interaction effect
    l_interaction <- stats::glm(cc ~ ., data = Int_data, family = 'binomial')
    param <- "int"
    if(param %in% rownames(stats::coef(summary(l_interaction)))){
      gamma <- summary(l_interaction)$coeff[param, 1]
      se <- summary(l_interaction)$coeff[param, 2]
      z <- summary(l_interaction)$coeff[param, 3]
      p <- summary(l_interaction)$coeff[param, 4]
    }else{
      warning("Singular model matrix")
      gamma <- NA
      se <- NA
      z <- NA
      p <- NA
    }


    gamma <- flip * gamma
  }
  return(list(interaction_effect = gamma, standard_error = se, z = z, pval = p))
}
