#' Histogram by genotype
#'
#' @description
#' This tool creates three histogram plots. One per genotype. Additionally, outliers are colored red (by default subjects that are in the top and bottom 2.5% of the distribution)
#' and blue lines are added to indicate the mean and (by default) one standard deviation in each direction.
#'
#' @param qt A numeric vector.
#' @param g An integer vector.
#' @param bins An integer.
#' @param trait_name A string.
#' @param title A string.
#' @param outlier_quantiles A vector with length 2.
#' @param sd_lines A vector with length 2.
#' @returns
#' A histgram plot
#'
#' @examples
#' n_val <- 50000L
#' geno_vec <- sample(c(0, 1, 2), size = n_val, replace = TRUE)
#' qt_vec <- rnorm(n_val) * (1.3^geno_vec) + 1 * geno_vec
#' hist_by_gen(qt_vec, geno_vec)
#' @export hist_by_gen
hist_by_gen <- function(qt, g, bins = 100, trait_name = 'qt trait', title = '', outlier_quantiles = c(0.025, 0.975), sd_lines = c(1,1)) {
  outlier <- NULL
  g <- round(g)
  D <- as.data.frame(cbind(qt,g))
  D$g_factor <- factor(D$g,levels=0:2, labels = c('Non-carriers', 'Heterozygotes', 'Homozygotes'))
  quant_vals <- as.numeric(stats::quantile(D$qt, outlier_quantiles))
  D$outlier <- as.factor(0+(D$qt<quant_vals[1] | D$qt > quant_vals[2]))
  a <- sd_lines[1]
  b <- sd_lines[2]
  vertical_lines <- data.frame(g_factor = levels(D$g_factor), mean = c(mean(D$qt[D$g == 0]), mean(D$qt[D$g == 1]), mean(D$qt[D$g == 2])),
                         lower_sd <- c(mean(D$qt[D$g == 0]) - a * stats::sd(D$qt[D$g == 0]), mean(D$qt[D$g == 1]) - a * stats::sd(D$qt[D$g == 1]), mean(D$qt[D$g == 2]) - a * stats::sd(D$qt[D$g == 2])),
                         upper_sd <- c(mean(D$qt[D$g == 0]) + b * stats::sd(D$qt[D$g == 0]), mean(D$qt[D$g == 1]) + b * stats::sd(D$qt[D$g == 1]), mean(D$qt[D$g == 2]) + b * stats::sd(D$qt[D$g == 2])))

  ggplot2::ggplot(D, ggplot2::aes(x = qt, fill = outlier)) + ggplot2::geom_histogram(color = 'black', bins = bins) + ggplot2::theme_classic() +
    ggplot2::facet_grid(g_factor ~ . , scales = "free_y") + ggplot2::xlab(trait_name) + ggplot2::ggtitle(title) +
    ggplot2::scale_fill_manual(values=c('white', 'red'))+
    ggplot2::geom_vline(data=vertical_lines,ggplot2::aes(xintercept = mean),color='blue', size = 1)+
    ggplot2::geom_vline(data=vertical_lines,ggplot2::aes(xintercept = lower_sd),color='blue',linetype='dashed', size = 1)+
    ggplot2::geom_vline(data=vertical_lines,ggplot2::aes(xintercept = upper_sd),color='blue',linetype='dashed', size = 1)+
    ggplot2::theme(legend.position = "none")
}
