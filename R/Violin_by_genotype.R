#' Violin plot by genotype
#'
#' @description
#' This tool creates violin plots corresponding to each genotype.
#'
#' @param qt A numeric vector.
#' @param g An integer vector.
#' @param trait_name A string.
#' @param title A string.
#' @returns
#' A violin plot
#'
#' @examples
#' n_val <- 50000L
#' geno_vec <- sample(c(0, 1, 2), size = n_val, replace = TRUE)
#' qt_vec <- rnorm(n_val) * (1.3^geno_vec) + 1 * geno_vec
#' Viol.by.gen(qt_vec, geno_vec)
#' @export
Viol.by.gen <- function(qt, g, trait_name = 'qt trait', title = '') {
  g_factor <- NULL
  g <- round(g)
  D <- as.data.frame(cbind(qt, g))
  D$g_factor <- factor(D$g, levels = 0:2,
                       labels = c('Non-carriers',
                                  'Heterozygotes',
                                  'Homozygotes'))
  ggplot2::ggplot(D, ggplot2::aes(x = g_factor, y = qt)) +
    ggplot2::geom_violin()+ggplot2::theme_classic()+
    ggplot2::xlab('Genotype')+ggplot2::ylab(trait_name) +
    ggplot2::ggtitle(title)+ggplot2::geom_boxplot(width = 0.1)
}
