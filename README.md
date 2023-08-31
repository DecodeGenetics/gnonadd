
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnonadd

<!-- badges: start -->
<!-- badges: end -->

The goal of the `gnonadd` package is to simplify workflows with
non-additive analysis in genetic associations.

This includes e.g.

1)  Variance effects
2)  Correlation effects
3)  Interaction effects
4)  Dominance effects

## Included Functionality

The following is a non-comprehensive summary of the included functions:

- `alpha.calc` function to compute multiplicative variance effects
- `alpha.cond` function to do conditional analysis of variance effects
- `kappa_calc` function to compute correlation effects (gt/pheno/pheno)
- Correlation calibration
- `Var.assoc` Testing variance scores associations with data
- Dominance effect model implementation
- Interaction effect model (gt/gt/pheno) (genotype interaction)
  implementation
  - Pairwise genotype interaction implementation for list of genotypes
- Interaction effect model (gt/pheno/pheno) (environment interaction)
  implementation
  - Environment interaction cross of lists of phenotypes and genotypes
    (single outcome phenotype)
- Function to create traditional genetic score
- Function to create traditional genetic score with interaction effects
  as well
- Function to create traditional genetic score with interaction effects
  and dominance effects as well
- Function to create variance genetic score
- Summary visualizations
- Histograms by genotype

A more detailed description will be added with further development.

## Installation

You can install the latest version of the package via the `remotes`
package:

``` r
# Use remotes:
remotes::install_github("DecodeGenetics/gnonadd")
```

When the package is on CRAN, you should be able to install it with:

``` r
install.packages("gnonadd")
```
