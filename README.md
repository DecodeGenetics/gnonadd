
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnonadd

<!-- badges: start -->
<!-- badges: end -->

`gnonadd` is a package accompanying the paper [Complex effects of
sequence variants on lipid levels and coronary artery
disease](https://www.sciencedirect.com/science/article/pii/S0092867423009017?dgcid=coauthor)
published in Cell September 2023. The package is intended to properly
document the conducted analysis and aid researchers in studying various
non-additive models.

## What is in the package?

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

Please refer to the documentation for examples with simulated data.

## Installation

You can install the latest version of the package via the `remotes`
package:

``` r
# Use remotes:
remotes::install_github("DecodeGenetics/gnonadd")
```

The current version on CRAN can be installed with:

``` r
install.packages("gnonadd")
```

## Citing this package

For citing this package, please use the following source:

``` r
citation("gnonadd")
#> 
#> To cite gnonadd in publications use:
#> 
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Complex effects of sequence variants on lipid levels and coronary artery disease},
#>     author = {Audunn S. Snaebjarnarson et al.},
#>     journal = {Cell},
#>     year = {2023},
#>     volume = {186},
#>     issue = {19},
#>     pages = {4085-4099.e15},
#>     url = {https://www.sciencedirect.com/science/article/pii/S0092867423009017},
#>     doi = {https://doi.org/10.1016/j.cell.2023.08.012},
#>   }
```
