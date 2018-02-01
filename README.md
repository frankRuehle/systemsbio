# systemsbio #
Streamlined analysis package for omics-data

This package consists of modularized wrapper functions for multiple genomics analysis packages. 
It is designed for multi-omics analysis of expression, methylation and/or genotyping data. 
All modules can be used individually.

See the package [Vignette](./vignettes/systemsbio.Rmd) for documentation of the package and example analyses.

## Installation

1. Install (if you haven't already) the `devtools` package via CRAN:

        install.packages(c("devtools", "rstudioapi"))
  
2. Install `systemsbio` from GitHub via `devtools`:

        devtools::install_github("systemsbio", username="frankRuehle", build_vignettes=TRUE)


## First steps
Start by loading the `systemsbio` package and read the package vignette.

    library(systemsbio)
    vignette("systemsbio")


## Author
Frank Ruehle