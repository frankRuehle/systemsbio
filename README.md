# systemsbio #
Streamlined analysis package for omics-data

This package consists of modularized wrapper functions for multiple genomics analysis packages. It is designed for multi-omics analysis of expression, methylation and/or genotyping data. All modules can also be used individually.

See the package [Vignette](inst/doc/systemsbio-vignette.html) for documentation of the package and example analyses.

## Installation

1. Install (if you haven't already) a working development environment:
    * **Windows**: Install [Rtools](http://cran.r-project.org/bin/windows/Rtools).
    * **Mac**: Install [Xcode](https://itunes.apple.com/de/app/xcode/id497799835).
    * **Linux**: Install a compiler for your distribution. For instance, for Ubuntu this would be `sudo apt-get install     r-base-dev`. Further instructions can be found at [CRAN](http://cran.r-project.org/bin/linux).

2. Install (if you haven't already) the `devtools` package via CRAN:

        install.packages(c("devtools", "rstudioapi"))
  
3. Install `systemsbio` from GitHub via `devtools`:

        devtools::install_github("systemsbio", username="frankRuehle", build_vignettes=TRUE)


## First steps
Start by loading the `systemsbio` package and read the excellent documentation.

    library(systemsbio)
    vignette(systemsbio)


## Author
Frank Ruehle