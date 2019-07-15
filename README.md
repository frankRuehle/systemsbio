# systemsbio #
Streamlined analysis package for omics-data

This package consists of modularized wrapper functions for multiple genomics analysis packages. 
It is designed for multi-omics analysis of expression, methylation and/or genotyping data. 
All modules can be used individually.

See the package Vignette [(version without figures)](./vignettes/systemsbio.Rmd) for documentation of the package and example analyses. The full version of the vignette <a href=".doc/systemsbio.html">here</a> is available by calling `vignette("systemsbio")` as soon as the package is installed. A flow chart of the package content is given [here](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=Pipeline_Systems_Biology.html#Uhttps%3A%2F%2Fdrive.google.com%2Fuc%3Fid%3D12uDryY6msteXpXoty8qFtpvZjfTxz6MR%26export%3Ddownload).

## Installation

1. Install (if you haven't already) the `devtools` package via CRAN (on Windows you will also need Rtools):

        install.packages(c("devtools", "rstudioapi"))
        library(devtools)

2. Install `systemsbio` from GitHub via `devtools`:

        devtools::install_github("frankRuehle/systemsbio", build_vignettes=TRUE)


## First steps
Start by loading the `systemsbio` package and read the package vignette.

    library(systemsbio)
    vignette("systemsbio")


## Author
Frank Ruehle