####### Organise packages

##### Attach packages
# either from CRAN or Bioconductor
attach_package <- function(pkg.cran=NULL, pkg.bioc=NULL, source.bioc="http://bioconductor.org/biocLite.R") {
  
  # install CRAN libraries
  if(!is.null(pkg.cran)) {
    for (p in 1:length(pkg.cran)) {  
      if (!is.element(pkg.cran[p], installed.packages()[,1])) {
        install.packages(pkg.cran[p], dependencies=TRUE)}
      require(pkg.cran[p], character.only = TRUE) 
    }
  }
  
  # install Bioconductor libraries
  if(!is.null(pkg.bioc)) {
    source(source.bioc)
    for (p in 1:length(pkg.bioc)) {  
      if (!is.element(pkg.bioc[p], installed.packages()[,1])) {
        biocLite(pkg.bioc[p])}
      require(pkg.bioc[p], character.only = TRUE) 
    }
  }
  
} # end function definition

