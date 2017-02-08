

#' Organise temporarily needed packages
#' 
#' Attach required packages which may be detached after function run.
#' 
#' Help functions to organise packege loading. Some functions require huge annotation packages which shall 
#' not be loaded for the full session. These packages are attached for an individual function and detached 
#' after function run.
#' 
#' 
#' @param pkg.cran Character vector with names of CRAN packages to attach
#' @param pkg.bioc Character vector with names of Bioconductor packages to attach
#' @param source.bioc character with Bioconductor source site (e.g. \code{http://bioconductor.org/biocLite.R})
#' @param pkg Character vector with names of packages to be detached
#' @param character.only a logical indicating whether name can be assumed to be a character string.
#' 
#' @author Frank Ruehle
#' 
#' @export attach_package 
#' @export detach_package  




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


##### Detach packages
# It is possible to have multiple versions of a package loaded at once (for example, if you have a development 
# version and a stable version in different libraries). To detach guarantee that all copies are 
# detached, use this function.
detach_package <- function(pkg, character.only = FALSE)
{
  for(p in pkg) {
    if(!character.only)
    {
      p <- deparse(substitute(p))
    }
    search_item <- paste("package", p, sep = ":")
    while(search_item %in% search())
    {
      detach(search_item, unload = TRUE, character.only = TRUE)
    }
  } # end for-loop
} # end function definition


