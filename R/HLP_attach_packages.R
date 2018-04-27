

#' Organise temporarily needed packages
#' 
#' Attach required packages which may be detached after function run.
#' 
#' Help functions to organise package loading. Some functions require huge annotation packages which shall 
#' not be loaded for the full session. These packages are attached for an individual function and detached 
#' after function run. Packages are installed if necessary.
#' 
#' @param pkg.cran Character vector with names of CRAN packages to attach
#' @param pkg.bioc Character vector with names of Bioconductor packages to attach
#' @param suppressUpdatesBioc logical. If \code{FALSE}, the user is asked whether old bioconductor packages should 
#' be updated. If \code{TRUE}, the user is not prompted to update old packages. 
#' @param source.bioc character with Bioconductor source site (e.g. \code{http://bioconductor.org/biocLite.R})
#' @param pkg Character vector with names of packages to be detached
#' 
#' @return character vector with package names which have not yet been attached to workspace before call of
#'         \code{attach_package}. This vector may be used to detach packages not needed any more with \code{detach_package}. 
#' 
#' @author Frank Ruehle
#' 
#' @export attach_package 
#' @export detach_package  




##### Attach packages
# either from CRAN or Bioconductor
attach_package <- function(pkg.cran=NULL, pkg.bioc=NULL, suppressUpdatesBioc=T, source.bioc="http://bioconductor.org/biocLite.R") {
 
  not.yet.attached.pkg.cran <- NULL
  not.yet.attached.pkg.bioc <- NULL
  
  # install CRAN libraries
  if(!is.null(pkg.cran)) {
    not.yet.attached.pkg.cran <- pkg.cran[!(pkg.cran %in% loadedNamespaces())]
    for (p in 1:length(pkg.cran)) {  
      if (!is.element(pkg.cran[p], installed.packages()[,1])) {
        install.packages(pkg.cran[p], dependencies=TRUE)}
      require(pkg.cran[p], character.only = TRUE) 
    }
  }
  
  # install Bioconductor libraries
  if(!is.null(pkg.bioc)) {
    not.yet.attached.pkg.bioc <- pkg.bioc[!(pkg.bioc %in% loadedNamespaces())]
    source(source.bioc)
    for (p in 1:length(pkg.bioc)) {  
      if (!is.element(pkg.bioc[p], installed.packages()[,1])) {
        biocLite(pkg.bioc[p], suppressUpdates = suppressUpdatesBioc)}
      require(pkg.bioc[p], character.only = TRUE) 
    }
  }

  pkg2detach <- c(unique(not.yet.attached.pkg.cran, not.yet.attached.pkg.bioc))
   
  return(pkg2detach) 
} # end function definition


##### Detach packages
# It is possible to have multiple versions of a package loaded at once (for example, if you have a development 
# version and a stable version in different libraries). To detach guarantee that all copies are 
# detached, use this function.
detach_package <- function(pkg)
{
  for(p in pkg) {
   
    search_item <- paste("package", p, sep = ":")
    while(search_item %in% search())
    {
      detach(search_item, unload = TRUE, character.only = TRUE, force = TRUE)
    }
  } # end for-loop
  R.utils::gcDLLs() # Identifies and removes DLLs of packages already unloaded
} # end function definition


