
#' Organise temporarily needed packages
#' 
#' Attach required packages which may be detached after function run.
#' 
#' Help functions to organise package loading. Some functions require huge annotation packages which shall 
#' not be loaded for the full session. The function \code{attach_package} attaches these packages 
#' for an individual function run and installes them from CRAN and/or Bioconductor if necessary. 
#' \code{detach_package} can be used to detach these packages after function run. DLLs of unloaded packages
#' are removed as well.
#' 
#' @describeIn attach_package Attach and install packages.
#' 
#' @param pkg.cran Character vector with names of CRAN packages to attach
#' @param pkg.bioc Character vector with names of Bioconductor packages to attach
#' @param update_bioc_pkgs logical. When \code{FALSE}, \code{BiocManager::install()} does not attempt to 
#' update old packages. When TRUE, update old packages according to ask. 
#' @param pkg Character vector with names of packages to be detached
#' 
#' @return character vector with package names which have not yet been attached to workspace before call of
#'         \code{attach_package}. This vector may be used to detach packages not needed any more with 
#'         \code{detach_package}. \code{detach_package} invisibly returns the set of identified stray DLLs 
#'         reported by \code{R.utils::gcDLLs}.
#' 
#' @author Frank Ruehle
#' 
#' @export attach_package 
#' @export detach_package  

attach_package <- function(pkg.cran=NULL, pkg.bioc=NULL, update_bioc_pkgs=F) {
 
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
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    for (p in 1:length(pkg.bioc)) {  
      if (!is.element(pkg.bioc[p], installed.packages()[,1])) {
        BiocManager::install(pkg.bioc[p], update = update_bioc_pkgs)}
      require(pkg.bioc[p], character.only = TRUE) 
    }
  }

  pkg2detach <- c(unique(not.yet.attached.pkg.cran, not.yet.attached.pkg.bioc))
   
  return(pkg2detach) 
} 


#' @describeIn attach_package Remove packages not needed anymore
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

# It is possible to have multiple versions of a package loaded at once (for example, if you have a development 
# version and a stable version in different libraries). To detach guarantee that all copies are 
# detached, use this function.

