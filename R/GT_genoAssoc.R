#' Association analysis with GenABEL
#'
#' \code{genoAssoc} takes genotype data in GenABEL gwaa format and performs association analysis.
#' 
#' The qtscore()-function from GenABEL package is used for association analysis of a gwaa object.
#' The trait of interest is given in \code{trait.name} either as single trait (e.g. "affection")
#' or as formula (e.g. "affection01 ~ Age + sex") if association shall be adjusted for covariates.
#' All traitnames in \code{trait.name} must exist within the phenotype data of the gwaa object. 
#' Use phdata() for checking which trait names ara avalable. If \code{times} is set to more than 1,
#' empirical p-values are calculated using \code{times} permutations.
#' A table and a Manhattan plot with association results are stored in \code{projectfolder}.
#'  
#'
#' @param gwaa gwaa object from GenABEL
#' @param projectfolder character containing path to output folder (will be generated if not existing).
#' @param projectname character used as suffix for output files.
#' @param trait.name character indicating trait(s) of interest.
#' @param trait.type character with data type of analysed trait. Either "gaussian", "binomial" or "guess" 
#' (later option guesses trait type).
#' @param times If more than one, the number of permutations to be used for empirical p-values.
#' @param quiet boolean. Do not print warning messages.
#' @param TopHitsReported numeric with number of "top" hits to describe. Use nsnps(gwaa.object) for
#' all SNPs
#' 
#' @return Association results as scan.gwaa object. Intermediary results and plots are stored in \code{projectfolder} 
#' as side effects.
#' 
#' @author Frank Ruehle
#' 
#' @export genoAssoc


genoAssoc <- function(gwaa, 
                    projectfolder = "GT/Assoc_GenAbel",
                    projectname = NULL, 
                    trait.name = "affection01 ~ Age + sex", 
                    trait.type = "binomial", # "gaussian" or "binomial" 
                    times = 1,   
                    quiet=FALSE, 
                    TopHitsReported = 1000
                    ) {

  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }  

  projectname <- if (!is.null(projectname) && !grepl("_$", projectname)) {paste0(projectname, "_")} else {""}
  

  ## install/load required packages from CRAN and Bioconductor
  pkg.bioc <- NULL
  pkg.cran <- c("GenABEL")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
 ###### association analysis (additive model)  
 attach(phdata(gwaa))
   
   if(grepl("~", trait.name)) {traits <- as.formula(trait.name)
                              } else {traits <- get(trait.name)}
   
      gwaa.qt <- qtscore(formula= traits, data= gwaa, quiet= quiet, trait.type= trait.type, times= times, details = TRUE)

 detach(phdata(gwaa))
 

  cat("\n\nlambda = ")
  print(lambda(gwaa.qt))
  
  if(times == 1) {cat("\nAssociation results:\n\n")
    } else {cat(paste0("\nEmpirical association results: (", times, ") permutations:\n\n"))}
  
  print(gwaa.qt)
  AssocResults <- descriptives.scan(gwaa.qt, sort="Pc1df", top=TopHitsReported)
  AssocResults <- data.frame(SNP=rownames(AssocResults), AssocResults)
  
  result.filename <- file.path(projectfolder, paste0(projectname, "AssocResults.txt"))
  cat(paste("\nWrite first", TopHitsReported, "results to"), result.filename)  
  write.table(AssocResults, row.names=F, quote=F, sep="\t", file=result.filename)
  
  manhattan.filename <- file.path(projectfolder, paste0(projectname, "ManhattanPlot.png"))
  cat(paste("\nWrite Manhattan plot to"), manhattan.filename)  
  png(filename=manhattan.filename, width = 297, height = 210, units = "mm", res=600)
  
    plot(gwaa.qt, df="Pc1df",col=c("black", "gray"), 
         main=paste0("Manhattanplot, GC corrected p-values: ", trait.name, if(times>1) {paste(", ", times, " permutations.")}))
    abline(h= -log10(5*1e-8), col="red") # significance line
    abline(h= -log10(1e-5), col="blue") # suggestive significance line
  dev.off()


return(gwaa.qt)


}