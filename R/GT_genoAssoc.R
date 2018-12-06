#' Association analysis with GenABEL
#'
#' \code{genoAssoc} takes genotype data in GenABEL \code{gwaa} format and performs association analysis.
#' 
#' GenABEL package is used for association analysis of a \code{gwaa} object.
#' The trait of interest is given in \code{trait.name} either as single trait (e.g. \code{affection})
#' or as formula (e.g. \code{affection01 ~ Age + sex}) if association shall be adjusted for covariates.
#' All traitnames in \code{trait.name} must exist within the phenotype data of the \code{gwaa} object. 
#' Use \code{phdata(gwaa)} for checking which trait names are avalable. 
#' Either linear or logistic regression or fast score test can be used for association as given in \code{fun.assoc}. 
#' In case of regression, the mode of inheritance is given in \code{gtmode}. 
#' If the fast score test from GenABELs \code{qtscore}-function is used without covariates, the function is equivalent
#' to Armitage TREND Test. If covariates are used, the trait is analysed using \code{GLM}. 
#' Set \code{times} to more than 1, to calculate empirical p-values for the score test using \code{times} permutations.
#' A table and a Manhattan plot with association results are stored in \code{projectfolder}.
#'  
#'
#' @param gwaa \code{gwaa} object from GenABEL
#' @param projectfolder character containing path to output folder (will be generated if not existing).
#' @param projectname character used as suffix for output files.
#' @param trait.name character indicating trait(s) of interest.
#' @param trait.type character with data type of analysed trait. Either "gaussian", "binomial" or "guess" 
#' (latter option guesses trait type).
#' @param fun.assoc character with association function from GenABEL package to be used. 
#' Either \code{mlreg} for linear or logistic regression or \code{qtscore} for fast score test.
#' @param gtmode character with mode of inheritance if regression is used. 
#' Either \code{additive}, \code{dominant}, \code{recessive} or \code{overdominant}. 
#' @param times If more than one, the number of permutations to be used for calculationg empirical p-values.
#' Relevant for score test only.
#' @param quiet logical Do not print warning messages.
#' @param TopHitsReported numeric with number of top hits to report. Use nsnps(gwaa.object) for
#' all SNPs.
#' 
#' @return \code{scan.gwaa} object with Association results. Intermediary results and plots are stored 
#' in \code{projectfolder} as side effects.
#' 
#' @author Frank Ruehle
#' 
#' @export genoAssoc





genoAssoc <- function(gwaa, 
                    projectfolder = "GT/Assoc_GenAbel",
                    projectname = NULL, 
                    trait.name, 
                    trait.type = "binomial",  
                    fun.assoc = "mlreg", 
                    gtmode = "additive", 
                    times = 1,   
                    quiet=FALSE, 
                    TopHitsReported = 1000
                    ) {

  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }  

  projectname <- if (!is.null(projectname) && !grepl("_$", projectname)) {paste0(projectname, "_")} else {""}
  

  ## install/load required packages from CRAN and Bioconductor
  pkg.bioc <- NULL
  pkg.cran <- c("GenABEL")
  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
 ###### association analysis   
 attach(phdata(gwaa))
   
   if(grepl("~", trait.name)) {traits <- as.formula(trait.name)
                              } else {traits <- get(trait.name)}
   
    if (fun.assoc == "qtscore") { # SCORE TEST
        if(times == 1) {cat("\nAssociation results from fast score test:\n\n")
        } else {cat(paste0("\nEmpirical association results from fast score test: (", times, ") permutations:\n\n"))}
        result.assoc <- qtscore(formula= traits, data= gwaa, quiet= quiet, trait.type= trait.type, times= times, details = TRUE)
    }
  
    if (fun.assoc == "mlreg") { # REGRESSION
      cat(paste("\nPerforming regression analysis applying", gtmode, "model\n\n"))
      result.assoc <- mlreg(formula= traits, data= gwaa, trait.type= trait.type, gtmode= gtmode)
    }  
      
 detach(phdata(gwaa))
 
  print(result.assoc)
  
  cat("\n\nlambda = ")
  print(lambda(result.assoc))
  
  
  AssocResults <- descriptives.scan(result.assoc, sort="Pc1df", top=TopHitsReported)
  AssocResults <- data.frame(SNP=rownames(AssocResults), AssocResults)
  
  result.filename <- file.path(projectfolder, paste0(projectname, "AssocResults.txt"))
  cat(paste("\nWrite first", TopHitsReported, "results to"), result.filename)  
  write.table(AssocResults, row.names=F, quote=F, sep="\t", file=result.filename)
  
  manhattan.filename <- file.path(projectfolder, paste0(projectname, "ManhattanPlot.png"))
  cat(paste("\nWrite Manhattan plot to"), manhattan.filename)  
  png(filename=manhattan.filename, width = 297, height = 210, units = "mm", res=600)
  
    plot(result.assoc, df="Pc1df",col=c("black", "gray"), 
         main=paste0("Manhattanplot, GC corrected p-values: ", trait.name, if(times>1) {paste(", ", times, " permutations.")}))
    abline(h= -log10(5*1e-8), col="red") # significance line
    abline(h= -log10(1e-5), col="blue") # suggestive significance line
  dev.off()


return(result.assoc)

}