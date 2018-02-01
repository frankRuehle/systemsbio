#' Reading methylation data from Illumina array HumanMethylation450k
#' 
#' Create an RGChannelSet object from Illumina IDAT files and add phenotype data.
#' 
#' This functions uses the \code{minfi} package to load methylation data from Illumnia 
#' Human Methylation 450k microarrays. 
#' If the \code{SampleSheetRegExp} argument is NULL, the function finds all two-color IDAT files in the 
#' directory given by \code{baseDir} including all subdirectories. Two-color IDAT files consist of pairs of 
#' files with names ending in \code{_Red.idat} or \code{_Grn.idat}. If \code{SampleSheetRegExp} is not NULL,
#' the corresponding sample sheet is identified via regular expression (file must be located in base) 
#' and all listed arrays are loaded. 
#' Susequently, additional phenotype data given in \code{covarfile} is added to the object.
#' 
#' @param baseDir character with path to base directory containing IDAT files for all arrays.
#' @param SampleSheetRegExp dataframe with sample sheet or character with regular expression for identifying 
#' the corresponding sample sheet within the base directory or its subdirectories. If NULL, no sample information file is loaded.
#' @param sampleColumn Name of sample column in sample sheet
#' @param covarfile character with filepath to a covariates file if given. NULL otherwise.
#' @param covarsampleID character with name of sample column in (optional) covariates file.   
#' 
#'    
#' @return RGChannelSet object
#' 
#' @author Frank Ruehle
#' 
#' @export
#' 
#' 
#' @note methylation data classes (from \code{minfi} manual) :
#' \itemize{
#'   \item RGChannelSet: raw data from the IDAT files; this data is organized at the probe (not CpG locus) level. 
#'   This data has two channels: Red and Green. This data type contains control probe data.
#'   \item MethylSet: data organized by the CpG locus level, but not mapped to a genome. 
#'   This data has two channels: Meth (methylated) and Unmeth (unmethylated).
#'   \item RatioSet: data organized by the CpG locus level, but not mapped to a genome. The data has at least 
#'   one of two channels: Beta and/or M (logratio of Beta). It may optionally include a CN channel (copy number).
#'   \item GenomicMethylSet: like a MethylSet, but mapped to a genome.
#'   \item GenomicRatioSet: like a RatioSet, but mapped to the genome.
#'   }





read_Illu_meth450k <- function(baseDir, 
                   SampleSheetRegExp = ".csv", 
                   sampleColumn = "Sample_Name", 
                   covarfile = NULL, 
                   covarsampleID  = "ID"  
                   ) {

  

  # load required libraries
  pkg.cran <- c("plyr")
  pkg.bioc <- c("minfi")  
  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  
  
# Loading SampleSheet
SSheetUnique <- FALSE
if(is.null(SampleSheetRegExp)) {
  cat("\nNo SampleSheet supplied. Methylation data is loaded without SampleSheet information.\n")
  SampleSheetMT <- NULL
        } else {

      if (is.data.frame(SampleSheetRegExp)) {
        SampleSheetMT <- SampleSheetRegExp # no loading of sample sheet necessary.
        SSheetUnique <- TRUE
        } else {
          
      cat("\nLooking for SampleSheet in",baseDir, "with pattern:",SampleSheetRegExp, "\n")
      LookupSSheet <- grep(SampleSheetRegExp, list.files(baseDir), value=T)

        if(length(LookupSSheet)>1) {
            cat("\n\n", length(LookupSSheet), "possible SampleSheets found by regular expression. Redefine expression:", SampleSheetRegExp)
            cat("\n", LookupSSheet, "\n")
            cat("\nMethylation data is loaded without SampleSheet information.\n")
            SampleSheetMT <- NULL
              } else {
              
              if(length(LookupSSheet)==0) {
                cat("\nno SampleSheet found by regular expression. Redefine expression:", SampleSheetRegExp, "\n")
                cat("\nMethylation data is loaded without SampleSheet information.\n")
                SampleSheetMT <- NULL
                  } else {
                  
                    cat("\nLoading", LookupSSheet, "\n")
                    SSheetUnique <- TRUE
                    SampleSheetMT <- read.metharray.sheet(base = baseDir, pattern = SampleSheetRegExp)   # read Sample Sheet  
                          }
                      }
            }
      }



# Read IDAT-files to RGChannelSet object:
  cat("\nReading methylation data from",baseDir,"\n")
  if(!is.null(SampleSheetMT)) {baseDir <- NULL} # don't give base AND targets argument at once!
  RGset <- read.metharray.exp(base = baseDir, targets = SampleSheetMT, recursive = T)
  
  
  # add phenotype data
  if(SSheetUnique) {
    sampleNames(RGset) <- pData(RGset)[,sampleColumn] # row order is congruent
  
    # load and append covar files
      if (!is.null(covarfile)) {
        cat("\nAdding covariates from", covarfile, "\n")
        covar <- read.delim(file=covarfile, header=T, sep="\t", stringsAsFactors =T, 
                           na.strings = c("NA", "-9", "", " "), fileEncoding="latin1") 
        covar[,covarsampleID] <- sub("@.*$", "", covar[,covarsampleID])  # remove possible "@cohortID"-suffix from SNPZone
        covar[,covarsampleID] <- sub(" ", "", covar[,covarsampleID])    # remove possible whitespace from sample names
    
        # join() is used for merging because merge() does not preserve the roworder 
        names(covar)[names(covar)==covarsampleID] <- sampleColumn
        pData(RGset) <- plyr::join(pData(RGset), covar, by=sampleColumn, type="left")
        
  
      # after joining, sampleNames are set to numbers, therefore renaming.
      sampleNames(RGset) <- pData(RGset)[,sampleColumn]
  
      }
  }
  

  # Detaching libraries not needed any more
  detach_package(unique(pks2detach))

show(RGset)

return(RGset)

}

