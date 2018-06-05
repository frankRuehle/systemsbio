#' Processing and data transformation of DESeqDataSet 
#'
#' Process count matrix of DESeqDataSet and optionally convert it to a SummarizedExperiment object.
#' 
#' The count matrix of the input object is filtered for low abundancy transcripts, normalized for 
#' library size factor (size factors are first calculated if necessary) and is transformed
#' according to the function given in \code{transform.function}. Parameters for the selected 
#' \code{transform.function} can be passed via the \code{...} -argument (see help for the 
#' respective function). The \code{DESeq2} package recommends \code{varianceStabilizingTransformation}
#' for further processing of the read count matrix. If you prefer to analyse reads with limma, you may 
#' apply \code{voom} transformation as recommended in the \code{limma} documentation. The \code{voom}
#' transformation can be complemented by the same between-array normalization methods as would be used 
#' for microarrays (e.g. quantile normalization) via the \code{normalize.method} parameter of \code{voom}.
#' The function will return either the processed input object or an \code{SummarizedExperiment} object
#' if \code{return_SummarizedExperiment} is set to \code{TRUE}. In this case, the design parameter of 
#' the initial \code{DESeqDataSet} is omitted.
#'
#' @param dds \code{DESeqDataSet}
#' @param min_rowsum numeric. Minimum rowsum of countmatrix. All rows with \code{rowSums < min_rowsum}
#' are removed from the count matrix. 
#' @param norm_sizeFactors logical. If \code{TRUE} normalize for library size. This may already be included
#' in the data transformation (e.g. for \code{rlog} and \code{varianceStabilizingTransformation}).
#' @param transform.function function name to be used for data transformation, 
#' e.g. \code{DESeq2::rlog} or \code{DESeq2::rlogTransformation} for log2-transformation, 
#' \code{DESeq2::varianceStabilizingTransformation}, \code{fpkm}, \code{fpm}, \code{voom}. 
#' @param return_SummarizedExperiment logical. If \code{True} function returns a \code{SummarizedExperiment} object. 
#' If FALSE the processed input object is returned (may be \code{matrix}, \code{DESeqDataSet} or \code{DESeqTransform}).
#' @param ... further parameter for \code{transform.function}
#' 
#' @return \code{SummarizedExperiment} or processed input object (\code{matrix}, \code{DESeqDataSet} or 
#' \code{DESeqTransform} depending on applied data transformation)
#' 
#' @author Frank Ruehle
#' 
#' @export 



process_dds <- function(dds,
                       min_rowsum = 10,
                       norm_sizeFactors = FALSE,
                       transform.function = NULL,
                       return_SummarizedExperiment = FALSE,
                       ...) {
  

  # load required packages. Packages are not detached afterwards
  pkg.bioc <- c("DESeq2", "SummarizedExperiment", "limma")
  pks2detach <- attach_package(pkg.bioc=pkg.bioc)
  
  if(is.null(sizeFactors(dds))) {dds <- DESeq2::estimateSizeFactors(dds)}
  
  
    # Pre-filtering
    cat("\nApply pre-filtering for rowSums >=", min_rowsum)
    dds <- dds[rowSums(counts(dds)) >= min_rowsum,]
    
    # countmat <- assay(dds)
   
    if(norm_sizeFactors) {
        if(any(grepl(deparse(substitute(transform.function)), c("varianceStabilizingTransformation", "rlog", "fpkm", "fpm")))) {
          warning("Your data transformation already included normalization of size factors. 
                   Therefore, norm_sizeFactors is set to FALSE!\n")
        } else {
      cat("\nApply normalization by library size factors")
      assay(dds) <- sweep(assay(dds), 2, DESeq2::sizeFactors(dds),"/")
        }
    }

 
    # prepare SummarizedExperimen object
      sumexp <-  SummarizedExperiment::SummarizedExperiment(assay= assay(dds), # count data
                                                            colData=colData(dds), # sample data
                                                            rowData=rowData(dds), # feature data
                                                            metadata=metadata(dds)) # experiment data
     
    
 
    # transform input object
    if(!is.null(transform.function)) {
      cat("\napply data transformation:", deparse(substitute(transform.function)))
 
      if(any(grepl(deparse(substitute(transform.function)), c("voom")))) { # if function need matrix as imput
        dds <- transform.function(assay(dds), ...) } else {
        
        dds <- transform.function(dds, ...)
      }

      if(class(dds) %in% c("DESeqDataSet", "DESeqTransform")) { # adapt assay matrix of SummarizedExperiment object
        assay(sumexp) <- assay(dds) } else {
          if(class(dds) %in% c("EList")) { # e.g. for voom
            assay(sumexp) <- limma::getEAWP(dds)$exprs} else {
              assay(sumexp) <- dds # if transformation output is matrix , e.g. fpkm
          }
        }
    }

       
   # get return object          
   if(return_SummarizedExperiment) {
     result <- sumexp
   } else {
     result <- dds
   }

  cat("\nreturn", class(result), "object.")
  return(result)

} # end function definition

