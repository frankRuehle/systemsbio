#' Transform DESeqDataSet to SummarizedExperiment
#'
#' Process count matrix of DESeqDataSet to SummarizedExperiment to be used as expression matrix.
#' 
#' The count matrix of the input object is filtered for low abundancy transcripts, normalized for 
#' library size factor (size factors are first calculated if necessary) and is transformed
#' according to the function given in \code{transform.function}. The design parameter of the initial \code{DESeqDataSet}
#' is omitted
#'
#' @param dds \code{DESeqDataSet}
#' @param min_rowsum numeric. Minimum rowsum of countmatrix. All rows with \code{rowSums < min_rowsum}
#' are removed from the count matrix. 
#' @param transform.function function name to be used for data transformation, 
#' e.g. \code{DESeq2::rlog} or \code{DESeq2::rlogTransformation} for log2-transformation or \code{DESeq2::varianceStabilizingTransformation}. 
#' @param norm_sizeFactors logical. If TRUE normalize for library size. This may already be included
#' in the data transformation (e.g. for \code{rlog} and \code{varianceStabilizingTransformation}).
#' @param export return either a \code{"SummarizedExperiment"} or a \code{"DESeq"} object 
#' (i.e. either \code{"DESeqDataSet"} or \code{"DESeqTransform"} object).
#' @param ... further parameter for \code{transform.function}
#' 
#' @return \code{SummarizedExperiment} or \code{DESeqDataSet} object
#' 
#' @author Frank Ruehle
#' 
#' @export 



process_dds <- function(dds,
                       min_rowsum = 10,
                       transform.function = NULL,
                       norm_sizeFactors = FALSE,
                       export="SummarizedExperiment",
                       ...) {
  

  # load required packages. Packages are not detached afterwards
  pkg.bioc <- c("DESeq2", "SummarizedExperiment")
  pks2detach <- attach_package(pkg.bioc=pkg.bioc)
  
  if(is.null(sizeFactors(dds))) {dds <- DESeq2::estimateSizeFactors(dds)}
  
  
    # Pre-filtering
    cat("\nApply pre-filtering for rowSums >=", min_rowsum)
    dds <- dds[rowSums(counts(dds)) >= min_rowsum,]
    
    # countmat <- assay(dds)
  

    # transform input object
    if(!is.null(transform.function)) {
       cat("\napply data transformation")
       dds <- transform.function(dds, ...)
    }
  
    if(norm_sizeFactors) {
    cat("\nApply normalization by library size factors")
    assay(dds) <- sweep(assay(dds), 2, DESeq2::sizeFactors(dds),"/")
    }
  
  if(export=="SummarizedExperiment") {
    dds <-  SummarizedExperiment::SummarizedExperiment(assay=assay(dds), # count data
                                                        colData=colData(dds), # sample data
                                                        rowData=rowData(dds), # feature data
                                                        metadata=metadata(dds)) # experiment data
    } 
    
  cat("\nreturn", class(dds), "object.")
  return(dds)

} # end function definition

### to be added: DESeq::fpm, DESeq::fpkm
