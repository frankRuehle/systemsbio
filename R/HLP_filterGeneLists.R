
#' Filter gene lists
#' 
#' Gene list are filtered for selected columns considering desired direction and/or transformation.
#' 
#' Function filters dataframe or character directing to dataframe for up to two categories (e.g. p-value and foldchange).
#' Columns for filter criteria must be included in dataframe \code{genes}.
#' Values in filter criteria may be transformed e.g. to absolute values by \code{filtercat.function = abs}, 
#' which is needed for log foldchanges. No tranformation if \code{filtercat.function = identity}.
#' 
#' 
#' @param genes dataframe or character with path directing to table with genelist and columns to filter for.
#' @param newheader NULL if \code{genes} already supplied with header. Character vector with new header otherwise.
#' @param filtercat1 column name of first category to filter \code{genes} (e.g. p-values). Skipped if NULL.
#' @param filtercat1.decreasing (boolean): direction to order and filter \code{filtercat1}
#' @param filtercat1.function select transforming function for \code{filtercat1} (no quotes). 
#'          e.g. \code{abs} for absolute values, \code{identity} for no transformation
#' @param filtercat1.threshold Threshold for \code{filtercat1} or 'top123' for (123) top Hits
#' @param filtercat2 column name of second category to filter \code{genes} (e.g. effect size). Skipped if NULL.
#' @param filtercat2.decreasing (boolean): direction to order and filter \code{filtercat2}
#' @param filtercat2.function select transforming function for \code{filtercat2} (no quotes). E.g. abs for foldchanges.
#' @param filtercat2.threshold Threshold for \code{filtercat2} or 'top123' for top Hits
#' 
#' @return dataframe filtered for desired criteria and sorted for \code{filtercat1} (if not Null).
#' 
#' @author Frank Ruehle
#' 
#' @export 


filterGeneLists <- function(genes,
                            newheader=NULL,
                            filtercat1 = "adj.P.Val",
                            filtercat1.decreasing = FALSE,
                            filtercat1.function = abs,
                            filtercat1.threshold = 0.05,
                            filtercat2 = "logFC",
                            filtercat2.decreasing = TRUE,
                            filtercat2.function = abs,
                            filtercat2.threshold = log2(1.5)) {
  
  
  # read file if 'genes' is character string with file path
  if(is.character(genes)) {
    cat("\n\nReading:", genes, "\n")
    genes <- read.table(genes, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
  }
  if(!is.null(newheader)){# if 'newheader' is defined, it is used as names(genes)
    cat("\nNew header added to input file:", newheader, "\n")
    names(genes) <- newheader
  }  
  
  
  
  # apply filtercat2  
  if (!is.null(filtercat2)) {
    if(filtercat2 %in% names(genes)) {
      genes <- genes[order(filtercat2.function(genes[,filtercat2]), decreasing=filtercat2.decreasing),]
      
      if(is.character(filtercat2.threshold) & grepl("top", filtercat2.threshold, ignore.case=T)) { # select top hits
        topthreshold <- min(nrow(genes), as.numeric(sub("top", "", filtercat2.threshold, ignore.case=T)))
        cat("\nTop", topthreshold, "entries selected according to", filtercat2, "(decreasing =",filtercat2.decreasing, ")\n")
        genes <- genes[1:topthreshold,]
      } else {
        cat("\nFiltering Data for", filtercat2, ", decreasing =",filtercat2.decreasing, ", threshold =", filtercat2.threshold, "\n")
        if(filtercat2.decreasing==FALSE) { # less than threshold
          genes <- genes[!is.na(genes[,filtercat2]) & filtercat2.function(genes[,filtercat2]) < filtercat2.threshold, ]
        } else { # greater than threshold
          genes <- genes[!is.na(genes[,filtercat2]) & filtercat2.function(genes[,filtercat2]) > filtercat2.threshold, ]
        }
      } 
    }}
  
  
  # apply filtercat1
  if (!is.null(filtercat1)) {
    if(filtercat1 %in% names(genes)) {
      genes <- genes[order(filtercat1.function(genes[,filtercat1]), decreasing=filtercat1.decreasing),]
      
      if(is.character(filtercat1.threshold) & grepl("top", filtercat1.threshold, ignore.case=T)) { # select top hits
        topthreshold <- min(nrow(genes), as.numeric(sub("top", "", filtercat1.threshold, ignore.case=T)))
        cat("\nTop", topthreshold, "entries selected according to", filtercat1, "(decreasing =",filtercat1.decreasing, ")\n")
        genes <- genes[1:topthreshold,]
      } else {
        cat("Filtering Data for", filtercat1, ", decreasing =",filtercat1.decreasing, ", threshold =", filtercat1.threshold, "\n")
        if(filtercat1.decreasing==FALSE) { # less than threshold
          genes <- genes[!is.na(genes[,filtercat1]) & filtercat1.function(genes[,filtercat1]) < filtercat1.threshold, ] 
          
        } else { # greater than threshold
          genes <- genes[!is.na(genes[,filtercat1]) & filtercat1.function(genes[,filtercat1]) > filtercat1.threshold, ]
        }
      } 
    }}
  
  
  return(genes)
  
} # end function definition 'filterGeneLists'


