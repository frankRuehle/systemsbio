
#' Basic gene annotation
#'
#' Annotate gene symbols and EntrezIDs using an organism-specific annotation package. 
#' 
#' This function adds annotation of either gene symbol or EntrezID to \code{data}
#' using the corresponding annotation package of the organism specified in \code{org}.
#' 
#' 
#' @param data dataframe or character with filepath to data to be loaded.
#' @param Symbol.column character with column name of Gene Symbols in \code{data}.
#'                if NULL, SYMBOLs are derived from \code{Entrez.column} using the annotation package for \code{org}.
#' @param Entrez.column character with column name of ENTREZ IDs in \code{data}.
#'                if NULL, ENTREZ IDs are derived from \code{Symbol.column} using the annotation package for \code{org}.
#'                Either \code{Symbol.column} or \code{Entrez.column} must be specified.
#' @param org character with species name ("human", "mouse", "rat").
#' 
#' 
#' @return dataframe with added annotation column.
#' 
#' @author Frank Ruehle
#' 
#' @export 


basicAnno <- function(data, 
                      Symbol.column = "SYMBOL",
                      Entrez.column = "ENTREZID",
                      org ="human"
                      ) {



  # Load genomic annotation package (includes OrgDb, TxDb and GODb)
  annopkg <- switch(org, human = "org.Hs.eg.db", # "Homo.sapiens", 
                    hsapiens = "org.Hs.eg.db", # "Homo.sapiens", 
                    mouse = "org.Mm.eg.db", # "Mus.musculus",
                    rat= "org.Rn.eg.db") # "Rattus.norvegicus")
  pks2detach <- attach_package(pkg.bioc= c(annopkg, "GenomicFeatures", "rtracklayer"), pkg.cran = "plyr")
  
  
  
  ### read file if 'data' is character string with file path 
  # column names as they are, not checked/modified for valid R names!
    if(is.character(data) & length(data)==1) {
    cat("\n\nReading expression data:", data, "\n")
      data <- read.table(data, header=T, sep="\t", na.strings = c("", " ", "NA"), check.names =F) 
  }
  
  
  
 # check if dedicated annotation columns exist. Assign NULL otherwise.
  if(!is.null(Symbol.column)) {
    if(!(Symbol.column %in% colnames(data))) {
      warning(cat("\n", Symbol.column, "not found in colnames of data object!\n"))
      Symbol.column <- NULL}
  }
    
  if(!is.null(Entrez.column)) {
    if(!(Entrez.column %in% colnames(data))) {
      warning(cat("\n", Entrez.column, "not found in colnames of data object!\n"))
      Entrez.column <- NULL}
  }
  
  # check if no annotation columns exists.
  if(is.null(c(Symbol.column, Entrez.column))) {
    stop("\nEither 'Symbol.column' or 'Entrez.column' must be specified!\n")
  }
  
   
  
  
  ## convert SYMBOLs to ENTREZ IDs if necessary (no additional rows introduced)
  if(is.null(Entrez.column) & !is.null(Symbol.column)) {
    cat(paste("\nConverting SYMBOLs to ENTREZIDs using",annopkg,"annotation package.\n"))
    new.entrezids <- select(get(annopkg), keys = as.character(data[,Symbol.column]), keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")) 
    names(new.entrezids)[names(new.entrezids)=="SYMBOL"] <- Symbol.column
    data <- plyr::join(data, new.entrezids, by=Symbol.column, type="left", match="first") 
  }
  
  ## convert ENTREZ IDs to SYMBOLs if necessary (no additional rows introduced)
  if(is.null(Symbol.column) & !is.null(Entrez.column)) {
    cat(paste("\nConverting ENTREZ IDs to SYMBOLs using",annopkg,"annotation package.\n"))
    new.symbols <- select(get(annopkg), keys = as.character(data[,Entrez.column]), keytype = "ENTREZID", columns = c("SYMBOL", "ENTREZID"))
    names(new.symbols)[names(new.symbols)=="ENTREZID"] <- Entrez.column
    data <- plyr::join(data, new.symbols, by=Entrez.column, type="left", match="first") 
  }


  
  # # Detaching libraries not needed any more
  # detach_package(unique(pks2detach))
  
  return(data)
  
} # end of function definition

