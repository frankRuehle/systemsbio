
## Description
# Annotate gene symbols and EntrezIDs using an organism-specific annotation package. 

## Usage 
basicAnno <- function(data, 
                      Symbol.column = "SYMBOL",
                      Entrez.column = "ENTREZID",
                      org ="human"
                      ) {


  ## Arguments
  # data: dataframe or character with filepath to data to be loaded.
  # Symbol.column: character with column name of Gene Symbols in 'data'.
  #                if NULL, SYMBOLs are derived from 'Entrez.column' using the annotation package for 'org'.
  # Entrez.column: character with column name of ENTREZ IDs in 'data'.
  #                if NULL, ENTREZ IDs are derived from 'Symbol.column' using the annotation package for 'org'.
  #                Either 'Symbol.column' or 'Entrez.column' must be specified.
  # org: character with species name ("human", "mouse", "rat").
  
  
  ## Details
  # If either gene symbols or EntrezIDs are missing in 'data', this annotation is
  # added using the corresponding annotation package of the organism specified in 'org'.
  

  ## value:
  # dataframe with added annotation column.
  
  
  ## Author(s) 
  # Frank R?hle 
  


  # Load genomic annotation package (includes OrgDb, TxDb and GODb)
  annopkg <- switch(org, human = "Homo.sapiens", 
                    mouse = "Mus.musculus",
                    rat="Rattus.norvegicus")
  attach_package(pkg.bioc=annopkg, pkg.cran = "plyr")
  
  
  
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


  
  # Detaching libraries not needed any more
  # detach_package(c(annopkg))
  
  return(data)
  
} # end of function definition

