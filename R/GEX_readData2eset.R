
#' Create ExpressionSet object from expression data
#' 
#' ExpressionSet object assembled from expression, feature and phenotype data.
#' 
#' Expression data, phenotype data and optional feature data are used to make an ExpressionSet object.
#' Only numeric columns in \code{exprsData} are used. If no feature data is given, non-numeric data columns 
#' from \code{exprsData} are used as feature data (if any). If dataframes don't contain rownames as identifier, 
#' the column \code{ProbeID} is used as rownames in \code{exprsData} and \code{featureData} and column \code{sampleColumn} in \code{phenoData}.
#' Rownames of \code{exprsData} must match with rownames of \code{featureData} (if given). Column names of \code{exprsData}
#' must match with rownames of \code{phenoData}. Characters "/" and "-" in sample names or group names are 
#' replaced by ".". If either gene symbols or EntrezIDs are missing in feature data, this annotation is
#' added using the corresponding annotation package if available (species name is derived from \code{annotation}). 
#' The optional function given in \code{transform} is used to transform the expression data.
#' 
#' @param exprsData dataframe or character with filepath to expression data to be loaded.
#'            If no rownames are included, a column containing ProbeIDs is expected and used as feature names.
#' @param phenoData dataframe or character with filepath to sample sheet / phenotype data to be loaded.
#' @param featureData dataframe or character with filepath to feature data to be loaded.
#'              If NULL, non-numeric data columns from 'exprsData' are used as feature data (if any).
#' @param ProbeID character with name of the column in 'exprsData' and 'featureData' that contains 
#'          identifiers used to uniquely identify each probe. 
#' @param sampleColumn character with name of sample column in 'phenoData'.
#' @param groupColumn character with name of group column in 'phenoData'.
#' @param experimentData MIAME object with optional experiment description.
#' @param annotation character string specifying array annotation package (e.g. "Humanv4"). 
#' character() if not available.
#' @param transform optional function definition to transform expression data (e.g. log2 for logarithm of base 2).
#'            if NULL, no data transformation is performed.
#'            
#' @return ExpressionSet object
#' 
#' @author Frank Ruehle
#' 
#' @export readData2eset




## Usage 
readData2eset <- function(exprsData, phenoData, featureData=NULL,
                          ProbeID = "PROBE_ID",
                          sampleColumn = "Sample_Name", 
                          groupColumn  = "Sample_Group",  
                          experimentData = MIAME(),
                          annotation= character(),
                          transform = log2
                          ) {

  
  
  
  
  # load required packages. limma and Biobase are not detached afterwards
  pkg.bioc <- c("limma", "Biobase")
  pks2detach <- attach_package(pkg.bioc=pkg.bioc)
  
  
### Expression data 
  # read file if 'exprsData' is character string with file path 
  # column names as they are, not checked/modified for valid R names!
  if(is.character(exprsData) & length(exprsData)==1) {
    cat("\n\nReading expression data:", exprsData, "\n")
    exprsData <- read.table(exprsData, header=T, sep="\t", na.strings = c("", " ", "NA"), check.names =F) 
    }
  

  exprscolumns <- sapply(exprsData, is.numeric)  # boolean vector with TRUE for numeric data columns
  annodata <- exprsData[, !exprscolumns, drop=F] # non-numeric data columns
      
  colnames(exprsData) <- sub(".AVG_Signal", "", colnames(exprsData)) # optional remove ".AVG_Signal" from column names
  exprsData <- as.matrix(exprsData[,exprscolumns]) # make matrix with numeric data
  
  if (is.null(rownames(exprsData))) {rownames(exprsData) <- annodata[,ProbeID]} # add rownames if not existing
  exprsData <- exprsData[,sort(colnames(exprsData))] # sort columns (Samples)
  
  if(!is.null(transform)) {exprsData <- transform(exprsData)} # transform expression data



### Feature data  
  if(!is.null(featureData)) {
    if(is.character(featureData) & length(featureData)==1) {
      cat("\n\nReading feature data:", featureData, "\n")
      featureData <- read.table(featureData, header=T, sep="\t", na.strings = c("", " ", "NA"), check.names =F) 
      }
    
  } else { # if no separate feature file loaded, use annotation columns from 'exprsData'
    featureData <- annodata
  }
  
  if(.row_names_info(featureData, type=1L) <0) { # if just automatically generated rownames, use ProbeID as rownames
    rownames(featureData) <- featureData[,ProbeID]
  }
  
  # order features in featureData if necessary
  featureData <- featureData[match(rownames(featureData), rownames(exprsData)),]

  # annotate with Gene symbols or ENTREZIDs with annotation package.
  if(length(annotation)>0) {  # define organism from annotation (if given)
    organism <- NULL 
      if(grepl("rat", annotation, ignore.case=T))   {organism <- "rat"}
      if(grepl("mouse", annotation, ignore.case=T)) {organism <- "mouse"}
      if(grepl("human", annotation, ignore.case=T)) {organism <- "human"}
        
    Symbol.column.name <- grep("SYMBOL", colnames(featureData), ignore.case=T, value=T)[1] # guess column names from featureData
      if(is.na(Symbol.column.name)) {Symbol.column.name <- NULL} # Null if no grep result
    Entrez.column.name <- grep("ENTREZ", colnames(featureData), ignore.case=T, value=T)[1] # guess column names from featureData
      if(is.na(Entrez.column.name)) {Entrez.column.name <- NULL} # Null if no grep result
    
    if(!is.null(organism)) {
      featureData <- basicAnno(data=featureData, Symbol.column = Symbol.column.name, Entrez.column = Entrez.column.name, org=organism)
    }
  } 
  
  # check if feature names are congruent
  if(any(rownames(featureData) != rownames(exprsData))) {stop("\nFeature names of 'exprsData' and 'featureData' don't match\n")}
  
  
  
#### Pheno data
  if(is.character(phenoData) & length(phenoData)==1) {
    cat("\n\nReading phenotype data:", phenoData, "\n")
    phenoData <- read.table(phenoData, header=T, sep="\t", na.strings = c("", " ", "NA"), check.names =F) 
    }
  
  # if(.row_names_info(phenoData, type=1L) <0) { # if just automatically generated rownames, use sampleColumn as rownames
    rownames(phenoData) <- phenoData[,sampleColumn]
  # } # if condition not necessary
  
  phenoData <- phenoData[sort(rownames(phenoData)),] # sorting rows (samples)
  
  
  # check if feature names are congruent
  if(any(rownames(phenoData) != colnames(exprsData))) {stop("\nSample names of 'exprsData' and 'phenoData' don't match\n")}
  
  

#### make ExpressionSet 
  eset <- ExpressionSet(assayData = exprsData,
                            phenoData = AnnotatedDataFrame(phenoData), 
                            experimentData = experimentData,
                            annotation = annotation)
  
  if (length(featureData)>0) {
    fData(eset) <- featureData # add feature data to eset
    }
  
  # purify sample names and group names from potential operator characters.
  for(op in c("/", "-")) {
      sampleNames(eset) <- gsub(op, ".", sampleNames(eset))
      pData(eset)[,sampleColumn] <- gsub(op, ".",  pData(eset)[,sampleColumn])
      pData(eset)[,groupColumn] <- gsub(op, ".",  pData(eset)[,groupColumn])
  } 
  

  

return(eset)
} # end function definition


