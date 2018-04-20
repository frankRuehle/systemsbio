

#' Reading expression data from Illumina arrays 
#' 
#' Create an ExpressionSetIllumina object from Illumina GenomeStudio output files.
#' 
#' Unnormalized expression data from an Illumina GenomeStudio project is read into a \code{ExpressionSetIllumina} object.
#' If a covariate file is given, covariates are included into the phenotype data of the object. 
#' If a paired sample design was chosen, the variable indicating sample pairing must be included either in
#' the sample sheet file or in the covariates file. In default settings, the expression data is log2-transformed
#' but not normalized.
#' The feature data of the ExpressionSetIllumina object is supplemented by annotation data given in 
#' the respective Illumina annotation package (if available \code{REANNOTATED}-columns used). 
#' The respective chip-type is given in \code{exprchip}.
#' 
#' 
#' @param dataFile character with filepath to SampleProbeProfile from GenomeStudio (tab-delimited txt-file) 
#'           containing columns with \code{AVG_Signal}, \code{Detection Pval}, \code{BEAD_STDERR} und \code{Avg_NBEADS}.
#'           Additional columns are ignored except for \code{"PROBE_ID"}, \code{"SYMBOL"} and \code{"ENTREZ_GENE_ID"}.
#' @param qcFile character with filepath to \code{ControlProbeProfile} from GenomeStudio (tab-delimited txt-file)
#'         containing columns with \code{AVG_Signal} und \code{Detection Pval} 
#' @param sampleSheet character with filepath to Sample sheet from GenomeStudio project (csv-file). The header is expected in row 8.
#' @param ProbeID character string with name of the column in dataFile that contains identifiers used to uniquely identify each probe 
#' @param skip number of header lines to skip at the top of \code{dataFile}. 
#' @param controlID character string specifying the column in qcFile that contains the identifiers that 
#'            uniquely identify each control probe. If the \code{ControlGeneProfile} from GenomeStudio is used, 
#'            you will need to set \code{controlID="TargetID"} 
#' @param qc.skip number of header lines to skip at the top of qcFile 
#' @param qc.columns list defining the column headings in qcFile which correspond to the matrices stored in the 
#'             QCInfo slot of the final ExpressionSetIllumina object
#' @param sampleColumn Name of sample column in sample sheet
#' @param groupColumn Name of group column in sample sheet
#' @param exprchip character string specifying expression chip type (e.g. \code{"HumanHT-12 v4"}, \code{"MouseWG-6 v2"}, \code{"MouseRef-8 v2"}, \code{"RatRef-12 v1"}.
#' @param org character string specifying organism: \code{"human"}, \code{"rat"}, \code{"mouse"}.
#' @param covarfile character with filepath if a covariates file is given. NULL otherwise.
#' @param covarsampleID character with name of sample column in (optional) covariates file.   
#' @param method_norm character with normalisation method. Options are \code{"quantile"}, \code{"qspline"}, 
#' \code{"vsn"}, \code{"rankInvariant"}, \code{"median"} and \code{"none"}
#' @param transform character with data transformation method. Options are \code{"none"}, \code{"log2"}, \code{"neqc"}, \code{"rsn"} and \code{"vst"}.
#' @param fields2Add character vector with names of Illumina mappings to add to feature data (chip type dedicated in \code{exprchip}).
#' Remark: works only if feature data is \code{"PROBE_ID"}, but this identifier may be non-unique and may therefore throw an error when loading.
#'
#'
#' @return Annotated (and log-transformed) ExpressionSetIllumina object
#' 
#' @author Frank Ruehle
#' 
#' @export 


## Usage 
read_Illu_expr_array <- function(
      dataFile, 
      qcFile, 
      sampleSheet,
      ProbeID = "ProbeID", 
      skip = 0, 
      controlID= "ProbeID", 
      qc.skip = 0, 
      qc.columns = list(exprs = "AVG_Signal", Detection = "Detection Pval"),
      sampleColumn = "Sample_Name", 
      groupColumn  = "Sample_Group",  
      exprchip= NULL,
      org= NULL,
      covarfile = NULL, 
      covarsampleID = "ID",   
      method_norm = "none", 
      transform= "log2",  
      fields2Add= NULL 
      ) {


  
  # load required libraries
  pkg.cran <- NULL
  pkg.bioc <- c("beadarray", "Biobase") # pkg.bioc not detached afterwards
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  # determine name of Illumina array annotation package and install/load it
  if(!is.null(exprchip)) {
    if(is.null(org)) {stop("\nOrganism not specified!")}
    org <- tolower(org)
    if(!(org %in% c("human", "mouse", "rat"))) {stop("\norg needs to be either human, mouse or rat")}
    exprchip <- tolower(exprchip)
    exprchipBez     <- sub("-.*", "", exprchip)
    exprchipVersion <- sub(".*v", "v", exprchip)
    ArrayAnnotation.GEX <- switch(org, human = paste("Human", exprchipVersion, sep=""), 
                                  mouse = paste("Mouse", exprchipVersion, sep=""),
                                  rat = paste("Rat", exprchipVersion, sep=""))
    chipannopkg <- paste0("illumina", ArrayAnnotation.GEX, ".db")
    
    attach_package(pkg.bioc=chipannopkg)
    
    #cat("\nArray Annotation package:", chipannopkg)
    cat("\nSee keytypes(", chipannopkg, ") for available annotation data.")
  }

# Reading Bead-Summary Data
cat("\nRead bead data. Illumina annotation package:", chipannopkg, "\n")
  eset = readBeadSummaryData(dataFile = dataFile, qcFile = qcFile, sampleSheet = sampleSheet, sep="\t",
                             columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR",
                                            nObservations = "Avg_NBEADS", Detection="Detection Pval"),
                             ProbeID = ProbeID, skip = skip, 
                             controlID = controlID, qc.skip = qc.skip, qc.columns = qc.columns,
                             illuminaAnnotation = ArrayAnnotation.GEX,
                             annoCols = c("PROBE_ID", "SYMBOL", "ENTREZ_GENE_ID")) #

  # Comment: the following warning messages are observed frequently:
  # Warning messages:
  #   1: In readQC(file = qcFile, sep = qc.sep, skip = qc.skip, columns = qc.columns,  :
  #                  controlIDs non-unique: 4 repeated entries have been removed.
  #   2: In grep(colnames(exprs(BSData)), samples, fixed = TRUE) :
  #    argument 'pattern' has length > 1 and only the first element will be used
  #
  # @1: this is due to a property of the Illumina array which has 4 identically named QC probes.
  # @2: the first Sample ID from colnames(exprs(BSData)) is sufficient to identify the column of sample IDs in the sample sheet
  

# Reading covariates file (optional)
if (!is.null(covarfile)) {
  cat("\nAdding covariates from", covarfile, "\n")
  covar <- read.delim(file=covarfile, header=T, sep="\t", stringsAsFactors =T, 
                      na.strings = c("NA", "-9", "", " "), fileEncoding="latin1") 
  covar[,covarsampleID] <- sub("@.*$", "", covar[,covarsampleID])  # remove "@cohort"-suffix from SNPZone database if present
  covar[,covarsampleID] <- sub(" ", "", covar[,covarsampleID]) # remove whitespace in IDs if present
  
  # join() is used for merging because merge() does not preserve the roworder (even when sort=F)
  names(covar)[names(covar)==covarsampleID] <- sampleColumn # for join() bost dataframes must have the same column name
  pData(eset) <- plyr::join(pData(eset), covar, by=sampleColumn, type="left")
 
  # after joining, sampleNames are set to numbers, therefore renaming.
  sampleNames(eset) <- pData(eset)[,sampleColumn]
}

pData(eset)[,sampleColumn] <- factor(pData(eset)[,sampleColumn])
pData(eset)[,groupColumn]  <- factor(pData(eset)[,groupColumn])
if("Sentrix_ID" %in% names(pData(eset))) {pData(eset)$Sentrix_ID <- factor(pData(eset)$Sentrix_ID)}
if("Sentrix_Position" %in% names(pData(eset))) {pData(eset)$Sentrix_Position <- factor(pData(eset)$Sentrix_Position)}
#if (!is.null(matchvar)) {pData(eset)[,matchvar] <- factor(pData(eset)[,matchvar])}

# name row names by sample names, otherwise rownames from first column of pData
sampleNames(eset) <- pData(eset)[,sampleColumn]


## examples of data accession
# exprs(eset)[1:5,1:5]
# se.exprs(eset)[1:5,1:5]
# head(fData(eset))
# table(fData(eset)[,"Status"])
# pData(eset)
# channelNames(eset)
# sampleNames(eset)
# eset[1:10, 1:5] # Row indexing addresses regular probes only, i.e. no quality control probes  
# 


# transforming expression values (default is log2-transformation)
cat("\nnormalising (", method_norm, ") and/or transforming (", transform, ") expression data.\n")
eset <- normaliseIllumina(eset, method=method_norm, transform=transform)



###### add illumina annotation to feature data (if given)
cat("\nAdd Illumina annotation to feature data:", fields2Add, "\n")
if(!is.null(fields2Add)) {
  eset <- addFeatureData(eset, toAdd = fields2Add)
}


# If 'REANNOTATED' columns are available for SYMBOL and ENTREZID, these are used instead the original chip annotation.
if("SYMBOLREANNOTATED" %in% names(fData(eset))) {
  cat("\nSYMBOLREANNOTATED used for SYMBOL")
  if("SYMBOL" %in% names(fData(eset))) {names(fData(eset))[names(fData(eset))=="SYMBOL"] <- "SYMBOL_chipannot"}
  names(fData(eset))[names(fData(eset))=="SYMBOLREANNOTATED"] <- "SYMBOL"
  }
if("ENTREZREANNOTATED" %in% names(fData(eset))) {
  cat("\nENTREZREANNOTATED used for ENTREZID")
  if("ENTREZID" %in% names(fData(eset))) {names(fData(eset))[names(fData(eset))=="ENTREZID"] <- "ENTREZID_chipannot"}
  names(fData(eset))[names(fData(eset))=="ENTREZREANNOTATED"] <- "ENTREZID"
}



# remove redundant feature data columns
cols2remove <- c("Row.names")
cols2keep <- !(colnames(fData(eset)) %in% cols2remove)
fData(eset) <- fData(eset)[,cols2keep]
# "ProbeID" from QC probes (and regular probes)
# "PROBE_ID" from regular probes


# Detaching libraries not needed any more
detach_package(c(pkg.cran, chipannopkg))


return(eset)

}


