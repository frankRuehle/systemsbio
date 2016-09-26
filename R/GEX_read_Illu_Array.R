

## Description
# Reading expression data from Illumina arrays 

## Usage 
read_Illu_Array <- function(
      dataFile, 
      qcFile, 
      sampleSheet,
      ProbeID = "PROBE_ID", 
      skip = 0, 
      controlID= "ProbeID", 
      qc.skip = 0, 
      qc.columns = list(exprs = "AVG_Signal", Detection = "Detection Pval"),
      sampleColumn = "Sample_Name", 
      groupColumn  = "Sample_Group",  
      exprchip= "HumanHT-12 v4",
      org= "human",
      covarfile = NULL, 
      covarsampleID = "ID",   
      matchvar = NULL, 
      method_norm = "none", 
      transform= "log2",  
      fields2Add= NULL 
      ) {


   
  ## Arguments
  # dataFile: character with filepath to SampleProbeProfile from GenomeStudio (tab-delimited txt-file) 
  #           containing columns with AVG_Signal, Detection Pval, BEAD_STDERR und Avg_NBEADS
  # qcFile: character with filepath to ControlProbeProfile from GenomeStudio (tab-delimited txt-file)
  #         containing columns with AVG_Signal und Detection Pval 
  # sampleSheet: character with filepath to Sample sheet from GenomeStudio project (csv-file). The header is expected in row 8.
  # ProbeID: character string with name of the column in dataFile that contains identifiers used to uniquely identify each probe 
  # skip: number of header lines to skip at the top of dataFile. 
  # controlID: character string specifying the column in qcFile that contains the identifiers that 
  #            uniquely identify each control probe 
  # qc.skip: number of header lines to skip at the top of qcFile 
  # qc.columns: list defining the column headings in qcFile which correspond to the matrices stored in the 
  #             QCInfo slot of the final ExpressionSetIllumina object
  # sampleColumn: Name of sample column in sample sheet
  # groupColumn: Name of group column in sample sheet
  # exprchip: character string specifying expression chip type (e.g. "HumanHT-12 v?", "MouseWG-6 v?", "MouseRef-8 v?")
  # org: character string specifying organism: "human", "rat", "mouse".
  # covarfile: character with filepath if a covariates file is given. NULL otherwise.
  # covarsampleID: character with name of sample column in (optional) covariates file.   
  # matchvar: NULL if unpaired study design. In  paired study design, character with variable name 
  #           indicating pairing of samples. This column must be contained either in sampleSheet or covariates file. 
  # method_norm: character with normalisation method. Options are "quantile", "qspline", "vsn", "rankInvariant", "median" and "none"
  # transform: character with data transformation method. Options are "none", "log2", "neqc", "rsn" and "vst".
  # fields2Add: character vector with names of Illumina mappings to add to feature data (chip type dedicated in 'exprchip').
  
  
  ## Details
  # Unnormalized expression data from an Illumina GenomeStudio project is read into a ExpressionSetIllumina object.
  # If a covariate file is given, covariates are included in the phenotype data of the object. 
  # If a paired sample design was chosen, the variable indicating sample pairing must be included either in
  # the sample sheet file or in the covariates file. In default settings, the expression data is log2-transformed
  # but not normalized.
  # The feature data of the ExpressionSetIllumina object is supplemented by annotation data given in 
  # the respective Illumina annotation package (if available "REANNOTATED"-columns used). 
  # The respective chip-type is given in 'exprchip'.
  
   
  ## value:
  # Annotated (and log-transformed) ExpressionSetIllumina object
  
  
  ## Author(s) 
  # Frank R?hle 
  
  
  
  # load required libraries
  pkg.cran <- NULL
  pkg.bioc <- c("beadarray", "Biobase") # pkg.bioc not detached afterwards
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  # determine name of Illumina array annotation package and install/load it
  exprchip <- tolower(exprchip)
  exprchipBez     <- sub("-.*", "", exprchip)
  exprchipVersion <- sub(".*v", "v", exprchip)
  ArrayAnnotation.GEX <- switch(org, human = paste("Human", exprchipVersion, sep=""), 
                                mouse = paste("Mouse", exprchipVersion, sep=""))
  chipannopkg <- paste0("illumina", ArrayAnnotation.GEX, ".db")
  
  attach_package(pkg.bioc=chipannopkg)
  
  


# Reading Bead-Summary Data
cat("\nRead bead data. Illumina annotation package:", chipannopkg, "\n")
  eset = readBeadSummaryData(dataFile = dataFile, qcFile = qcFile, sampleSheet = sampleSheet, sep="\t",
                             columns = list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR",
                                            nObservations = "Avg_NBEADS", Detection="Detection Pval"),
                             ProbeID = ProbeID, skip = skip, 
                             controlID = controlID, qc.skip = qc.skip, qc.columns = qc.columns,
                             illuminaAnnotation = ArrayAnnotation.GEX,
                             annoCols = c("PROBE_ID","SYMBOL")) 

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
if (!is.null(matchvar)) {pData(eset)[,matchvar] <- factor(pData(eset)[,matchvar])}

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


