#' QTL Analysis for expression and methylation data.
#'
#' \code{wrapMatrixEQTL} uses MatrixEQTL package for QTL analysis of expression or methylation data with SNP genotypes.
#' 
#' Expression or methylation data can be processed for QTL analysis. Coordinates 
#' for genes and SNPs are either given in \code{SNPfile_transpAddCoded_name} or \code{genepos}, respectively,
#' or are downloaded from biomaRt. Model type and optional covariates can be selected for analysis. 
#' Analysis can be perfomed for either cis or trans QTLs or both simulatanously. 
#'
#'
#' @param inputset Either ExpressionSet, SummarizedExperiment, DESeqDataSet, MethylSet or
#' data.frame containing summarized methylation Island data from COHCAP.avg.by.island-function (package COHCAP).
#' @param SNPfile_transpAddCoded_name character with path to genotype file (transposed and additive coded!).
#' @param SNPfile_tfam character with path to tfam file with sample information. 
#' If NULL, sample names are expected as header line of \code{SNPfile_transpAddCoded_name}. 
#' If given, \code{SNPfile_transpAddCoded_name} is read as no-header file and sample names are added 
#' from \code{SNPfile_tfam} as column names. Check in advance if you need \code{SNPfile_tfam} or not.
#' @param covariates_file_name character with path to covariates file.
#' @param covarsampleID Character with column name of sample IDs in covar file. 
#' @param covar2adjust Character vector with column names of covariates to adjust QTL analysis. Omitted if NULL.
#' @param projectfolder Character containing path to output folder (will be generated if not existing).
#' @param projectname Character used as suffix for output files.
#' @param sampleColumn Character with column name of sample IDs in input dataset.
#' @param useModel model to use (modelANOVA or modelLINEAR or modelLINEAR_CROSS). 
#' Set useModel = modelLINEAR to model the effect of the genotype as additive linear and test for its 
#' significance using t-statistic. Set useModel = modelANOVA to treat genotype as a categorical variables 
#' and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3.
#' Set otherwise like this: options(MatrixEQTL.ANOVA.categories=4). Set useModel = modelLINEAR_CROSS to add 
#' a new term to the model equal to the product of genotype and the last covariate; the significance of this 
#' term is then tested using t-statistic.
#' @param errorCovariance Numeric error covariance matrix. Use numeric() for homoskedastic independent errors.
#' @param QTLtype Character with "both" for calculating cis and trans QTL separately, 
#' "all" for no separation of cis and trans, "cis" for cis QTLs only.
#' @param cisDist Numeric maximal baisepair distance for cis gene-SNP pairs.
#' @param pvOutputThreshold_cis Numeric significance threshold p-value for cis QTL tests.
#' @param pvOutputThreshold_tra Numeric significance threshold p-value for trans QTL tests.
#' @param pvOutputThreshold_all Numeric significance threshold p-value for all QTL tests (cis and trans mixed).
#' @param genepos dataframe with 4 columns (geneid, chr, left, right) or character with path to gene position 
#' file. If \code{NULL}, genepos file is generated with biomaRt.
#' @param genemart biomaRt object to be used for updating gene positions.
#' @param ensembl_filter Character with filter name to search in genemart.
#' @param updateSNPpos Boolean. Shall SNP-positions be updated via biomaRt?
#' @param snpmart biomaRt object to be used for updating SNP positions.
#'  
#' @return list containing analysis parameter and QTL results. Intermediary results and plots 
#' are stored in \code{projectfolder} as side effects.
#' 
#' @author Frank Ruehle
#' 
#' @export 


wrapMatrixEQTL <- function(inputset, 
                      SNPfile_transpAddCoded_name, 
                      SNPfile_tfam = NULL, 
                      covariates_file_name = NULL, 
                      covarsampleID = "IID",
                      covar2adjust= NULL,   
                      projectfolder= "QTL",
                      projectname= NULL,  
                      sampleColumn = "Sample_Name",

                      useModel = modelLINEAR, 
                      errorCovariance = numeric(),  
                      
                      QTLtype = "both",  
                      cisDist = 1e6, 
                      pvOutputThreshold_cis = 1e-6, 
                      pvOutputThreshold_tra = 1e-6, 
                      pvOutputThreshold_all = 1e-6, 
                      genepos = NULL, 
                      genemart =  useMart("ENSEMBL_MART_ENSEMBL", host = "feb2014.archive.ensembl.org", dataset="hsapiens_gene_ensembl"), # GRCh37.p13, latest hg19 annotation
                      ensembl_filter = "illumina_humanht_12_v4",
                      updateSNPpos=FALSE,                             
                      snpmart =  useMart("ENSEMBL_MART_SNP", host = "feb2014.archive.ensembl.org", dataset="hsapiens_snp") # GRCh37.p13, latest hg19 annotation
                      ) {


  # Create output directory if not yet existing  
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T)}

  projectNameSuffix <- if (!is.null(projectname)) {paste0(projectname, "_")} else {""}
  covarNameSuffix <- ""
  
  
## install/load required packages from CRAN and Bioconductor
pkg.bioc <- c("biomaRt")
pkg.cran <- c("MatrixEQTL", "GenABEL", "plyr")
pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)



# SNP genotypes: transposed and additively coded ped file
if (is.null(SNPfile_transpAddCoded_name) || !file.exists(SNPfile_transpAddCoded_name)) {
  stop("\n\nno SNP-file for QTL analysis\n")
} 


# Covariates (optional)
if (!(is.null(covariates_file_name) || is.null(covar2adjust) || !file.exists(covariates_file_name))) {
  cov_exists =TRUE
  covarNameSuffix <- if (!is.null(covar2adjust)) {paste("_Covar", paste(covar2adjust, collapse="_"), sep="_")} else {""}

  cat("\nLoading covariates file: ", covariates_file_name, "\n")
  covar_file <- read.table(file=covariates_file_name, header=T)
  # covar_file <- covar_file[,-1]  # remove FID column
  covar_file[,covarsampleID] <- gsub("@[[:digit:]]*$", "", covar_file[,covarsampleID])  # remove cohort-suffix from SNPZone if necessary
  covar_file <- covar_file[,c(covarsampleID,covar2adjust)]
  names(covar_file)[names(covar_file)==covarsampleID] <- "id"
  
  # MatrixEQTL allows just numeric values in covar-file
  for (cov in covar2adjust) {
    if (!is.numeric(covar_file[,cov])) {
      cat("\nCovariate", cov, "is converted to numeric by factor level.\n")
      covar_file[,cov] <- as.numeric(factor(covar_file[,cov]))
      }
    }
    
  } else {
    cov_exists =FALSE
    cat("\n\nno covariates file given or no covariates selected. QTL analysis performed without covariates.\n")
    covariates_file_name <- character()
  } 
  
  
  
########## Input Files:
# SNP-file and covar-file must have the same samples as in expression-file or methylation-file!! 


### 1.) SNP-file: We need an additiv coded and transponed genotype file. Can be generated by
# GenABELs export.plink(geno.data, transpose = T, export012na = T) function or by PLINK 1.9:
# plink --file myfile --out myfile_transpAddCoded --recode A-transpose
cat("\n\nRead additively coded transposed SNP file: ", SNPfile_transpAddCoded_name, "\n")

  if(is.null(SNPfile_tfam)) {
    cat(paste("\nHeader expected in", SNPfile_transpAddCoded_name))
    SNPfile <- read.table(file=SNPfile_transpAddCoded_name, header=T, stringsAsFactors = FALSE, na.strings = c("NA", "", " ")) 
 
    
     } else {
    cat(paste("\nNo header expected in", SNPfile_transpAddCoded_name, ". Sample names loaded from:", SNPfile_tfam))
    SNPfile <- read.table(file=SNPfile_transpAddCoded_name, header=F, stringsAsFactors = FALSE, na.strings = c("NA", "", " ")) 
    tfam <- read.table(file=SNPfile_tfam, header=F, stringsAsFactors = FALSE, na.strings = c("NA", "", " "))
 
           if(nrow(tfam) == (ncol(SNPfile)-4) ) {
                  colnames(SNPfile) <- c("CHR", "SNP", "(C)M", "POS", paste(tfam[,1], tfam[,2], sep="_")) 
           } else {if(nrow(tfam) == (ncol(SNPfile)-3) ) {
             colnames(SNPfile) <- c("CHR", "SNP", "POS", paste(tfam[,1], tfam[,2], sep="_")) 
           } else {stop("SNPfile_transpAddCoded_name does not fit to SNPfile_tfam!")}
           } 
             
     }
    
  snpspos <- SNPfile[, c("CHR", "SNP", "POS")] # SNP location data is included in SNPfile
  
  names(SNPfile) <- gsub(".*_", "", names(SNPfile))   # FID have been added to IIDs by plink. They are removed.
  names(SNPfile) <- gsub("[@|\\.][[:digit:]]*$", "", names(SNPfile))  # Purify sample names from optional suffixes (@ or . followed by digits)
  col2remove <- names(SNPfile)  %in%  c("CHR", "X.C.M", "(C)M", "POS", "COUNTED", "ALT")  # needless columns
  SNPfile <- SNPfile[,!col2remove] # remove needless columns
  
  ###### 2.) Load SNP location file (map-file, coordinates optionally updated with biomaRt)
  #cat("\n\nRead SNPs location file: ", SNPfile_tfam, "\n")
  #snpspos <- read.table(SNPfile_tfam, header = F, stringsAsFactors = FALSE, na.strings = c("NA", "", " "))
  #if(ncol(snpspos) >= 4) {snpspos <- snpspos[,-3]}  # if map-file contains linkage column
  snpspos[,2] <- ifelse(snpspos[,2]==".", paste0("chr", snpspos[,1], ":", snpspos[,3]), snpspos[,2]) # if no rsid existing
  snpspos <- data.frame(snpid=snpspos[,2], chr=paste0("chr", snpspos[,1]), pos=snpspos[,3])
  
  if (updateSNPpos) {
    # Updating SNP positions using the supplied Biomart object. Attribute and Filter are the same for human and mouse. 
    cat("\nUpdating SNP positions with biomaRt: \n")
    snppos_martdata <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'),  
                             filters = 'snp_filter', values = snpspos$snpid, mart = snpmart) 
    dupSnps <- duplicated(snppos_martdata[,1]) # Remove duplikates from Biomart query
    snppos_martdata_nodup <- snppos_martdata[!dupSnps,]
    patchchrSnps <- grepl("CHR_", snppos_martdata_nodup[,2]) # Remove Patch-Chromosome entries from Biomart query
    snppos_martdata_nodup <- snppos_martdata_nodup[!patchchrSnps,]
    # mergen Biomart-Pos with snps pos(map-file). SNP count and SNP order in snpspos must not be changed!
    names(snppos_martdata_nodup)[names(snppos_martdata_nodup)=="refsnp_id"] <- "snpid"
    snpspos <- plyr::join(snpspos, snppos_martdata_nodup, by="snpid", type="left") 
    snpspos <- data.frame(snpid=snpspos$snpid, chr=paste0("chr", snpspos$chr_name), pos=snpspos$chrom_start)
  }
  
  snpspos$pos[is.na(snpspos$pos)] <- 0
  levels(snpspos$chr) <- c(levels(snpspos$chr), 0)
  snpspos$chr[is.na(snpspos$chr)] <- 0
  write.table(snpspos, file=file.path(projectfolder, paste0(projectNameSuffix,"snpspos.txt")), quote = F, sep = "\t",  row.names = F)
  



## 3a: use expressionSet and store SNP-file with relevant samples
if(class(inputset)== "ExpressionSet") {
  if (!file.exists(file.path(projectfolder, "eQTL"))) {dir.create(file.path(projectfolder, "eQTL"))} # create subdirectory if not yet existing
  qtl=c("eQTL") # filename prefix
  
  samplelist.GEX <- pData(inputset)[,sampleColumn]  # Expression Samples
  samples.eQTL   <- colnames(SNPfile) %in% as.character(samplelist.GEX)
  samplenames.eQTL <- colnames(SNPfile)[samples.eQTL]
  
  cat("\n",length(samplenames.eQTL), "samples used for expression QTL analysis: \n", samplenames.eQTL, "\n")
  write.table(samplenames.eQTL, file.path(projectfolder, "eQTL", paste0(projectNameSuffix, "SampleList_eQTL.txt")), quote = F, sep = "\t",  row.names = F, col.names = F)
  SNP_file_name <- file.path(projectfolder, "eQTL", paste0(projectNameSuffix, "SNPfile_eQTL.traw") )
  cat("\nwrite SNP-file for expressionQTL to: ", SNP_file_name,"\n") 
  samples.eQTL[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
  write.table(SNPfile[,samples.eQTL], file=SNP_file_name, quote = F, sep = "\t",  row.names = F) 
  
  ids <- featureNames(inputset[which(fData(inputset)[,"Status"] == "regular"),]) # feature names of "regular" probes
  data_file <- exprs(inputset[ids, pData(inputset)[,sampleColumn] %in% samplenames.eQTL])  # QC probes removed
  
      if (is.null(genepos)) {
        # Gene positions are downloaded using the supplied biomaRt object.
        # Ensemble filer is given in function call, e.g. Illumina probes.
        # Unfortunately, biomaRt does not find all Illumina probes, but this is not mandatory.
        
        cat("\n\nDownloading gene positions for expression file from biomaRt\n")
        geneposEnsemblData <- getBM(attributes = c(ensembl_filter, 'chromosome_name', 'start_position', 'end_position'), 
                                    filters = ensembl_filter, values = rownames(data_file), mart = genemart)
        colnames(geneposEnsemblData) <- c("filtered_ID", "chromosome_name", "start_position", "end_position") # otherwise colnames = biomart description field
        # Duplicated entries and entries on chromosome patches are removed. 
        geneposEnsemblData_purif <- geneposEnsemblData[!duplicated(geneposEnsemblData[,1]),] 
        geneposEnsemblData_purif <- geneposEnsemblData_purif[!grepl("[_HG]", geneposEnsemblData_purif$chromosome_name),] # catch PATCH chromosomes
        
        genepos <- data.frame(geneid=geneposEnsemblData_purif[,"filtered_ID"], 
                              chr=paste0("chr", geneposEnsemblData_purif$chromosome_name), 
                              left=geneposEnsemblData_purif$start_position, right=geneposEnsemblData_purif$end_position)
        write.table(genepos, file=file.path(projectfolder, "eQTL", paste0(projectNameSuffix,"genepos_eQTL.txt")), quote = F, sep = "\t",  row.names = F)
        
          } else {
            
            if(is.character(genepos)) {
              genepos <- read.table(file=genepos, header=T, sep="\t") 
              }           
            }
  
  feature.pos  <- genepos
                            
  if(cov_exists) { # covar file
    covar_file_eQTL <- covar_file[covar_file$id %in% as.character(samplenames.eQTL),] # select samples given in "samplelist.GEX" only
    covar_file_eQTL <- t(covar_file_eQTL)
    covar_file_name <- file.path(projectfolder, "eQTL", "eQTL.covar")
    cat("Writing covariates for eQTL samples to: ", covar_file_name, "\n")
    write.table(covar_file_eQTL, file=covar_file_name, quote = F, sep = "\t",  col.names = F)
  }
  
 }

 
## 3b: use SummarizedExperiment or DESeqDataSet and store SNP-file with relevant samples
  if(class(inputset) %in% c("SummarizedExperiment", "DESeqDataSet")) {
    if (!file.exists(file.path(projectfolder, "eQTL"))) {dir.create(file.path(projectfolder, "eQTL"))} # create subdirectory if not yet existing
    qtl=c("eQTL") # filename prefix
    
    samplelist.GEX <- colData(inputset)[,sampleColumn]  # Expression Samples
    samples.eQTL   <- colnames(SNPfile) %in% as.character(samplelist.GEX)
    samplenames.eQTL <- colnames(SNPfile)[samples.eQTL]
    
    cat("\n",length(samplenames.eQTL), "samples used for expression QTL analysis: \n", samplenames.eQTL, "\n")
    write.table(samplenames.eQTL, file.path(projectfolder, "eQTL", paste0(projectNameSuffix, "SampleList_eQTL.txt")), quote = F, sep = "\t",  row.names = F, col.names = F)
    SNP_file_name <- file.path(projectfolder, "eQTL", paste0(projectNameSuffix, "SNPfile_eQTL.traw") )
    cat("\nwrite SNP-file for expressionQTL to: ", SNP_file_name,"\n") 
    samples.eQTL[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
    write.table(SNPfile[,samples.eQTL], file=SNP_file_name, quote = F, sep = "\t",  row.names = F) 
    
    ids <- rownames(inputset)
    data_file <- assay(inputset[ids, colData(inputset)[,sampleColumn] %in% samplenames.eQTL])  # QC probes removed
    
    if (is.null(genepos)) {
      # Gene positions are downloaded using the supplied biomaRt object.
      # Ensemble filer is given in function call, e.g. Illumina probes.
      # Unfortunately, biomaRt does not find all Illumina probes, but this is not mandatory.
      
      cat("\n\nDownloading gene positions for expression file from biomaRt\n")
      geneposEnsemblData <- getBM(attributes = c(ensembl_filter, 'chromosome_name', 'start_position', 'end_position'), 
                                  filters = ensembl_filter, values = rownames(data_file), mart = genemart)
      colnames(geneposEnsemblData) <- c("filtered_ID", "chromosome_name", "start_position", "end_position") # otherwise colnames = biomart description field
      # Duplicated entries and entries on chromosome patches are removed. 
      geneposEnsemblData_purif <- geneposEnsemblData[!duplicated(geneposEnsemblData[,1]),]
      geneposEnsemblData_purif <- geneposEnsemblData_purif[!grepl("[_HG]", geneposEnsemblData_purif$chromosome_name),] # catch PATCH chromosomes
      
      genepos <- data.frame(geneid=geneposEnsemblData_purif[,"filtered_ID"], 
                            chr=paste0("chr", geneposEnsemblData_purif$chromosome_name), 
                            left=geneposEnsemblData_purif$start_position, right=geneposEnsemblData_purif$end_position)
      write.table(genepos, file=file.path(projectfolder, "eQTL", paste0(projectNameSuffix,"genepos_eQTL.txt")), quote = F, sep = "\t",  row.names = F)
      
    } else {
      
      if(is.character(genepos)) {
      genepos <- read.table(file=genepos, header=T, sep="\t") 
      } 
    }
    
    feature.pos  <- genepos # genepos given as dataframe
    
    if(cov_exists) { # covar file
      covar_file_eQTL <- covar_file[covar_file$id %in% as.character(samplenames.eQTL),] # select samples given in "samplelist.GEX" only
      covar_file_eQTL <- t(covar_file_eQTL)
      covar_file_name <- file.path(projectfolder, "eQTL", "eQTL.covar")
      cat("Writing covariates for eQTL samples to: ", covar_file_name, "\n")
      write.table(covar_file_eQTL, file=covar_file_name, quote = F, sep = "\t",  col.names = F)
    }
 }
  
  

## 3c: use Methylset and store SNP-file with relevant samples
if(class(inputset) %in% c("MethylSet", "RatioSet", "GenomicMethylSet", "GenomicRatioSet")) {
  if (!file.exists(file.path(projectfolder, "mQTL"))) {dir.create(file.path(projectfolder, "mQTL"))} # create subdirectory if not yet existing
  qtl=c("mQTL") # filename prefix
  
  samplelist.MT  <- pData(inputset)[,sampleColumn]  # Methylation Samples
  samples.mQTL <- colnames(SNPfile) %in% as.character(samplelist.MT)
  samplenames.mQTL <- colnames(SNPfile)[samples.mQTL]
  
  cat("\n",length(samplenames.mQTL), "samples used for methylation QTL analysis: \n", samplenames.mQTL, "\n")
  write.table(samplenames.mQTL, file.path(projectfolder, "mQTL", paste0(projectNameSuffix, "SampleList_mQTL.txt")), quote = F, sep = "\t",  row.names = F)
  SNP_file_name <- file.path(projectfolder, "mQTL", paste0(projectNameSuffix,"SNPfile_mQTL.traw") )
  cat("\nwrite SNP-file for methylationQTL to: ", SNP_file_name, "\n")
  samples.mQTL[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
  write.table(SNPfile[,samples.mQTL], quote = F, sep = "\t",  row.names = F, file= SNP_file_name)
 
    data_file <- getBeta(inputset, type = "Illumina")
    # gene symbols attached to rownames (in methylpos, too!):
    rownames(data_file) <- paste(rownames(data_file), mcols(inputset)$SYMBOL, sep="_")
    data_file <- data_file[, colnames(data_file) %in% samplenames.mQTL] # restrict to samples from samplenames.mQTL
    
    methylpos= data.frame(geneid=paste(mcols(inputset)$Name, mcols(inputset)$SYMBOL, sep="_"), chr= mcols(inputset)$chr,    
                          left=mcols(inputset)$pos, right=mcols(inputset)$pos)
    
    feature.pos <- methylpos
                              
  if(cov_exists) { # covar file
    covar_file_mQTL <- covar_file[covar_file$id %in% as.character(samplenames.mQTL),] # select samples given in "samplelist.MT" only
    covar_file_mQTL <- t(covar_file_mQTL)
    covar_file_name <- file.path(projectfolder, "mQTL", "mQTL.covar")
    cat("Writing covariates for mQTL samples to: ", covar_file_name, "\n")
    write.table(covar_file_mQTL, file=covar_file_name, quote = F, sep = "\t",  col.names = F)
  }
}


## 3d: use Methylset and store SNP-file with relevant samples
# Data frame from COHCAP.avg.by.island() of average beta (or percentage methylation) values across differentially 
# methylated sites within a differentially methylated CpG island.
  if(class(inputset)== "data.frame") { 
  if (!file.exists(file.path(projectfolder, "mQTL_island"))) {dir.create(file.path(projectfolder, "mQTL_island"))} # create subdirectory if not yet existing
  qtl=c("mQTL_island") # filename prefix
  
  annotationCols2remove <- c("island", "gene") # inputset contains annotation columns
  samplelist.MTisland  <- colnames(inputset)[!(colnames(inputset) %in% annotationCols2remove)] 
  samples.mQTLisland <- colnames(SNPfile) %in% as.character(samplelist.MTisland)
  samplenames.mQTLisland <- colnames(SNPfile)[samples.mQTLisland]
  
  cat("\n",length(samplenames.mQTLisland), "samples used for methylation QTL analysis with CpG islands: \n", samplenames.mQTLisland, "\n")
  write.table(samplenames.mQTLisland, file.path(projectfolder, "mQTL_island", paste0(projectNameSuffix, "SampleList_mQTL_cpgIslands.txt")), quote = F, sep = "\t",  row.names = F)
  SNP_file_name <- file.path(projectfolder, "mQTL_island", paste0(projectNameSuffix, "SNPfile_mQTL_cpgIslands.traw") )
  cat("\nwrite SNP-file for methylationQTL to: ", SNP_file_name, "\n")
  samples.mQTLisland[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
  write.table(SNPfile[,samples.mQTLisland], quote = F, sep = "\t",  row.names = F, file= SNP_file_name)

  ## Methylationfile with CpG-islands (from COHCAP.avg.by.island())  
    # COHCAPdiff.samplewise nehmen --> filtered.cpgislands.allgroups (ungefiltert)
    data_file <- inputset[,colnames(inputset) %in% samplenames.mQTLisland] # restrict to samples from samplenames.mQTLisland
    data_file <- as.matrix(data_file)
    rownames(data_file) <- paste(inputset[,"island"], inputset[,"gene"], sep="_")
    
    islandPosWithoutCHR <- sub("^.*:", "", inputset[,"island"])   # format: "chr1:10000-15000"
    
    methylpos.island <- data.frame(geneid=paste(inputset[,"island"], inputset[,"gene"], sep="_"), 
                                   chr = sub(":.*", "", inputset[,"island"]),    
                                   left = as.numeric(sub("-.*", "", islandPosWithoutCHR)), 
                                   right = as.numeric(sub("^.*-", "", islandPosWithoutCHR)))
    feature.pos <- methylpos.island
                              
  if(cov_exists) {
    covar_file_mQTL_island <- covar_file[covar_file$id %in% as.character(samplenames.mQTLisland),] # select samples given in "samplelist.MTisland" only
    covar_file_mQTL_island <- t(covar_file_mQTL_island)
    covar_file_name <- file.path(projectfolder, "mQTL_island", "mQTL_island.covar")
    cat("Writing covariates for mQTL islands samples to: ", covar_file_name, "\n")
    write.table(covar_file_mQTL_island, file=covar_file_name, quote = F, sep = "\t",  col.names = F)
  }
}


###########  
rm(SNPfile) # SNPfile not needed any more

          
          



## Start QTL analysis

  pvOutputThreshold      <- switch(QTLtype, "both"=pvOutputThreshold_tra, "cis"=0, "all"=pvOutputThreshold_all)
  pvOutputThreshold.cis  <- switch(QTLtype, "both"=pvOutputThreshold_cis, "cis"=pvOutputThreshold_cis, "all"=0)
  
  # output filenames: (covariates are added to output filename, if applied)
  output_file_name  <- switch(QTLtype, "both" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_trans",covarNameSuffix,".txt")), 
                              "cis" = "", 
                              "all" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_all",covarNameSuffix,".txt")))
  
  output_file_name.cis <- switch(QTLtype, "both" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_cis",covarNameSuffix,".txt")), 
                                 "cis" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_cis",covarNameSuffix,".txt")), 
                                 "all" = "")
  
  
  #### preparation for MatrixEQTL
  
    cat("\n\nReading SNP-file for", qtl, "Analysis\n")
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"      # the TAB character
    snps$fileOmitCharacters = "NA" # denote missing values;
    snps$fileSkipRows = 1          # one row of column labels
    snps$fileSkipColumns = 1       # number of columns of row labels 
    snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
    snps$LoadFile(SNP_file_name)
    
    cat("\nReading gene quantitative data file for", qtl, "Analysis\n")
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"      # the TAB character
    gene$fileOmitCharacters = "NA" # denote missing values;
    gene$fileSkipRows = 1          # one row of column labels
    gene$fileSkipColumns = 0       # number of columns of row labels
    gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
    # gene$LoadFile(expression_file_name)      # if file with expression data to be loaded first. 
    gene$CreateFromMatrix(data_file)
    gene$ResliceCombined(sliceSize = 2000)     
    
    cat("\nReading covariates for", qtl, "Analysis\n")
    cvrt = SlicedData$new()
    cvrt$fileDelimiter = "\t"      # the TAB character
    cvrt$fileOmitCharacters = "NA" # denote missing values;
    cvrt$fileSkipRows = 1          # one row of column labels
    cvrt$fileSkipColumns = 1       # one column of row labels
    if(length(covariates_file_name)>0) {cvrt$LoadFile(covar_file_name)}  # if no covariates available covariates_file_name=character()
  

  cat("\n\nRun", qtl, "analysis\n")
  
  resultQTL = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name.cis,
    pvOutputThreshold.cis = pvOutputThreshold.cis,
    snpspos = snpspos,
    genepos = feature.pos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",  # sonst TRUE for histogram
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  

  ## Results:
  cat("\nAnalysis done in: ", resultQTL$time.in.sec, " seconds", "\n");
  cat("\nResults stored in", file.path(projectfolder, qtl), "\n\n")
  ## Plot qq-plot of cis and trans p-values
  png(filename=file.path(projectfolder, qtl, paste0(projectNameSuffix, "qq_plot_", QTLtype,".png")), width = 210 , 
      height = 210, units = "mm", res=600)
      plot(resultQTL)
  dev.off()
  
  
  # Detaching libraries not needed any more
  detach_package(unique(pks2detach))
  
 return(resultQTL) 
  
} 


