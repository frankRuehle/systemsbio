#' QTL Analysis for expression and/or methylation data.
#'
#' \code{matrixQTL} uses MatrixEQTL package for QTL analysis of an ExpressionSet or/and a MethylationSet object with SNP genotypes.
#' 
#' Expression and methylation data can be processed simulatanously within a single function call. Coordinates 
#' for genes and SNPs are either given in \code{SNPfile_transpAddCoded_name} or \code{genepos}, respectively,
#' or are downloaded from biomaRt. Model type and optional covariates can be selected for analysis. 
#' Analysis can be perfomed for either cis or trans QTLs or both simulatanously. 
#'
#'
#' @param eset ExpressionSet object.
#' @param mset MethylationSet object.
#' @param islandset data from summarized methylation Island (output from COHCAP.avg.by.island).
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
#' @param sampleGEX Character with column name of sample IDs in phenotype data of \code{eset}.
#' @param sampleMT Character with column name of sample IDs in phenotype data of \code{mset} or \code{islandset}.
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
#' @param genepos Character with path to gene position file for expression data. If NULL, genepos file is 
#' generated with biomaRt from eset probe IDs.
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


matrixQTL <- function(eset = NULL, 
                      mset = NULL, 
                      islandset = NULL, 
                      SNPfile_transpAddCoded_name, 
                      SNPfile_tfam = NULL, 
                      covariates_file_name = NULL, 
                      covarsampleID = "IID",
                      covar2adjust= NULL,   
                      projectfolder= "QTL",
                      projectname= NULL,  
                      sampleGEX = "Sample_Name",
                      sampleMT =  "Sample_Name",
                      
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

  
if(all(is.null(eset), is.null(mset), is.null(islandset))) {stop("\nNo quantitative data found.")}
  
## install/load required packages from CRAN and Bioconductor
pkg.bioc <- c("biomaRt")
pkg.cran <- c("MatrixEQTL", "GenABEL", "plyr")
attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)


## Check if required files are given
# quantitative data (expression or methylation)
if(is.null(eset) && is.null(mset) && is.null(islandset)) {
  stop("\n\nno expression or methylation file for QTL analysis\n")
  }

# SNP genotypes: transposed and additively coded ped file
if (is.null(SNPfile_transpAddCoded_name) || !file.exists(SNPfile_transpAddCoded_name)) {
  stop("\n\nno SNP-file for QTL analysis\n")
} 


# Covariates (optional)
if (is.null(covariates_file_name) || is.null(covar2adjust) || !file.exists(covariates_file_name)) {
  cat("\n\nno covariates file given or no covariates selected. QTL analysis performed without covariates.\n")
  covariates_file_name <- character()
    } 
  covarNameSuffix <- if (!is.null(covar2adjust)) {paste("_Covar", paste(covar2adjust, collapse="_"), sep="_")} else {""}

projectNameSuffix <- if (!is.null(projectname)) {paste0(projectname, "_")} else {""}

        

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
  






## store SNP-file with same samples as in eset or methylation files, respectively. 
if(!is.null(eset)) {
  if (!file.exists(file.path(projectfolder, "eQTL"))) {dir.create(file.path(projectfolder, "eQTL"))} # create subdirectory if not yet existing
  
  samplelist.GEX <- pData(eset)[,sampleGEX]  # Expression Samples
  samples.eQTL   <- colnames(SNPfile) %in% as.character(samplelist.GEX)
  samplenames.eQTL <- colnames(SNPfile)[samples.eQTL]
  
  cat("\n",length(samplenames.eQTL), "samples used for expression QTL analysis: \n", samplenames.eQTL, "\n")
  write.table(samplenames.eQTL, file.path(projectfolder, "eQTL", paste0(projectNameSuffix, "SampleList_eQTL.txt")), quote = F, sep = "\t",  row.names = F, col.names = F)
  SNP_file_name_GEX <- file.path(projectfolder, "eQTL", paste0(projectNameSuffix, "SNPfile_eQTL.traw") )
  cat("\nwrite SNP-file for expressionQTL to: ", SNP_file_name_GEX,"\n") 
  samples.eQTL[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
  write.table(SNPfile[,samples.eQTL], file=SNP_file_name_GEX, quote = F, sep = "\t",  row.names = F) 
  }

if(!is.null(mset)) {
  if (!file.exists(file.path(projectfolder, "mQTL"))) {dir.create(file.path(projectfolder, "mQTL"))} # create subdirectory if not yet existing
  
  samplelist.MT  <- pData(mset)[,sampleMT]  # Methylation Samples
  samples.mQTL <- colnames(SNPfile) %in% as.character(samplelist.MT)
  samplenames.mQTL <- colnames(SNPfile)[samples.mQTL]
  
  cat("\n",length(samplenames.mQTL), "samples used for methylation QTL analysis: \n", samplenames.mQTL, "\n")
  write.table(samplenames.mQTL, file.path(projectfolder, "mQTL", paste0(projectNameSuffix, "SampleList_mQTL.txt")), quote = F, sep = "\t",  row.names = F)
  SNP_file_name_MT <- file.path(projectfolder, "mQTL", paste0(projectNameSuffix,"SNPfile_mQTL.traw") )
  cat("\nwrite SNP-file for methylationQTL to: ", SNP_file_name_MT, "\n")
  samples.mQTL[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
  write.table(SNPfile[,samples.mQTL], quote = F, sep = "\t",  row.names = F, file= SNP_file_name_MT)
  }

if(!is.null(islandset)) {
  if (!file.exists(file.path(projectfolder, "mQTL_island"))) {dir.create(file.path(projectfolder, "mQTL_island"))} # create subdirectory if not yet existing
  
  annotationCols2remove <- c("island", "gene") # islandset contains annotation columns
  samplelist.MTisland  <- colnames(islandset)[!(colnames(islandset) %in% annotationCols2remove)] 
  samples.mQTLisland <- colnames(SNPfile) %in% as.character(samplelist.MTisland)
  samplenames.mQTLisland <- colnames(SNPfile)[samples.mQTLisland]
  
  cat("\n",length(samplenames.mQTLisland), "samples used for methylation QTL analysis with CpG islands: \n", samplenames.mQTLisland, "\n")
  write.table(samplenames.mQTLisland, file.path(projectfolder, "mQTL_island", paste0(projectNameSuffix, "SampleList_mQTL_cpgIslands.txt")), quote = F, sep = "\t",  row.names = F)
  SNP_file_name_MTisland <- file.path(projectfolder, "mQTL_island", paste0(projectNameSuffix, "SNPfile_mQTL_cpgIslands.traw") )
  cat("\nwrite SNP-file for methylationQTL to: ", SNP_file_name_MTisland, "\n")
  samples.mQTLisland[colnames(SNPfile)=="SNP"] <- TRUE # column "SNP" selected for saving together with sample columns
  write.table(SNPfile[,samples.mQTLisland], quote = F, sep = "\t",  row.names = F, file= SNP_file_name_MTisland)
}

rm(SNPfile) # SNPfile not needed any more

          
          






########### 3.) Expression/Methylation file and Genepos-file 
if(!is.null(eset)) {
  
ids <- featureNames(eset[which(fData(eset)[,"Status"] == "regular"),]) # feature names of "regular" probes
expression_file <- exprs(eset[ids, pData(eset)[,sampleGEX] %in% samplenames.eQTL])  # QC probes removed

  if (is.null(genepos)) {
  # Gene positions are downloaded using the supplied biomaRt object.
  # Ensemble filer is given in function call, e.g. Illumina probes.
  # Unfortunately, biomaRt does not find all Illumina probes, but this is not mandatory.
    
  cat("\n\nDownloading gene positions for expression file from biomaRt\n")
  geneposEnsemblData <- getBM(attributes = c(ensembl_filter, 'chromosome_name', 'start_position', 'end_position'), 
                              filters = ensembl_filter, values = rownames(expression_file), mart = genemart)
  # Duplicated entries and entries on chromosome patches are removed. 
  geneposEnsemblData_purif <- geneposEnsemblData[!duplicated(geneposEnsemblData[,1]),]
  geneposEnsemblData_purif <- geneposEnsemblData_purif[!grepl("[_HG]", geneposEnsemblData_purif$chromosome_name),] # catch PATCH chromosomes
  
  genepos <- data.frame(geneid=geneposEnsemblData_purif[,ensembl_filter], 
                       chr=paste0("chr", geneposEnsemblData_purif$chromosome_name), 
                       left=geneposEnsemblData_purif$start_position, right=geneposEnsemblData_purif$end_position)
  write.table(genepos, file=file.path(projectfolder, "eQTL", paste0(projectNameSuffix,"genepos_eQTL.txt")), quote = F, sep = "\t",  row.names = F)
  
  } else {
    genepos <- read.table(file=genepos, header=T, sep="\t") 
      }
}


## Methylationfile (analogous to expression data)  
if(!is.null(mset)) {

methylation_file <- getBeta(mset, type = "Illumina")
# gene symbols attached to rownames (in methylpos, too!):
rownames(methylation_file) <- paste(rownames(methylation_file), mcols(mset)$SYMBOL, sep="_")
methylation_file <- methylation_file[, colnames(methylation_file) %in% samplenames.mQTL] # restrict to samples from samplenames.mQTL

methylpos= data.frame(geneid=paste(mcols(mset)$Name, mcols(mset)$SYMBOL, sep="_"), chr= mcols(mset)$chr,    
                      left=mcols(mset)$pos, right=mcols(mset)$pos)

}

## Methylationfile with CpG-islands (from COHCAP.avg.by.island())  
if (!is.null(islandset))  {
# COHCAPdiff.samplewise nehmen --> filtered.cpgislands.allgroups (ungefiltert)
methylation_file_island <- islandset[,colnames(islandset) %in% samplenames.mQTLisland] # restrict to samples from samplenames.mQTLisland
methylation_file_island <- as.matrix(methylation_file_island)
rownames(methylation_file_island) <- paste(islandset[,"island"], islandset[,"gene"], sep="_")

islandPosWithoutCHR <- sub("^.*:", "", islandset[,"island"])   # format: "chr1:10000-15000"

methylpos.island <- data.frame(geneid=paste(islandset[,"island"], islandset[,"gene"], sep="_"), 
                             chr = sub(":.*", "", islandset[,"island"]),    
                             left = as.numeric(sub("-.*", "", islandPosWithoutCHR)), 
                             right = as.numeric(sub("^.*-", "", islandPosWithoutCHR)))

}



########## 4.) Covariates file:
if (length(covariates_file_name)>0) {
        cat("\nLoading covariates file: ", covariates_file_name, "\n")
        covar_file <- read.table(file=covariates_file_name, header=T)
        # covar_file <- covar_file[,-1]  # remove FID column
        covar_file[,covarsampleID] <- gsub("@[[:digit:]]*$", "", covar_file[,covarsampleID])  # remove cohort-suffix from SNPZone if necessary
        covar_file <- covar_file[,c(covarsampleID,covar2adjust)]
        names(covar_file)[names(covar_file)==covarsampleID] <- "id"
        
        # MatrixEQTL allows just figures in covar-file
        for (cov in covar2adjust) {
            if (!is.numeric(covar_file[,cov])) {
            cat("\nCovariate", cov, "is converted to numeric by factor level.\n")
            covar_file[,cov] <- as.numeric(factor(covar_file[,cov]))
                }
          }
        
        
        if(!is.null(eset)) {
          covar_file_eQTL <- covar_file[covar_file$id %in% as.character(samplenames.eQTL),] # nur die in "samplelist.GEX" aufgef?hrten Samples ausgew?hlt
          covar_file_eQTL <- t(covar_file_eQTL)
          covariates_file_name_GEX <- file.path(projectfolder, "eQTL", "eQTL.covar")
          cat("Writing covariates for eQTL samples to: ", covariates_file_name_GEX, "\n")
          write.table(covar_file_eQTL, file=covariates_file_name_GEX, quote = F, sep = "\t",  col.names = F)
          }

        if(!is.null(mset)) {
          covar_file_mQTL <- covar_file[covar_file$id %in% as.character(samplenames.mQTL),] # nur die in "samplelist.MT" aufgef?hrten Samples ausgew?hlt
          covar_file_mQTL <- t(covar_file_mQTL)
          covariates_file_name_MT <- file.path(projectfolder, "mQTL", "mQTL.covar")
          cat("Writing covariates for mQTL samples to: ", covariates_file_name_MT, "\n")
          write.table(covar_file_mQTL, file=covariates_file_name_MT, quote = F, sep = "\t",  col.names = F)
          }
        
        if(!is.null(islandset)) {
          covar_file_mQTL_island <- covar_file[covar_file$id %in% as.character(samplenames.mQTLisland),] # nur die in "samplenames.mQTLisland" aufgef?hrten Samples ausgew?hlt
          covar_file_mQTL_island <- t(covar_file_mQTL_island)
          covariates_file_name_MT_island <- file.path(projectfolder, "mQTL_island", "mQTL_island.covar")
          cat("Writing covariates for mQTL islands samples to: ", covariates_file_name_MT_island, "\n")
          write.table(covar_file_mQTL_island, file=covariates_file_name_MT_island, quote = F, sep = "\t",  col.names = F)
        }

      }







# Initialisation of qtl loop for the required QTL analysis (eQTL, mQTL and/or mQTL_island)
# Results are saved in different repositories.

helpE <- if (is.null(eset))       {NULL} else {"eQTL"}  
helpM <- if (is.null(mset))       {NULL} else {"mQTL"}  
helpI <- if (is.null(islandset))  {NULL} else {"mQTL_island"}  
helpvec <- c(helpE, helpM, helpI)

returnResQTL <- list()  # initialisation of list which will be returned

for (qtl in helpvec) {
  
  SNP_file_name   <- switch(qtl, "eQTL"=SNP_file_name_GEX, "mQTL"=SNP_file_name_MT, "mQTL_island"=SNP_file_name_MTisland)
  data_file       <- switch(qtl, "eQTL"=expression_file, "mQTL"=methylation_file, "mQTL_island"=methylation_file_island)
  feature.pos     <- switch(qtl, "eQTL"=genepos, "mQTL"=methylpos, "mQTL_island"=methylpos.island)
  if (length(covariates_file_name)==0) {covar_file_name <- character()} else {
      covar_file_name <- switch(qtl, "eQTL"=covariates_file_name_GEX, "mQTL"=covariates_file_name_MT, "mQTL_island"=covariates_file_name_MT_island)}
  
  pvOutputThreshold      <- switch(QTLtype, "both"=pvOutputThreshold_tra, "cis"=0, "all"=pvOutputThreshold_all)
  pvOutputThreshold.cis  <- switch(QTLtype, "both"=pvOutputThreshold_cis, "cis"=pvOutputThreshold_cis, "all"=0)
  
  # output filenames: (covariates are added to output filename, if applied)
  output_file_name  <- switch(QTLtype, "both" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_trans",covarNameSuffix,".txt")), 
                              "cis" = "", 
                              "all" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_all",covarNameSuffix,".txt")))
  
  output_file_name.cis <- switch(QTLtype, "both" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_cis",covarNameSuffix,".txt")), 
                                 "cis" = file.path(projectfolder, qtl, paste0(projectNameSuffix, qtl,"_result_cis",covarNameSuffix,".txt")), 
                                 "all" = "")
  
  
  #### preparation for matrixQTL
  
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
    # gene$LoadFile(expression_file_name)      # wenn Datei mit Expressionsdaten erst geladen werden soll
    gene$CreateFromMatrix(data_file)
    gene$ResliceCombined(sliceSize = 2000)     # n?tig weil obiges fileSliceSize-arg bei CreateFromMatrix offenbar nicht greift
    
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
  
  returnResQTL[[qtl]] <- resultQTL
  
  ## Results:
  cat("\nAnalysis done in: ", resultQTL$time.in.sec, " seconds", "\n");
  cat("\nResults stored in", file.path(projectfolder, qtl), "\n\n")
  ## Plot qq-plot of cis and trans p-values
  png(filename=file.path(projectfolder, qtl, paste0(projectNameSuffix, "qq_plot_", QTLtype,".png")), width = 210 , 
      height = 210, units = "mm", res=600)
      plot(resultQTL)
  dev.off()
  
  
} # end qtl-loop
  

 if(length(returnResQTL)==1) {returnResQTL <- returnResQTL[[1]]} 
 return(returnResQTL) 
  
} 


