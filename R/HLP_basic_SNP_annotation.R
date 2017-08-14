#' SNP annotation using biomaRt and/or manufacturer data
#'
#' \code{basic_SNP_annotation} adds annotation data to  SNP IDs using biomaRt and/or a 
#' manufacturer's annotation file while keeping the order of input SNPs.
#'
#' This function uses the SNP ID column from a given dataframe as input for adding annotation data.
#' All annotation data is added in addiotional columns and does not change the order of input
#' SNP IDs. Since biomaRt queries of large datasets (e.g. from SNP arrays) are prone to servive
#' malfunction, \code{basic_SNP_annotation} devides the data in chunks of feasable size given in
#' \code{max.SNPs.per.biomaRt.call}.
#' 
#' Before biomaRt data are merged to input data, data columns containing multiple entries per entries
#' (given in \code{biomaRt.attributes.summarized}) are collapsed seperated by a semicolon. Data columns 
#' given in \code{biomaRt.attributes.groupColumns} are considered as grouping variables.
#' 
#' If a \code{annotationFile} is specified, all included data is merged to the input dataframe.
#' The \code{annotationFile} may be supplied directly as dataframe or as character containing a 
#' file path. In latter case, the file is automatically loaded.
#'
#' @param data dataframe containing SNP IDs.
#' @param max.SNPs.per.biomaRt.call numeric. Number of SNP IDs to be querierd in biomaRt at once.
#' @param snpmaRt biomaRt object to be used for annotation. If NULL, biomaRt annotation is skipped.
#' @param data.SNP.columnName character with column name of SNP IDs in \code{data} or "row.names".
#' @param biomaRt.SNP.columnName character with attribute name for SNP IDs of the biomaRt object.
#' @param biomaRt.filter character with filter name to be used in biomaRt query.
#' @param biomaRt.attributes.groupColumns character vector with attribute names to be queried in biomaRt. 
#' @param biomaRt.attributes.summarized character vector with further attribute names, which will be 
#'         summarized according to the attributes in \code{biomaRt.attributes.groupColumns} (separated by ";").
#' @param annotationFile dataframe or character with path to dataframe containing annotation data by the 
#'         assay manufacturer. If NULL, annotation is skipped.
#' @param lines2skip.start Numeric with number of rows to skip when loading \code{annotationFile} or 
#'                         regular expression for character string to identify corresponding row number 
#'                         to be skipped, e.g. \code{[Assay]} in Illumina annotation files.
#' @param lines2skip.end Numeric with number of rows to read when loading \code{annotationFile} or 
#'                        regular expression for character string to identify corresponding row number to be read,
#'                        e.g. \code{[Controls]} in Illumina annotation files as start of annotation of control probes.
#'                        All rows from that number on (incl. \code{lines2skip.end}) are skipped. Negative and other 
#'                        invalid values are ignored.
#' @param annofile.SNP.columnName character with column name of SNP IDs in \code{annotationFile}.
#' @param annofile.columns Optional character vector with column names of \code{annotationFile} to be included.
#' If NULL, all columns of \code{annotationFile} are merged
#'
#' @return input dataframe annotated with biomaRt and/or manufacturer data in additional columns 
#' (starting with "SNPMart_" or "Annofile_", respectively). Order of entries within the dataframe 
#' remains unchanged.
#' 
#' @author Frank Ruehle
#' 
#' @export 


basic_SNP_annotation <- function(data,
                                 max.SNPs.per.biomaRt.call = 10000,
                                 data.SNP.columnName = "SNP",
                                 snpmaRt = useMart("ENSEMBL_MART_SNP", host = "feb2014.archive.ensembl.org", dataset="hsapiens_snp"), # GRCh37.p13, latest hg19 annotation
                                 biomaRt.SNP.columnName = "refsnp_id", 
                                 biomaRt.filter = "snp_filter",
                                 biomaRt.attributes.groupColumns  = c("refsnp_id", "chr_name", "chrom_start"),
                                 biomaRt.attributes.summarized = c("ensembl_gene_stable_id", "ensembl_type"),
                                 annotationFile = NULL,
                                 lines2skip.start = 0,
                                 lines2skip.end = -1,
                                 annofile.SNP.columnName ="Name",
                                 annofile.columns = NULL 
                                 ) {
  
  # biomaRt.SNP.columnName has to be included in biomaRt attributes
  if(!biomaRt.SNP.columnName %in% c(biomaRt.attributes.groupColumns, biomaRt.attributes.summarized)) {
    cat("\n", biomaRt.SNP.columnName, "has been added to biomaRt.attributes.groupColumns\n")
    biomaRt.attributes.groupColumns <- c(biomaRt.SNP.columnName, biomaRt.attributes.groupColumns)}
  
  ## install/load required packages from CRAN and Bioconductor
  pkg.bioc <- c("biomaRt")
  pkg.cran <- c("plyr")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  # read data if necessary
  if(is.character(data) && length(data)==1) {
    cat("\n\nReading data:", data, "\n")
    data <- read.table(data, header=T, sep="\t", na.strings = c("", " ", "NA"), check.names =T, stringsAsFactors = F) 
  }
   
  # add temporary index column
  data$temp.index <- 1:nrow(data)
    
  if(!is.null(snpmaRt)) {
    # define data subsets for biomaRt query
   if (nrow(data) <= max.SNPs.per.biomaRt.call) {
      start <- 1
      stop <- nrow(data)
        } else {
            start <- seq(from = 1, to = nrow(data), by = max.SNPs.per.biomaRt.call)
            stop <- c(seq(from = max.SNPs.per.biomaRt.call, to = nrow(data), by = max.SNPs.per.biomaRt.call),  nrow(data))
   }
 
  cat("\nAnnotate SNP-IDs with data from biomaRt\n")
  snp.anno <- data.frame()
  for(i in 1:length(start)){
    cat("biomaRt query", i, "of", length(start), "\n")
    snp.temp <- getBM(attributes = c(biomaRt.attributes.groupColumns, biomaRt.attributes.summarized), 
                      filters = biomaRt.filter, 
                      values = if(data.SNP.columnName == "row.names") {rownames(data[start[i]:stop[i],])
                                } else {data[start[i]:stop[i], data.SNP.columnName]}, 
                      mart = snpmaRt)
    snp.anno <- rbind(snp.anno, snp.temp)
  }
  
  snp.anno[snp.anno==""] <- NA 
  snp.anno[snp.anno==" "] <- NA 
  snp.anno[snp.anno=="NA"] <- NA 
  
    # remove entries with Patch chromosomes and set NA to 0
  if("chr_name" %in% names(snp.anno)) {
    snp.anno <- snp.anno[!grepl("H", snp.anno$chr_name),]
    } 
  
  # remove entries with missing rs_ids
  snp.anno <- snp.anno[!is.na(snp.anno[,biomaRt.SNP.columnName]),]
  snp.anno <- snp.anno[!duplicated(snp.anno[,biomaRt.SNP.columnName]),] # remove duplicate entries from biomart query (only first entry used)
  

  
  # format biomaRt results by collapsing columns given in biomaRt.attributes.summarized
  if(!is.null(biomaRt.attributes.summarized)) { 
    snp.anno <- ddply(snp.anno, biomaRt.attributes.groupColumns, colwise(paste, collapse=";"))
  }
   
  names(snp.anno) <- paste0("SNPMart_", names(snp.anno))
   
  # merge with data
  data <- merge(data, snp.anno, by.x= data.SNP.columnName, by.y=paste0("SNPMart_", biomaRt.SNP.columnName), all.x=T)
  
  # rename column rows if overwritten and remove temporara row.name columns
  if(data.SNP.columnName == "row.names") {
    rownames(data) <- data$Row.names
    data <- data[, !(names(data) %in% c("Row.names"))]
  }
} # end if snpmaRt
  
  if(!is.null(annotationFile)) {
    
    # read annotation file if necessary
    if(is.character(annotationFile) && length(annotationFile)==1) {
      cat("\nReading Annotation File:", annotationFile, "\n")
      separator <- if(grepl("\\.csv", annotationFile)) {","} else {"\t"}
   
      if(is.character(lines2skip.start)) {
        lines2skip.start <- grep(lines2skip.start, readLines(annotationFile), ignore.case = T)
        if (length(lines2skip.start==0)) {stop("\nCharacter in lines2skip.start not found!\n")}
        if (length(lines2skip.start>=2)) {stop(paste("\nCharacter in lines2skip.start found multiple times:", print(lines2skip.start), "\n"))}
      }
      
      if(is.character(lines2skip.end)) {
        lines2skip.end <- grep(lines2skip.end, readLines(annotationFile), ignore.case = T)
        if (length(lines2skip.end==0)) {stop("\nCharacter in lines2skip.end not found!\n")}
        if (length(lines2skip.end>=2)) {stop(paste("\nCharacter in lines2skip.end found multiple times:", print(lines2skip.end), "\n"))}
        }
      
    annotationFile <- read.table(annotationFile, header=T, sep= separator, na.strings = c("", " ", "NA"), check.names =F, 
                                   stringsAsFactors = F, skip = lines2skip.start, nrow= lines2skip.end-2-lines2skip.start) 
    }
    
    cat("\nAnnotate SNP-IDs with data from annotation file")
    annotationFile <- annotationFile[!is.na(annotationFile[,annofile.SNP.columnName]), ] # remove entries without SNP identifier
    annotationFile <- annotationFile[!duplicated(annotationFile[,annofile.SNP.columnName]), ] # remove duplicated SNP entries (affects control probes)
    
    if (!is.null(annofile.columns)) {
      if(!annofile.SNP.columnName %in% annofile.columns) {annofile.columns <- c(annofile.SNP.columnName, annofile.columns)}
      annotationFile <- annotationFile[, annofile.columns, drop=F]
      } 
    
    # merge data with annotation file
    names(annotationFile) <- paste0("Annofile_", names(annotationFile))
    data <- merge(data, annotationFile, by.x= data.SNP.columnName, by.y=paste0("Annofile_", annofile.SNP.columnName), all.x=T)
    
    # rename column rows if overwritten and remove temporara row.name columns
    if(data.SNP.columnName == "row.names") {
    rownames(data) <- data$Row.names
    data <- data[, !(names(data) %in% c("Row.names"))]
    }
  }
 
 
  
  # re-order data and remove temporary index column
  data <- data[order(data$temp.index), !(names(data) %in% c("temp.index"))]
  
  return(data) 

} # end function definition
