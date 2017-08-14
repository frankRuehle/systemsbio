#' Read Genotype Data Using GenABEL.
#'
#' \code{readGeno} reads genotype data stored in ped/map files, adjusts format for GenABEL package and processes it for further analysis.
#'
#' If a GenABEL raw data file given in \code{existingRawFile} already exists, this file is directly loaded into a gwaa object. 
#' Otherwise this file has to be created first from the supplied ped and map file. For this, the map file in \code{mapFilename} is 
#' loaded, the linkage column (column 3) removed if present and a headerline is added. If an Array description file is supplied 
#' in \code{GTarrayDescriptionFile}, an extended map file is generated for GenABEL including strand and allele coding 
#' information. Otherwise strand information is set to unknown. If indicated in \code{updateSNPpos} SNP 
#' chromosomal postions can be updated according either to biomaRt (if a biomaRt object is given in \code{snpmart}) or 
#' to \code{GTarrayDescriptionFile}. Mind that GenABEL does not allow strings of characters as alleles in the 
#' corresponding ped-file (as for Indels), but just 1-character-alleles or "-". 
#' 
#' After loading the GenABEL raw file, phenotype information is added from the supplied covariates file \code{covarFilename}.
#' I.e. no phenotype or gender data is considered from initial ped-file. No white space allowed in phenotype entries.
#' No missing data allowed in gender column given in \code{sexvar}. All samples from geno file must have entries in pheno-file. 
#' For 1/2-coded binary phenotypes additional 0/1-coded phenotypes are generated (stored with suffix "01"). 
#' If gender is 1/2 (1=male, 2=female) coded in \code{sexvar}, a new variable \code{sex} 
#' is generated 0/1 coded (1=male, 0=female) as required by GenABEL. Optionally, numeric gonosome names can be converted to 
#' character names (e.g. "X", "Y"). If \code{removeNullPositions} is set TRUE, variants without valid coordinates are
#' removed from the generated \code{gwaa} dataset. 
#'
#' @param genoFilename Character with path to ped file.
#' @param mapFilename Character with path to map file.
#' @param covarFilename Character with path to covariates file (gender information is needed).
#' @param existingRawFile Character with path to an optionally already existing GenABEL Raw which will preferably loaded
#' and supersedes any other files. Omitted if NULL.  
#' @param covarsampleID Character with column name of sample IDs in covar file.
#' @param sexvar Character with column name indicating gender in covar file. 
#' @param projectfolder Character containing path to output folder (will be generated if not existing).
#' @param projectname Character used as suffix for output files.
#' @param organism Character with name of organism (e.g. "human").
#' @param Gonosomes2char Boolean. If TRUE, chromosome names 23, 24, 25, 26 will be converted to X, Y, XY, MT (for human only).
#' @param removeNullPositions Boolean. If TRUE, SNPs with Chr=0 or Pos=0 or NA will be removed from generated \code{gwaa} object
#' @param updateSNPpos Character with value "biomaRt" or "descriptionfile". For "biomaRt", chromosome and SNP positions 
#' in map-file will be updated by biomaRt. If "descriptionfile", they will be updated by the supplied \code{GTarrayDescriptionFile}.
#' If NULL chromosome and SNP positions are not updated.
#' @param snpmart biomaRt object to be used for updating SNP positions (or NULL). 
#' @param GTarrayDescriptionFile Optional character with path to Illumina array description file. If given, 
#' strand and allele coding data will be added to map file. 
#' @param GTarrayDescription.lines2skip.start Numeric with number of rows to skip when loading \code{annotationFile} or 
#'                         regular expression for character string to identify corresponding row number 
#'                         to be skipped, e.g. \code{[Assay]} in Illumina annotation files.
#' @param GTarrayDescription.lines2skip.end Numeric with number of rows to read when loading \code{annotationFile} or 
#'                        regular expression for character string to identify corresponding row number to be read,
#'                        e.g. \code{[Controls]} in Illumina annotation files as start of annotation of control probes.
#'                        All rows from that number on (incl. \code{GTarrayDescription.lines2skip.end}) are skipped.  
#' @param GTarrayDescription.colname.identifier Character with colnames for SNP identifier in Array description file.
#' @param GTarrayDescription.colname.coding Character with colnames for allele coding information
#' @param GTarrayDescription.colname.strand Character with colnames for strand information.
#' 
#' @return GenABEL gwaa object
#' 
#' @author Frank Ruehle
#' 
#' @export readGeno


readGeno <- function (genoFilename, 
                    mapFilename, 
                    covarFilename, 
                    existingRawFile = NULL,    
                    covarsampleID = "IID", 
                    sexvar = "Gender",      
                    projectfolder= "GT",
                    projectname = NULL,
                    organism = "human", 
                    Gonosomes2char = T,     
                    removeNullPositions = T, 
                    updateSNPpos = NULL,     
                    snpmart = useMart("ENSEMBL_MART_SNP", host = "feb2014.archive.ensembl.org", dataset="hsapiens_snp"), # GRCh37.p13, latest hg19 annotation
                    GTarrayDescriptionFile = NULL,
                      GTarrayDescription.lines2skip.start = "\\[Assay\\]",
                      GTarrayDescription.lines2skip.end = "\\[Controls\\]",
                      GTarrayDescription.colname.identifier = "Name",  
                      GTarrayDescription.colname.coding = "SNP",
                      GTarrayDescription.colname.strand = "RefStrand",
                      GTarrayDescription.colname.chromosome = "Chr",
                      GTarrayDescription.colname.position = "MapInfo"
                    ) {
                    
 
  
  # Create output directory if not yet existing
  if (!file.exists(projectfolder)) {dir.create(file.path(projectfolder), recursive=T)}

  ## install/load required packages from CRAN and Bioconductor
  pkg.bioc <- c("biomaRt")
  pkg.cran <- c("GenABEL")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  projectname <- if (!is.null(projectname)) {paste0(projectname, "_")} else {""}
  
  
  #### Modification of covariate file
  if(!is.null(covarFilename)) {
  # Reading as latin1-encoded for compability with Linux environment
  cat("\nProcessing covariates file for GenABEL\n")
  covar <- read.delim(file=covarFilename, header=T, sep="\t", stringsAsFactors =F, 
                      na.strings = c("NA", "-9", "", " "), fileEncoding="latin1") 
  covar[,covarsampleID] <- factor(gsub(" ", "", covar[,covarsampleID])) # remove whitespace from sample names
  
  # GenABEL needs a 0/1-coded column named "sex" (1=male, 0=female)! 
  if(!is.null(sexvar)){
    if(any(covar[,sexvar]==2)) {
      cat("\n0/1-coded variable sex (1=male, 0=female) generated from 1/2 coded variable", sexvar, "(1=male, 2=female)\n")
      if(sexvar=="sex") {covar$sex_unmod <- covar$sex}
      covar$sex <- covar[,sexvar]
      covar[!is.na(covar$sex) & covar$sex==2,"sex"] <- 0 # Females recoded from 2 to 0
    } else {covar$sex <- covar[,sexvar]}
  } else{stop("\n\nsex variable in covariates file needed!\n")}    
  

    # Covariate file must start with unique identifier in first column
  covar <- covar[,c(covarsampleID, setdiff(names(covar), covarsampleID)) ]
  
  # no white space allowed in character elements -> replaced by underline
  for(c in names(covar)) {
   if(is.character(covar[,c])) {
     covar[,c] <- gsub(" ", "_", covar[,c])
     }
  }
  
  
  modcovarFilename <- file.path(projectfolder, paste0(projectname, "GenABEL.cov"))
  cat("\nWriting covariates file to", modcovarFilename, "\n")
  write.table(covar, file= modcovarFilename, quote=F, row.names=F, sep="\t")
  } else {stop("Covariate file needed for affection and gender information!")}
  

  
  
#### Modification of map file for GenABEL conversion, if not done yet.
if (!is.null(existingRawFile) && file.exists(file.path(existingRawFile))) {
  cat("\n", file.path(existingRawFile), "already exists. Converting is omitted.\n")
  } else { 
 
  cat("\nReading map file", if(!is.null(GTarrayDescriptionFile)){"and Array description file"}, "\n")    

    # Introducing headerline and removing linkage column if present (required for GenABEL)
    mapfilemod <- read.table(mapFilename, header = FALSE, stringsAsFactors = F)
    if (length(mapfilemod)==4) {mapfilemod <- mapfilemod[,-3]}
    names(mapfilemod) <- c("chr", "name", "pos")
    
    # purify map file from optional characters such as "chr", "---", "null"
    mapfilemod$chr <- sub("chr", "", mapfilemod$chr, ignore.case = T)
    mapfilemod$chr[!is.na(mapfilemod$chr) & mapfilemod$chr=="---"] <- 0
    mapfilemod$pos[!is.na(mapfilemod$pos) & mapfilemod$pos=="null"] <- 0
    
    
    # Optionally replacing map file columns
    # if biomaRt and/or Array description file is supplied, chromosome and pos is updated accordingly.
    # if Array description file is supplied, strand and coding information is added to map file.
    # Strand information columns in Illumina manifest file: "IlmnStrand", "SourceStrand", "RefStrand."
     mapfilemod.anno <- basic_SNP_annotation(mapfilemod,
                   max.SNPs.per.biomaRt.call = 30000,
                   data.SNP.columnName = "name",
                   snpmaRt = if(is.null(updateSNPpos)){NULL} else {if(tolower(updateSNPpos) == "biomart") {snpmart} else {NULL}} ,
                   biomaRt.SNP.columnName = "refsnp_id", 
                   biomaRt.filter = "snp_filter",
                   biomaRt.attributes.groupColumns  = c("refsnp_id", "chr_name", "chrom_start"),
                   biomaRt.attributes.summarized = NULL,
                   annotationFile = GTarrayDescriptionFile,
                   annofile.SNP.columnName = GTarrayDescription.colname.identifier,
                   lines2skip.start = GTarrayDescription.lines2skip.start,
                   lines2skip.end = GTarrayDescription.lines2skip.end,
                   annofile.columns = NULL #c(GTarrayDescription.colname.identifier, GTarrayDescription.colname.strand, 
                                        #GTarrayDescription.colname.coding, GTarrayDescription.colname.chromosome,
                                        #GTarrayDescription.colname.position)
                  )
 
  
   if(!is.null(updateSNPpos)) {
     if(tolower(updateSNPpos) == "biomart") {
         if(!is.null(snpmart)) {
           cat("\nChromosome and SNP Positions from biomaRt are used for map file")
           mapfilemod <- mapfilemod.anno[, c("SNPMart_chr_name", "name", "SNPMart_chrom_start")]
           names(mapfilemod) <- c("chr", "name", "pos")
           } else {cat("\no biomaRt specified!")}
       }
    }
 
    if(!is.null(GTarrayDescriptionFile)) {
      if(all(c(paste0("Annofile_", GTarrayDescription.colname.strand), 
               paste0("Annofile_", GTarrayDescription.colname.coding)) %in% names(mapfilemod.anno))) {
        cat("\nAdding strand and allele coding data from description file to map file")
        mapfilemod <- data.frame(mapfilemod, 
                               strand = mapfilemod.anno[,paste0("Annofile_", GTarrayDescription.colname.strand)],
                               coding = mapfilemod.anno[,paste0("Annofile_", GTarrayDescription.colname.coding)])
        } else {cat("\nno strand and/or allele coding information found.")}
      if(!is.null(updateSNPpos)) {
        if(tolower(updateSNPpos) == "descriptionfile") {
          if(all(c(paste0("Annofile_", GTarrayDescription.colname.chromosome), 
                   paste0("Annofile_", GTarrayDescription.colname.position)) %in% names(mapfilemod.anno))) {
                cat("\nChromosome and SNP Positions from description file are used for map file")
              mapfilemod$chr <- mapfilemod.anno[, paste0("Annofile_", GTarrayDescription.colname.chromosome)]
              mapfilemod$pos <- mapfilemod.anno[, paste0("Annofile_", GTarrayDescription.colname.position)]
              } else {cat("\nno chromosome and/or BP information found.")}
           } 
        }
    }
    
    # if alternative annotation intoduces NAs during merging, they are set to 0.
    mapfilemod[is.na(mapfilemod$chr) | mapfilemod$chr=="NA" | mapfilemod$chr=="", "chr"] <- 0
    mapfilemod[is.na(mapfilemod$pos) | mapfilemod$pos=="NA" | mapfilemod$pos=="", "pos"] <- 0
    
    
 
filename.mapfilemod <- file.path(projectfolder, paste0(projectname, "GenABEL.map"))
cat("\nWriting modified map file to", filename.mapfilemod, "\n")
mapfilemod$pos <- format(mapfilemod$pos, scientific = FALSE) # otherwise some BP pos are stored in scientific format!
write.table(mapfilemod, filename.mapfilemod, quote=F, row.names=F, sep="\t")


##### Converting from Linkage-Format (ped/map-file) to GenABEL-Format gwaa.data, if not done yet
      cat("\nConverting SNP data to GenABEl format and saving to", file.path(projectfolder, paste0(projectname, "GenABEL.raw")), "\n")
      convert.snp.ped(pedfile = genoFilename,
                mapfile = filename.mapfilemod,
                mapHasHeaderLine=TRUE,
                format = "premakeped",
                traits = 1,  # How many traits are specified in the pedigree file? Usually 1 (affection) or 2 (affection and liability)
                wslash=FALSE, # assumed that alleles are separated with space.
                outfile = file.path(projectfolder, paste0(projectname, "GenABEL.raw")),
                strand = if("strand" %in% names(mapfilemod)) {"file"} else {"u"}) 
      # if strand = "file", strand information from extended map file is used. Otherwise strand = unknown.
      existingRawFile <- file.path(projectfolder, paste0(projectname, "GenABEL.raw"))
} # end of condition if GenABEL.raw already exists



# Laden der konvertierten Genotypen in ein gwaa-Object                
cat("\nLoading gwaa object from", file.path(projectfolder, "GenABEL.raw"), "and", modcovarFilename, "\n")
genos <- load.gwaa.data(phe = modcovarFilename, 
                        gen = existingRawFile, 
                        force = T, id = covarsampleID) 



#### transform 1/2 coded binary covars to 0/1 coded
cat("\nBinary 1/2 coded variables are converted to 0/1-coded.\n")

# which phenotypes are binary 1/2 coded?
binpheno12 <- apply(phdata(genos),2,function(x) { all(na.omit(x) %in% 1:2) })
binpheno12[sexvar] <- FALSE # sexvar needed separate coding (see above)
binpheno12 <- binpheno12[binpheno12==TRUE]

if (any(binpheno12)) { # is FALSE if length(binpheno12) == 0
  for (i in 1:length(binpheno12)) {
    genos <- add.phdata(genos, newph=phdata(genos)[,names(binpheno12)[i]]-1, name=paste0(names(binpheno12)[i],"01"))
  }
}




# Rename Chromosoms 23, 24, 25, 26, in X, Y, XY, MT (if human)
if(Gonosomes2char & grepl("(hsapiens)|(human)", tolower(organism))) {
cat("\nRenaming chromosomes 23, 24, 25, 26, in X, Y, XY, MT\n")
levels(genos@gtdata@chromosome)[levels(genos@gtdata@chromosome) == "23"] <- "X"
levels(genos@gtdata@chromosome)[levels(genos@gtdata@chromosome) == "24"] <- "Y"
levels(genos@gtdata@chromosome)[levels(genos@gtdata@chromosome) == "25"] <- "XY"
levels(genos@gtdata@chromosome)[levels(genos@gtdata@chromosome) == "26"] <- "mt"
}


### remove Chr=0 and pos=0
if(removeNullPositions) {
  cat("\nRemoving SNPs chr=0 or pos=0.")
  chromo0 <- chromosome(genos)==0
  cat("\n", sum(chromo0), "SNPs removed with chromosome = 0.")
  genos <- genos[,!chromo0]

  Position0 <- map(genos)==0
  cat("\n", sum(Position0), "SNPs removed with bp position = 0.")
  genos <- genos[,!Position0]
  }


# ### data access
# # sample data
# phdata(genos)
# nids(genos)
# idnames(genos)
# male(genos)
# 
# 
# # snp data
# gtdata(genos[1:10, 1:10])
# as.character(genos@gtdata)[1:10,1:10]
# as.numeric(genos@gtdata)[1:10,1:10]  # count of effector allele
# as.genotype(gtdata(genos[1:10, 1:10])) # format of "genetics" package
# as.hsgeno(gtdata(genos[1:10, 1:10]))  # format used by "haplo.stats"
# nsnps(genos)
# snpnames(genos)[1:10]
# chromosome(genos)[1:10]
# map(genos)[1:10]
# coding(genos)[1:10]  # First Allele Pos: reference allele, second Allele Pos: effector allele
# refallele(genos)[1:10]
# effallele(genos)[1:10]
# GenABEL::strand(genos)[1:10]  # strand() is masked by BiocGenerics
# 
# 




# Display summary data
# cat("\nSample Summary Data\n")
# print(idsum)
cat("\nMarker Summary Data\n")
print(descriptives.marker(genos))
cat("\nTrait Summary Data\n")
print(descriptives.trait(genos))


return(genos)

}



