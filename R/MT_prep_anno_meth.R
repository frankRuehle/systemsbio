#' Preprocessing and annotation of methylation array data
#' 
#' An RGChannelSet object is preprocessed, mapped to the genome and annotated by the respective annotation package. 
#' 
#' This function takes an RGChannelSet as input and transforms it to a MethylSet object using one of 3 possible preprocessing procedures 
#' from \code{minfi} package:
#' \itemize{
#'   \item \code{preprocessRaw} (no normalisiation)
#'   \item \code{preprocessIllumina} applies background subtraction as well as control normalization like in Illumina Genome studio
#'   \item \code{preprocessSWAN} applies SWAN normalization (Subset quantile Within-Array Normalization) method on a preprocessed 
#'                               MethylSet to reduce the technical variability of Type I and Type II Illumina probes. 
#'                               For this, Illumina preprocessed MethylSet is used if given, the raw preprocessed MethylSet otherwise.  
#'   }
#' If no preprocessing is selected, the input object is expected to be already preprocessed. Next, the MethylSet object is
#' mapped to the genome (Loci which cannot be mapped to a genomic position are dropped) and converted to a GenomicRatioSet. 
#' Finally, annotation data is added from the respective annotation package.
#' CpG sites overlapping with SNPs optionally can be removed from the dataset. If these sites are of special interest, 
#' the returned \code{GenomicRatioSet} can be specified to contains these CpG-sites only.
#'     
#' @param RGset RGChannelSet object with methylation data.
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param preprocessing character with preprocessing procedures to apply. Any of \code{"raw", "illumina", "swan"}.
#'                      If multiple procedures are selected, the one with higest priority is used: \code{"raw" < "illumina" < "swan"}.
#'                      If "swan" is selected, at least one more procedure must be selected, too.
#'                      If \code{NULL}, input object is expected to be already preprocessed.
#' @param addAnno logical. If TRUE, annotation data is added from annotation package.
#' @param SNPloci2drop character with SNP types, either \code{"CpG"} (SNPs overlapping the CpG site) and/or \code{"SBE"} 
#'                     (SNP at the single base extension of CpG-site). CpG sites overlapping these SNPs are removed from the dataset. 
#'                     If \code{NULL}, no CpG sites are removed. 
#' @param MAF2drop numeric with minimum minor allele frequency of SNPs to drop if \code{SNPloci2drop} is specified. 
#' @param return.EpiSNPs.only logical. If TRUE, the returned GenomicRatioSet contains only CpG-sites overlapping with SNP locations
#'                            (usefull to focus on SNPs influencing methylation processes).
#'   
#'   
#' @return GenomicRatioSet (containing Beta and M-values)
#' 
#' @author Frank Ruehle
#' 
#' @export


# Annotierung der CpG-sites:

prep_anno_meth <- function(RGset, 
                    projectfolder = "MT",
                    projectname = NULL,
                    preprocessing = c("illumina", "swan"),
                    addAnno = T, 
                    SNPloci2drop = c("CpG", "SBE"), 
                    MAF2drop = 0, 
                    return.EpiSNPs.only = F) {
 

  # load required packages. 
  annotation.package <- c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  pkg.bioc <- c("minfi", "limma")
  pks2detach <- attach_package(pkg.bioc = c(annotation.package, pkg.bioc))

  #orig_par <- par(no.readonly=T)      # make a copy of current settings
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T)}
  
  if (!is.null(projectname) && !grepl("_$", projectname)) {projectname <- paste0(projectname, "_")} 
  
 
  ## select desired preprocessing procedures and determine reference object for SWAN normalisation if applicable
  if(!is.null(preprocessing)) {
    
  objects.preprocess <- c("MSet.raw", "MSet.normIllumina", "MSet.normSwan")
  objects.preprocess <- factor(objects.preprocess, levels=objects.preprocess, ordered=T)
  objects.preprocess <- objects.preprocess[grepl( paste(preprocessing, collapse="|"), levels(objects.preprocess), ignore.case = T)]
  objects.preprocess <- sort(objects.preprocess, decreasing=T)


  if(objects.preprocess[1] == "MSet.raw" || (objects.preprocess[!objects.preprocess %in% "MSet.normSwan"][1] == "MSet.raw")) {
    cat("\nPrepare MethylSet from RGChannelSet by raw preprocessing\n")
    MSet <- preprocessRaw(RGset)
  } else {
    if(objects.preprocess[1] == "MSet.normIllumina" || (objects.preprocess[!objects.preprocess %in% "MSet.normSwan"][1] == "MSet.normIllumina")) {
      cat("\nPrepare MethylSet from RGChannelSet by Illumina normalisation\n")
      MSet <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 1) # as in GenomStudio.
      } 
  }
    
  if(objects.preprocess[1] == "MSet.normSwan") {
    if(length(objects.preprocess) == 1) {stop("Give an additional item in 'preprocessing' as necessary for SWAN normalisation")
      } else {
        cat(paste("\nPrepare SWAN normalised MethylSet from", objects.preprocess[!objects.preprocess %in% "MSet.normSwan"][1], "\n"))
        MSet <- preprocessSWAN(RGset, MSet)
      } 
    }

 } else {
   if(class(RGset) == "RGChannelSet") {stop("If no preprocessing is selected, input object must be preprocessed in advance!")}
   MSet <- RGset} # end if !is.null(preprocessing).
  
  

  
  # Mapping Ilumina methylation array data to the genome using an annotation package. 
  # This function silently drops loci which cannot be mapped to a genomic position.
  # Converts MethylSet object to GenomicMethylSet object
  cat("\nMapping Ilumina methylation array data to the genome using an annotation package.\n")
  cat("Loci which cannot be mapped to a genomic position are dropped.\n")
  MSet <- mapToGenome(MSet, mergeManifest = T) # adds features "assayStrand" "Name" "AddressA" "AddressB" "ProbeSeqA" "ProbeSeqB" "Type" "NextBase" "Color"
  # addSnpInfo adds information on which probes contain SNPs 
  MSet <- addSnpInfo(MSet, snpAnno = NULL) # adds features "Probe_rs" "Probe_maf" "CpG_rs" "CpG_maf" "SBE_rs" "SBE_maf"
  # Converting methylation data from methylation and unmethylation channels, to ratios (Beta and M-values)
  GRatioSet <- ratioConvert(MSet)



 #annotate GRatioSet with annotation object
  if(addAnno) {
  cat("\nGet annotations for GenomicRatioSet.\n")
  methanno <- minfi::getAnnotation(GRatioSet, what = "everything", lociNames = NULL, orderByLocation = FALSE, dropNonMapping = FALSE)
    # create SYMBOL column from methanno$UCSC_RefGene_Name, which containes unique genes only. Redundant Symbols from multiple transcripts are removed
    symbol <- strsplit(methanno$UCSC_RefGene_Name, split=";")
    symbol <- lapply(symbol, unique)
    symbol <- sapply(symbol, paste, collapse=";")
    methanno$SYMBOL <- symbol
  # add additional feature data from methanno
  annoColumns2add <- !(names(methanno) %in% names(mcols(GRatioSet)))
  mcols(GRatioSet) <- DataFrame(mcols(GRatioSet), methanno[rownames(GRatioSet), annoColumns2add])
  }

  # Overview of SNP types possibly influencing methylation measurement
  probe_SNPs <- !is.na(mcols(GRatioSet)$Probe_rs)   # SNP overlapping the probe (in Typ I probes a SBE_SNP may count as probe_SNP)
  CpG_SNPs <- !is.na(mcols(GRatioSet)$CpG_rs)       # SNPs overlapping the CpG site
  SBE_SNPs <- !is.na(mcols(GRatioSet)$SBE_rs)       # the single base extension of CpG-site (for Type II probes identical with CpG-SNP, for Type I probe 1 base behind CpG-Site)
  Assaytype.I <- mcols(GRatioSet)$Type=="I"         # 2 different probes for methylated or unmethylated (either both red or both green)
  Assaytype.II <- mcols(GRatioSet)$Type=="II"       # 1 probe type per CpG. Methylated = green; unmethylated = red. Majority of probes is type II.
  
  filenameVennDiagram <- file.path(projectfolder, paste0(projectname, "Venn_SNPtypes.png"))
  cat("print VENN Diagramm of Loci with interfering SNPs to:", filenameVennDiagram, "\n")
  png(filename=filenameVennDiagram, width = 210, height = 210, units = "mm", res=600)
  par(mfrow=c(2,1))
  vc1 <- vennCounts(cbind(probe_SNPs, CpG_SNPs, SBE_SNPs)) 
  vc2 <- vennCounts(cbind(Assaytype.II, probe_SNPs, CpG_SNPs, SBE_SNPs))  
  vennDiagram(vc1)
  vennDiagram(vc2)
  par(mfrow=c(1,1))
  dev.off()
  
  GRSet.epiSNP <- GRatioSet[CpG_SNPs,] # returned if return.EpiSNPs.only is TRUE
  


  # Removing epigenetic SNPs from main data object (if SNPloci2drop is set)
  # object needs to be of class 'GenomicMethylSet' or 'GenomicRatioSet' (not MethylSet)
  if(!is.null(SNPloci2drop)) {
    cat("\nLoci containing SNPs with MAF >",MAF2drop, "are removed:", SNPloci2drop,"\n")
    GRatioSet <- dropLociWithSnps(GRatioSet, snps = SNPloci2drop, maf = MAF2drop) # removing loci with SNPs, optionally based on their MAF.
    }
  
 
  # Detaching libraries not needed any more
  detach_package(unique(pks2detach))
  
  
  if (!return.EpiSNPs.only) {
  cat("\nreturning annotated GenomicRatioSet object. Unmapped probes and SNP loci are removed.\n")
  return(GRatioSet)
    } else {
      cat("\nreturning annotated GenomicRatioSet containing loci with SNP at CpG-site only.\n")
      return(GRSet.epiSNP) }

}



### Data access of MethylSet-Opjekt
# getMeth(MSet)[1:4,1:3]   # returns unlogged methylation channels
# getUnmeth(MSet)[1:4,1:3]
# getBeta(MSet, type = "Illumina")[1:4,1:3]  # gets "beta"-values which are values between 0 and 1
# getM(MSet)[1:4,1:3]    # gets "M-values" computed as logit(beta)
# getCN(MSet)[1:4,1:3]  #get copy number values which are defined as the sum of the methylation and unmethylation channel.
# getManifest(MSet)
# sampleNames(MSet)
# pData(MSet)[1:4,1:3]
# featureNames(MSet)[1:10]
# preprocessMethod(MSet)
# annotation(MSet)
# getIslandStatus(MSet)[1:10]

## Data access of GenomicMethylSet object:
# length(GMset)
# seqnames(GMset)
# ranges(GMset)
# names(GMset)
# strand(GMset)
# mcols(GMset, use.names=T) # Get or set the metadata columns. If use.names=TRUE and the metadata columns are not NULL, 
#                                 then the names of x are propagated as the row names of the returned DataFrame object. 
# seqinfo(GMset)
# seqlevels(GMset)
# seqlengths(GMset)
# isCircular(GMset)
# genome(GMset)
# seqlevelsStyle(GMset)
# mcols(GMset)[1:5,1:5]

