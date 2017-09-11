#' Finding differentially methylated CpG positions or regions
#' 
#' Use \code{minfi}-package to determine either differentially methylated positions (DMPs, \code{dmpFinder}) and/or 
#' differentially methylated regions (DMRs, \code{bumphunter}).
#' 
#' The \code{dmpFinder} function from \code{minfi} identifies differentially methylated CpG-sites using
#' linear regression for continuous phenotypes and F-test for categorical phenotypes. In case of a categorical
#' phenotype with many groups (e.g. experimental groups), F-test is applied over all groups. 
#' Tests are performed on M values, which are logit transformed Beta values 
#' (beta = Methylated allele intensity / (Unmethylated allele intensity + Methylated allele intensity + 100)).
#' Annotation data of \code{GRset} is added to the result table.  
#' 
#' Instead of looking for association between a single genomic location and a phenotype of interest, 
#' \code{bumphunter} looks for genomic regions that are differentially methylated (beta values). In the context of the 450k array, 
#' the algorithm first defines clusters of probes such that two consecutive probe locations in the cluster
#' are not separated by more than distance \code{mapGap}. Briefly, the algorithm first computes a t-statistic for beta values
#' at each genomic location, with optional smoothing. Then, it defines a candidate region to be a cluster of probes
#' for which all the t-statistics exceed a predefined threshold. To test for significance of the candidate regions, 
#' the algorithm uses permutations given in \code{nResamples}.
#' 
#' 
#' @param GRset GenomicRatioSet or MethylSet
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param phenotype character with phenotype for differential methylation analysis. Must be given within the 
#'                  phenotype data of \code{GRset}.
#' @param type.covar character with type of phenotype variable, either \code{continuous} or \code{categorical}.
#' @param returnResults character indicating differentially methylated elements to return. 
#'                      Either \code{"DMP"} for differentially methylated positions using \code{dmpFinder}
#'                      or \code{"DMR"} for differentially methylated regions using \code{bumphunter} or both.
#' @param qCutoff numeric with q-value (FDR) threshold for DMPs to be reported by \code{dmpFinder}.
#' @param bumpcutoff numeric start value to find the appropriate \code{bumphunter} cutoff.
#'                         A cutoff 0.1 corresponds to 10% difference on the Beta-values.
#'                         If the cutoff results in more than 10000 candidate DMRs, the cutoff is automatically increased
#'                         before applying permutation testing to avoid excessive computation time.
#' @param nResamples numeric, number of resamples to use when computing null distributions with \code{bumphunter}.
#' @param maxGap numeric with maximum location gap used to define clusters of probes for \code{bumphunter}.
#' 
#'   
#' @return list with up to 2 elements depending on selection in \code{returnResults}. 
#'         Result tables are also stored as side-effects in \code{projectfolder}
#' \itemize{
#'   \item \code{DMP} dataframe of CpG positions sorted by differential methylation p-value annotated by respective annotation package.
#'   \item \code{DMR} An object of class \code{bumps} with the components: \code{tab, coef, fitted, pvaluesMarginal, null, algorithm}.
#'                    See help(bumphunter) for details. The dataframe in \code{tab} is annotated with meta data for all CpG-sites residing
#'                    within the respective DMR.
#' }
#' 
#' @author Frank Ruehle
#' 
#' @export




diffMinfi <- function(GRset, 
                      projectfolder = "MT",
                      projectname = NULL,
                      phenotype, 
                      type.covar = "categorical", 
                      returnResults = c("DMP", "DMR"),
                      qCutoff= 1e-3,  
                      bumpcutoff = 0.1, 
                      nResamples = 1000, 
                      maxGap = 500
                      ) {

  
  # load required packages. 
  annotation.package <- c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  pkg.bioc <- c("minfi", "limma", "GenomicRanges")
  pks2detach <- attach_package(pkg.bioc = c(annotation.package, pkg.bioc))
  
  
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }

  
  if(class(GRset) == "GenomicRatioSet") { # GenomicRatioSet and MethylSet objects have different functions to access feature annotation
    funFeatureAnno <- getAnnotation} else {
      funFeatureAnno <- mcols
    }
  
  
  Mvalues <- getM(GRset) # get matrix with M values
  beta <- getBeta(GRset) # get matrix with beta values
  results <- list()

  ### dmpFinder
  # Input: MethylSet or a matrix (but not GenomicMethylSet)
  if("DMP" %in% returnResults) {
  
  cat("\nRunning dmpFinder on M values to identify differentially methylated positions (DMPs)\n")
  if (!file.exists(file.path(projectfolder, "minfi_DMP"))) {dir.create(file.path(projectfolder, "minfi_DMP")) }
  
  dmp <- dmpFinder(Mvalues, pheno=pData(GRset)[,phenotype], type= type.covar, qCutoff=qCutoff) 
  
  ## merging dmp with annotation from GRset to DataFrame 
  DMPannot <- merge(dmp, funFeatureAnno(GRset), by=0, all.x=T)  # by=0 means rownames. , use.names=T
  DMPannot <- DMPannot[order(DMPannot$pval),]
  names(DMPannot)[1] <- "CpG"
  
  dmp.filename <- file.path(projectfolder, "minfi_DMP", paste0(projectname,"_minfi_DMP_", type.covar, "_", phenotype, ".txt"))
  cat("\nWriting annotated dmpFinder results to:", dmp.filename, "\n")
  write.table(DMPannot, file=dmp.filename, sep="\t", quote=F, row.names = F)
  
  # Use the plotCpG function to plot methylation levels at individual positions (here: top 4 diff methylated, M and beta values)
  if(nrow(dmp)>=1) {
    plotdmp.filename <- gsub(".txt", "_topHits.png", dmp.filename)
    cat("\nPlotting methylation levels at top differentially methylated positions (M- and beta-values) to:", plotdmp.filename, "\n")
    
    cpgs <- rownames(dmp)[1:min(nrow(dmp),4)]
    png(filename=plotdmp.filename, width = 210, height = 210, units = "mm", res=600)
      par(mfrow=c(2,2))
      plotCpg(Mvalues, cpg=cpgs, pheno=pData(GRset)[,phenotype], ylab="M values")
      plotCpg(beta, cpg=cpgs, pheno=pData(GRset)[,phenotype], ylab="beta values")
    dev.off()
    par(mfrow=c(1,1))
  }

  results[["DMP"]] <- DMPannot
}


  ### Bumphunter to find differentially methylated regions (DMRs)
  if("DMR" %in% returnResults) {
    cat("\nRunning bumphunter on beta values to identify differentially methylated regions (DMRs)\n")
    if (!file.exists(file.path(projectfolder, "minfi_DMR"))) {dir.create(file.path(projectfolder, "minfi_DMR")) }
    
  # get column names for chromosome and basepair position
  columnChr <- grep("chr", names(funFeatureAnno(GRset)), ignore.case = T, value=T)[1]
  columnPos <- grep("(pos)|(bp)", names(funFeatureAnno(GRset)), ignore.case = T, value=T)[1]
  cat("\nColumn names used for defining regions: Chromosome:", columnChr, "; basepair position:", columnPos, "\n") 
    
  # Define your phenotype of interest:
  # Categorical phenotype with multiple levels is coerced to a single numeric column for design matrix
  if(is.numeric(pData(GRset)[,phenotype])) {pheno <- pData(GRset)[,phenotype]} else {
    cat("\nPhenotype", phenotype, "is converted to numeric by factor levels!\n")
    pheno <- as.numeric(factor(pData(GRset)[,phenotype]))
    print(unique(data.frame(before_conversion=pData(GRset)[,phenotype], converted=pheno)))  
    }
  
    designMatrix <- model.matrix(~ pheno)
  
  # Run the algorithm with B = 0 permutation on the Beta-values, with a medium difference
  # cutoff of e.g. 0.1 (which corresponds to 10% difference on the Beta-values):
  cat("\nTesting if cutoff value for bumphunter results in > 10000 candidate regions.\n")
  dmrs <- bumphunter(beta, design = designMatrix, chr=funFeatureAnno(GRset)[,columnChr], pos=funFeatureAnno(GRset)[,columnPos], 
                     maxGap=maxGap, cutoff = bumpcutoff, B=0, type="Beta")
  
  if(nrow(dmrs$table>=1)) { # otherwise no bumps to use.
      # If the number of candidate bumps is large (> 10000), increase the cutoff to reduce the number of candidate bumps.
      maxBumps <- 10000
      while(nrow(dmrs$table) > maxBumps) {
        bumpcutoff <- bumpcutoff + 0.05
        dmrs <- bumphunter(beta, design = designMatrix, chr=funFeatureAnno(GRset)[,columnChr], pos=funFeatureAnno(GRset)[,columnPos], 
                           maxGap=maxGap, cutoff = bumpcutoff, B=0, type="Beta")
      }

      # Once you have decided on the cutoff, run the algorithm with a large number of permutations, e.g. B = 1000
      cat("\nRunning bumphunter with cutoff:", bumpcutoff, "(", 100*bumpcutoff, "% difference on the Beta-values) and", nResamples, "permutations.\n")
      
      dmrs <- bumphunter(beta, design = designMatrix, chr=funFeatureAnno(GRset)[,columnChr], pos=funFeatureAnno(GRset)[,columnPos], 
                         maxGap=maxGap, cutoff = bumpcutoff, B= nResamples, type="Beta")
      
      
 
      ### annotation with data from annotation file
      cat("\nAnnotate bumphunter regions with data from CpG sites.\n")
      cpg.ranges <- granges(GRset)
        seq.columns <- grepl("seq", names(mcols(cpg.ranges)), ignore.case=T)
        mcols(cpg.ranges) <- mcols(cpg.ranges)[!(seq.columns)] # remove (long) sequence annotation data and rename forbidden meta column names 
        forbidden.meta.columnnames <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element")
        names(mcols(cpg.ranges))[names(mcols(cpg.ranges)) %in% forbidden.meta.columnnames] <- paste0(names(mcols(cpg.ranges))[names(mcols(cpg.ranges)) %in% forbidden.meta.columnnames], "_anno")
      
      dmr.ranges <- makeGRangesFromDataFrame(dmrs$table, keep.extra.columns = TRUE) 
      
      
      dmr.ranges.anno <- subsetByOverlaps.keepAllMeta(dmr.ranges, cpg.ranges)
      
      dmr.ranges.anno.df <- as.data.frame(dmr.ranges.anno)

      DMR.filename <- file.path(projectfolder, "minfi_DMR", paste0(projectname,"_minfi_DMR_", type.covar, "_", phenotype, ".txt"))
      cat("\nWriting bumphunter results to:", DMR.filename, "\n")
      write.table(dmr.ranges.anno.df[order(dmr.ranges.anno.df$fwer),], sep="\t", quote=F, row.names = F, file=DMR.filename)
  
  } else {dmr.ranges.anno <- "No DMRs found!"} # end of if condition 

  results[["DMR"]] <- dmr.ranges.anno
  
}

  # Detaching libraries not needed any more
  detach_package(unique(pks2detach))
  
return(results)

}






