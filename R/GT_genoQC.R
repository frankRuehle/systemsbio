#' Quality Control for genotype data
#'
#' \code{genoQC} takes genotype data in GenABEL gwaa format and performs quality control and PCA analysis.
#' 
#' The check.marker-function from GenABEL package is used for quality control of genotype data.
#' It is recommended to perform two round of quality control: first QC, remove samples with different genetic
#' substructure, second QC. Principal component analysis for detection of genetic substructure is done
#' if \code{PCA} = TRUE. The first 10 principal components are added to the covariates of the \code{gwaa} object.
#' Samples are assigned to clusters and colored accordingly in PCA plots. Sample assignment is done for
#' up to \code{maxCenters} cluster centers. All cluster sample lists are stored in a subfolder "ClusterLists".
#' The QC-purified gwaa object may be exported to PLINK-compatible file formats.
#'
#' @param gwaa gwaa object from GenABEL
#' @param projectfolder character containing path to output folder (will be generated if not existing).
#' @param projectname character used as suffix for output files.
#' @param trait.name character indicating column name with trait of interest in pheno data of \code{gwaa}. 
#' Needed for qq-plots as well as for PCA plot annotation. Omitted if NULL.
#' @param trait.type character with data type "gaussian" or "binomial" of \code{trait.name}.
#' @param export.genofile character or character vector with type(s) of QC-purified ped file to export into \code{projectfolder}.
#' Allowed values are "ped" for ped/map-file, "tped" for transposed ped file format or "add.tped" for 
#' transposed additive coded format. If NULL, no data is exported. 
#' @param p.level.hwe numeric cut-off p-value for HWE in check.markers. For first round of QC it is rcommended to
#' skip p-level cut-off, i.e. set \code{p.level.hwe} = 0.
#' @param hwe.id.subset Subset for HWE checks in check.markers (default means controls only 
#' if \code{trait.name} is 0/1-coded affection status).
#' @param maf numeric cut off for minor allele frequency to be used in check.markers.
#' @param checkX boolean. If TRUE, X-errors in \code{gwaa} are fixed by Xfix().
#' @param PCA boolean. IF TRUE, PCA analysis performed with genotype data.
#' @param maxCenters numeric with maximum count of reported clustering center if PCA is performed.
#' @param ... further parameter submitted to GenABEL's check.marker() function. See \code{?check.marker} for details.
#' 
#' @return list containing two objects. First the QC-purified GenABEL gwaa object whith all samples  
#' removed as recomended. Second an object of class check.marker containing the quality control information.
#' Intermediary results and plots are stored in \code{projectfolder} as side effects.
#' 
#' @author Frank Ruehle
#' 
#' @export genoQC


genoQC <- function(gwaa,  
                   projectfolder = "GT/QC",
                   projectname = "QC1",
                   trait.name= "affection01", 
                   trait.type = "binomial",
                   export.genofile = "ped", 
                   p.level.hwe = 0.05,
                   hwe.id.subset = T, # phdata(gwaa)[,"affection01"]==0,  
                   maf = 0.01,   
                   checkX = T,   
                   PCA = T, 
                   maxCenters = 5,
                   ...)
                {


                      
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }  
    
   ## install/load required packages from CRAN and Bioconductor
   pkg.bioc <- NULL
   pkg.cran <- c("GenABEL")
   pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
                    
   projectname <- if (!is.null(projectname)) {paste0(projectname, "_")} else {""}
                    

  ## qqPlot before QC
  if(!is.null(trait.name)) {
     an1 <- qtscore(phdata(gwaa)[,trait.name], gwaa, trait= trait.type)
    lambda_beforeQC <- round(lambda(an1)$estimate, digits=4)
    cat(paste("\nlambda before QC =", lambda_beforeQC, "\n\n"))
    png(filename=file.path(projectfolder, paste0(projectname, "qqPlot_before_QC.png")), width = 210 , 
        height = 210, units = "mm", res=600)
      estlambda(an1[, "P1df"], plot=TRUE, main=paste0("qqPlot before QC (lambda = ", lambda_beforeQC, ")" ))
    dev.off()
  } else {cat("\nNo trait.name specified. qq-plot is omitted.")}

   
  # Quality Control with GenABELs check.marker()
  qc <- check.marker(gwaa, 
                    p.level= p.level.hwe, 
                    hweidsubset = hwe.id.subset,
                    maf= maf,   
                    ...
                    )  

    ## check.marker() default values
    # fdrate = 0.2, # cut-off FDR level in check for Hardy-Weinberg Equilibrium. Applied only when p.level negative!
    # callrate = 0.95,  # cut-off SNP call rate
    # perid.call=0.95,  # cut-off individual call rate (maximum percent of missing genotypes in a person)
    # extr.call = 0.1,  # SNPs with this low call rate are dropped prior to main analysis
    # extr.perid.call = 0.1,  # people with this low call rate are dropped prior to main analysis
    # het.fdr = 0.01,   # FDR rate for unacceptably high individual heterozygosity. outliers who have increased average heterozygosity may be suggestive of contaminated DNA samples.
    # ibs.threshold = 0.95,  # threshold value for acceptable IBS. E.g. people with exteremly high (close to 1) 
    #                           IBS may indicate duplicated samples (or twins), simply high values of IBS may indicate relatives. 
    # ibs.mrk = 2000,  # How many random markers should be used to estimate IBS
    # ibs.exclude="both",  # "both", "lower" or "none" - whether both samples with IBS>ibs.threshold should be excluded, the one with lower call rate, or no check
    # odds = 1000,   # cut-off odds to decide whether marker/person should be excluded based on sex/X-linked marker data inconsistency
    # minconcordance = 2.0,  #  If "minconcordance" is > 1.0 only pairs of markers which are exactly the same, including NA pattern, are considered as redundant


   
  
  cat("\nSummary of Quality Control\n")
  print(summary(qc))
  
  # Continue with good samples and snps
  gwaaqc <- gwaa[qc$idok, qc$snpok]
  
  qc_SNPsRemoved <- setdiff(snpnames(gwaa), qc$snpok)
  qc_SamplesRemoved <- setdiff(idnames(gwaa), qc$idok)
  
  write.table(qc_SNPsRemoved, file=file.path(projectfolder, paste0(projectname, "SNPs_removed_in_QC.txt")), row.names=F, quote=F, col.names="SNPs_Removed_in_QC")
  write.table(qc_SamplesRemoved, file=file.path(projectfolder, paste0(projectname, "Samples_removed_in_QC.txt")), row.names=F, quote=F, col.names="Samples_Removed_in_QC")
  
  
  # any residual sporadic X-errors (male heterozygosity) are set to NA
  if (checkX) {
    gwaaqc <- Xfix(gwaaqc)
    }
  
  
  ### Finding genetic sub-structure
  if(isTRUE(PCA)) {
    
    cat("\nPrepare Principal component analysis")
    gwaaqc.gkin <- ibs(gwaaqc[, autosomal(gwaaqc)], weight="freq")
    # The numbers below the diagonal show the genomic estimate of kinship (aka 'genomic kinship' or 'genome-wide IBD'), 
    # the numbers on the diagonal correspond to 0.5 plus the genomic homozygosity, and the numbers above the diagonal
    # tell how many SNPs were typed successfully for both subjects (thus the IBD estimate is derived using this number of SNPs).
    
  
  # transform this matrix to a distance matrix
  gwaaqc.dist <- as.dist(0.5-gwaaqc.gkin)
  
  # Multidimensional Scaling. Ten principal components are computed and and attached to covariates.
  gwaaqc.mds <- cmdscale(gwaaqc.dist, k=10)
  for (pc in 1:10) {gwaaqc <- add.phdata(gwaaqc, newph=gwaaqc.mds[,pc], name=paste0("PC",pc, ".", sub("_$", "", projectname))) }
  

  
  # PCA plot with first 2 PCs:
  cat(paste("\nWrite PCA plot to"), file.path(projectfolder, paste0(projectname, "PCAplot.png")))  
  png(filename=file.path(projectfolder, paste0(projectname, "PCAplot.png")), width = 210 , 
      height = 210, units = "mm", res=600)
      plot(gwaaqc.mds[,1:2], main=paste0(projectname, "PCA plot"), xlab="PC1", ylab="PC2", 
           pch= if(!is.null(trait.name) && trait.type=="binomial") {phdata(gwaaqc)[,trait.name]} else {1})
      text(gwaaqc.mds[,1:2], sub("@.*$", "", rownames(gwaaqc.mds)), pos=1, cex=0.4) # remove SNPZone suffix from sample names if necessary
      
       if(!is.null(trait.name) && trait.type=="binomial") {
        legend(x="topleft", legend=paste(trait.name, unique(na.omit(phdata(gwaaqc)[,trait.name]))), 
               pch=unique(na.omit(phdata(gwaaqc)[,trait.name])), cex=1, pt.cex=1)
      }
    dev.off()
  
  # PCA plots with samples highlighted corresponding to clusters calculated by first 2 PCs.
  # Corresponding cluster sample lists are stored in a subfolder "ClusterLists".
  if (!file.exists(file.path(projectfolder, "ClusterLists"))) {dir.create(file.path(projectfolder, "ClusterLists")) }
    cat(paste("\nAssigning samples to clusters. Cluster sample lists are stored in"), file.path(projectfolder, "ClusterLists"))
    
    for (i in 2:maxCenters) {
      cat(paste("\nassign samples to", i, "clusters"))
      km <- kmeans(gwaaqc.mds[,1:2], centers=i, nstart=1000)
      
      png(filename=file.path(projectfolder, paste0(projectname, "PCAplot_",i,"_center.png")), width = 210 , 
          height = 210, units = "mm", res=600)
        plot(gwaaqc.mds[,1:2], main=paste("PCAplot",i,"center"), xlab="PC1", ylab="PC2", col=km$cluster, 
             pch= if(!is.null(trait.name) && trait.type=="binomial") {phdata(gwaaqc)[,trait.name]} else {1})
        text(gwaaqc.mds[,1:2], sub("@.*$", "", rownames(gwaaqc.mds)), pos=1, cex=0.4, col=km$cluster) # remove SNPZone suffix from sample names if necessary
        legend(x="topright", legend=paste("cluster", 1:i), col=1:i, pch="-", cex=1, pt.cex=1)
        
        if(!is.null(trait.name) && trait.type=="binomial") {
          legend(x="topleft", legend=paste(trait.name, unique(na.omit(phdata(gwaaqc)[,trait.name]))), 
                 pch=unique(na.omit(phdata(gwaaqc)[,trait.name])), cex=1, pt.cex=1)
        }
      dev.off()
    
        for (j in 1:i) {
          cl <- names(which(km$cluster==j))
          write.table(cl, file=file.path(projectfolder, "ClusterLists", paste0(projectname, "ClusterList_Center_",i,"_Cluster_",j,".txt")), 
                  row.names=F, quote=F, col.names=paste0("Center_",i,"_Cluster_",j))
          }
     }
  }
  
 
  ## qqPlot after QC
  if(!is.null(trait.name)) {
    an2 <- qtscore(phdata(gwaaqc)[,trait.name], gwaaqc, trait= trait.type)
    lambda_afterQC <- round(lambda(an2)$estimate, digits=4)
    cat(paste("\n\nlambda of QC purified dataset =", lambda_afterQC))
    png(filename=file.path(projectfolder, paste0(projectname, "qqPlot_after_QC.png")), width = 210 , 
        height = 210, units = "mm", res=600)
      estlambda(an2[, "P1df"], plot=TRUE, main=paste0("qqPlot after QC (lambda = ", lambda_afterQC, ")" ))
    dev.off()
  }
  
  # Display summary data
  # cat("\nSample Summary Data\n")
  # print(idsum)
  cat("\n\nMarker Summary Data\n")
  print(descriptives.marker(gwaaqc))
  cat("\nTrait Summary Data\n")
  print(descriptives.trait(gwaaqc))
  
  
   
  # Export QC-purified ped-file. 
  # no phenotype is stored in ped file, column 6 is all zero. Phenotypes are stored in phe-file instead. 
  if(!is.null(export.genofile)) {

    for (g in 1:length(export.genofile)) {
      if(export.genofile[g] == "ped") {transpose = F; export012na = F; file.suffix = ".ped"}
      if(export.genofile[g] == "tped") {transpose = T; export012na = F; file.suffix = ".tped"}
      if(export.genofile[g] == "add.tped") {transpose = T; export012na = T; file.suffix = "_AddCoded.tped"}
      
      cat(paste("\nExporting QC-purified", export.genofile[g], "file to: ", file.path(projectfolder, paste0(projectname,"export", file.suffix)),"\n"))
      
      export.plink(gwaaqc, filebasename = file.path(projectfolder, 
                        paste0(projectname,"export", if(export.genofile[g] == "add.tped"){"_AddCoded"})),
                   phenotypes = "all", transpose = transpose,  export012na = export012na)
      }
  }
  
 
  return(list(gwaa.data = gwaaqc, QC = qc))


}
















