

#' Quality Control procedures for expression data
#' 
#' Application of diverse quality control metrix for gene expression data.
#' 
#' 
#' Function implements several approaches for quality control of expression data, which
#' are stored in a newly generated 'project folder.
#' - Illumina internal control plots, box-plots and MA-plots (if applicable).
#' - Quality control as implemented in the \code{arrayQualityMetrics}-package (subfolder 'quality_metrics')
#' - Quality control routines implemented in \code{WGCNA}-package. This comprises a sample dendrogram
#'   with included outlier prediction as well as boolean vectors for sample and gene quality
#'   calculated by the \code{goodSamplesGenes}-function.
#'   
#'   
#' @param eset ExpressionSet or ExpressionSetIllumina object  
#' @param projectfolder character vector with name of desired output folder (will be generated if not exisiting).
#' @param projectname character vector containing prefix for names of output files.
#' @param groupColumn character vector with name of group column in phenotype data.   
#' @param phDendro character vector with phenotypes of \code{eset} to be displayed in sample dendrogram. 
#' @param flashClustMethod character vector with the agglomeration method to be used in \code{WGCNA} package. 
#'                   Can be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".      
#' @param cex.dendroLabels numeric. Controls size of labels sample in dendrogram.
#' @param threshold.Z.k numeric. Threshold for outlier prediction \code{WGCNA} package (outlier samples are highlighted 
#'               in red in sample dendrogram).
#'               
#'               
#' @return The returned value is a list of 3
#' - goodGenes: boolean vector indicating genes passing QC-criteria as implemented in \code{WGCNA}.
#' - goodSamples: boolean vector indicating samples passing QC-criteria as implemented in \code{WGCNA}.
#' - allOK: boolean. If TRUE all samples and all genes pass QC-criteria as implemented in \code{WGCNA}.
#' Several plots and files are generated as side-effects and stored in the output folder. 
#' 
#' @author Frank Ruehle
#' 
#' @export QC_expressionset




  ## Usage 
  QC_expressionset <- function(eset, 
                     projectfolder= "GEX/QC",
                     projectname="", 
                     groupColumn = "Sample_Group",
                     phDendro = "Sample_Group", 
                     flashClustMethod = "average", 
                     cex.dendroLabels = 0.6,
                     threshold.Z.k=-5
                     ) {

    
   
  orig_par <- par(no.readonly=T)      # make a copy of current settings
  
  
  # Creating QC directory if not yet existing
  if (!file.exists(projectfolder)) {dir.create(projectfolder, recursive=T) }
                             

  
  # load required libraries
  pkg.cran <- c("WGCNA", "flashClust", "gridSVG")
  pkg.bioc <- c("arrayQualityMetrics", "beadarray", "Biobase", "BiocGenerics")
  attach_package(pkg.cran, pkg.bioc)
  

####### run arrayQualityMetrics.
cat("\nPrepare arrayQualityMetrics.")  
# arrayQualityMetrics works for ExpressionSet, AffyBatch, NChannelSet, ExpressionSetIllumina, RGList, MAList.
# how to fix error message "Only one 'gridsvg' device may be used at a time"???

  aqmdir <- file.path(projectfolder, "quality_metrics")
  if(file.exists(aqmdir)) {aqmdir <- paste0(aqmdir, "_new")}
  
try(
  aqMetrics <- arrayQualityMetrics(eset, outdir= aqmdir, 
                                     force=T, spatial=T, do.logtransform=F, intgroup=groupColumn,
                                     reporttitle = paste("arrayQualityMetrics_report"))
  ,silent=F)

closeAllConnections() # if arrayQualityMetrics produces error and leaves connection open

# #### graph par settings for Boxplots (do not work. Why?)
# par(cex.lab=cex.dendroLabels)
# par(cex.axis=cex.dendroLabels)

# Boxplot with regular probes
cat("\nprepare boxplot of probes in:", file.path(projectfolder, paste0(projectname, "boxplot_probes.pdf")), "\n")  
if("Status" %in% names(fData(eset))) {
  regularIDs <- featureNames(eset[which(fData(eset)[,"Status"] == "regular"),])
  } else {regularIDs <- 1:nrow(fData(eset))}

tiff(file=file.path(projectfolder, paste0(projectname,"boxplot_probes.tiff")), width = 7016 , height = 4960, res=600, compression = "lzw")
  BiocGenerics::boxplot(eset[regularIDs,])
dev.off()

# Boxplot with number of observations
cat("\nprepare boxplot of number of observations in:", file.path(projectfolder, paste0(projectname, "boxplot_nObservations.pdf")), "\n")  
tiff(file=file.path(projectfolder, paste0(projectname, "boxplot_nObservations.tiff")), width = 7016 , height = 4960, res=600, compression = "lzw")
  BiocGenerics::boxplot(eset, what="nObservations")
dev.off()



# Boxplot of with Illumina Control probes. ERCC-probes are bundled, if present.
# Illumina QC probes as exported from Genome Studio:
#   Biotin: Indicator for successful secondary staining (biotin labelled oligonucleotides) 
#   CY3_HYB: Cy3-labelled hybridisation controls (low, med and high concentration mixed!)
#   ERCC: 92 External RNA Controls Consortium spike-in control ratio mixtures in gene expression experiments (if applied).
#   HOUSEKEEPING: probes for housekeeping genes to monitor the intactness of the biological sample.
#   LABELING: optional, if used labeling > background, oterwise labeling ~ background
#   LOW_STRINGENCY_HYB: perfect matching control probes (pm) AND control probes containing 2 mismatches (mm2) 
#   NEGATIVE: random sequences
#   regular: probes for all regular genes

if(class(eset)=="ExpressionSetIllumina") {
  if("Status" %in% names(fData(eset))) {
      cat("\nprepare boxplot of Illumina control probes in:", file.path(projectfolder, paste0(projectname, "boxplot_controlprofile.pdf")), "\n")  
      ERCClevels <- grepl("ERCC", levels(fData(eset)$Status)  )  # fData are character, no factors
      if (any(ERCClevels)) {levels(fData(eset)$Status)[ERCClevels] <- "ERCC"}
      tiff(file=file.path(projectfolder, paste0(projectname, "boxplot_controlprofile.tiff")), width = 4960 , height = 7016, res=600, compression = "lzw")
    print(BiocGenerics::boxplot(eset, probeFactor = "Status", scales=list(cex.lab=cex.dendroLabels, cex.axis=cex.dendroLabels)))
    dev.off()
  }
}
  
# MA-Plots for all Samples separately (ERCC-probes are bundled, if present)
cat("\nprepare MA-plots of all samples in:", file.path(projectfolder, paste0(projectname, "MAplots.pdf")), "\n")  
class(eset) <- "ExpressionSet"  # plotMA not working with ExpressionSetIllumina- Object
pdf(file=file.path(projectfolder, paste0(projectname, "MAplots.pdf")))
for (i in 1:length(sampleNames(eset))) {
  limma::plotMA(eset, i, status=fData(eset)$Status)
  }
# remark: if no column "Status" exists in fData, then status=NULL, i.e. all points are plotted in the default color, symbol and size
dev.off()


# restore graphical parameter settings
par(orig_par) # restore settings





################ Clustering with WGCNA
cat("\n\nClustering with WGCNA\n")  
enableWGCNAThreads()

# working object with expression data:
netdat_noNorm <- as.data.frame(t(exprs(eset)[regularIDs,]))

# Implementing an additional color row within the dendrogram for outlier detection 
adjMatrix <- adjacency(t(netdat_noNorm), type="distance", distFnc = "dist")   # matrix re-transponed
# Calculation of distance calculated by dist()-function and method 'euclidean' within the adjacency function
# Distance Matrix = 1-adjMatrix
# The argument type determines whether a correlation (type one of "unsigned", "signed", "signed hybrid"), 
# or a distance network (type equal "distance") will be calculated.
# Correlation and distance are transformed as follows: 
# for type = "unsigned", adjacency = |cor|^power; 
# for type = "signed", adjacency = (0.5 * (1+cor) )^power; 
# for type = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise; 
# for type = "distance", adjacency = (1-(dist/max(dist))^2)^power. ## dist is normalised by max(dist)
# power = if (type=="distance") 1 else 6



# this calculates the whole network connectivity
k=as.numeric(apply(adjMatrix,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying if their Z.k value is below the threshold
# threshold.Z.k=-5 # often -2.5    # defined in function call
# the color vector indicates outlyingness (red)
outlierColor=ifelse(Z.k<threshold.Z.k,"red","black") # Outlier highlighted in red



# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers. 
# We use the function flashClust that provides faster hierarchical clustering than the standard function hclust.
# sampleTree_noNorm = flashClust(dist(netdat_noNorm), method = flashClustMethod) # calculated with dist()
#  Calculated from Adjacency Matrix (normalised by max(dist)).
sampleTree_noNorm_adj = flashClust(as.dist(1-adjMatrix), method = "average") 


### Loading Clinical trait data
traitData = pData(eset)
datTraits = as.data.frame(traitData[,phDendro])  # choose phenotypes for Dendrogram.
names(datTraits) <- phDendro

traitColors = labels2colors(datTraits); 
datColors=data.frame(outlier=outlierColor,traitColors)  # add line for Outlier detection

# Plot the sample dendrogram and the colors underneath.
cat("\nprepare sample dendrogram with phenotypes in:", file.path(projectfolder, paste0(projectname, "SampleDendrogram_noNorm_adjacency.pdf")), "\n")  
tiff(file.path(projectfolder, paste0(projectname, "SampleDendrogram_noNorm_adjacency.tiff")), width = 7016 , height = 4960, res=600, compression = "lzw")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree_noNorm_adj, colors=datColors, addGuide = T, guideAll=T,
                    groupLabels = c("outlier", names(datTraits)), cex.dendroLabels = cex.dendroLabels,
                    main = "Sample dendrogram and trait heatmap - GEX no normalisation")
dev.off()



cat("\nQuality Control with WGCNA: We check for genes and samples with too many missing values\n")
gsg = goodSamplesGenes(netdat_noNorm, verbose = 3);
cat(paste("\n\nWGCNA quality control. All ok:",gsg$allOK, "\n\n"))


  # Optionally, print the gene and sample names that should be removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Remove genes:", paste(names(netdat_noNorm)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Remove samples:", paste(rownames(netdat_noNorm)[!gsg$goodSamples], collapse = ", ")));


disableWGCNAThreads()

# restore graphical parameter settings
par(orig_par) # restore settings


# Detaching libraries not needed any more
# detach_package(c(pkg.cran, pkg.bioc))


return(gsg)

}


