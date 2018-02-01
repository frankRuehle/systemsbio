#' Quality control of methylation array data
#' 
#' Multiple QC procedures of an RGChannelSet object containing methylation data.
#' 
#' This function takes an RGChannelSet as input and transforms it to a MethylSet object using 3 preprocessing procedures 
#' from \code{minfi} package:
#' \itemize{
#'   \item \code{preprocessRaw} (no normalisiation)
#'   \item \code{preprocessIllumina} applies background subtraction as well as control normalization like in Illumina Genome studio
#'   \item \code{preprocessSWAN} applies SWAN normalization method on a preprocessed MethylSet to reduce the technical variability
#'                               of Type I and Type II Illumina probes. For this, Illumina preprocessed MethylSet is used if given,
#'                               the raw preprocessed MethylSet otherwise.  
#'   }
#' Quality control functions are applied to all 3 MethylSet objects. 
#' The required annotation packages (currently hg19 only) \code{IlluminaHumanMethylation450kmanifest} and 
#' \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} are installed and/or loaded automatically. 
#' Required packages not already attached when running \code{QC_methylation} are detached afterwards.
#' The function generates the following plots using quality control functionalities from up to 3 packages:
#' \itemize{
#'   \item \code{minfi}
#'     \itemize{
#'       \item qcReport using array control probes (for RGChannelSet)
#'       \item Median intensity plots
#'       \item Predicted sex plots
#'       \item Multi-dimensional scaling (MDS) plots
#'       \item Histogram overview of Illumina and SWAN normalisation
#'      }
#'   \item \code{COHCAP} (needs much computational resources)
#'    \itemize{
#'     \item Cluster plot
#'     \item PCA plot
#'     \item histogram plot
#'      }
#'   \item \code{WGCNA}
#'    \itemize{
#'     \item clustered sample dendrogram
#'     \item indicated outlier detection
#'      }
#'    }
#' 
#' 
#' @param RGset RGChannelSet object with methylation data.
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param sampleColumn character with column name of Sample names in \code{RGset}.
#' @param groupColumn character with column name of group names in \code{RGset}.
#' @param preprocessing character with preprocessing procedures to apply. Any of \code{"raw", "illumina", "swan"}.
#'                      If "swan" is selected, at least one more procedure must be selected to.  
#' @param QC_procedures character with desired QC procedures to use. Any of \code{"minfi", "COHCAP", "WGCNA"}.
#' @param numPositionsMDSplot number of genomic positions with the most methylation variability to be used for 
#'                    calculating distance between samples in MDS plots (needed for \code{minfi} only).
#' @param methPlatform Annotation file to be used. Enter "450k-UCSC" for UCSC CpG Islands for 450k array probes, "
#'                     450k-HMM" for HMM CpG Islands for 450k array probes, "27k" for UCSC CpG Islands for 27k array probes
#'                     (needed for \code{COHCAP} only).
#' @param phDendro character vector with phenotypes of \code{RGset} to be displayed in sample dendrogram 
#'                 (needed for \code{WGCNA} only).
#' @param flashClustMethod character vector with the agglomeration method to be used in \code{WGCNA} package. 
#'                   Can be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#'                   (needed for \code{WGCNA} only). 
#' @param cex.dendroLabels numeric. Controls size of labels sample in dendrogram (needed for \code{WGCNA} only).
#' @param threshold.Z.k numeric. Threshold for outlier prediction in \code{WGCNA} package. Outlier samples are highlighted 
#'               in red in sample dendrogram (needed for \code{WGCNA} only).
#' 
#' 
#' @return no value returned. QC plots are stored as side-effects.
#' 
#' @author Frank Ruehle
#' 
#' @export
 
 

###### Methylation Quality Control

QC_methylation <- function(RGset, 
                  projectfolder = "MT/QC",
                  projectname = NULL,
                  sampleColumn = "Sample_Name",
                  groupColumn = "Sample_Group",
                  preprocessing = c("raw", "illumina", "swan"),
                  QC_procedures = c("minfi", "WGCNA"),
                  numPositionsMDSplot= 1000, 
                  methPlatform = "450k-UCSC",
                  phDendro= "Sample_Group", 
                  flashClustMethod = "average", 
                  cex.dendroLabels = 0.6,  
                  threshold.Z.k=-5)
        {

  
  # load required packages. 
  annotation.package <- c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  pkg.bioc <- c("minfi")
  pkg.cran <- NULL 
  if ("COHCAP" %in% QC_procedures) {pkg.bioc <- c(pkg.bioc, "COHCAP")}
  if ("WGCNA" %in% QC_procedures)  {pkg.bioc <- c(pkg.bioc, "impute")
                                    pkg.cran <- c(pkg.cran, "WGCNA", "flashClust")}
  pks2detach <- attach_package(pkg.bioc = c(annotation.package, pkg.bioc))
  pks2detach <- c(pks2detach, attach_package(pkg.cran = pkg.cran)) # the bioc package impute must be loaded before WGCNA
  
   
  orig_par <- par(no.readonly=T)      # make a copy of current settings
   if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T)}
  
   pd <- pData(RGset) # phenotype data
   if (!is.null(projectname) && !grepl("_$", projectname)) {projectname <- paste0(projectname, "_")} 
  
  
  # select desired preprocessing procedures and determine reference object for SWAN normalisation if applicable
   objects.preprocess <-  c("MSet.raw", "MSet.normIllumina", "MSet.normSwan")
   objects.preprocess <- objects.preprocess[grepl( paste(preprocessing, collapse="|"), objects.preprocess, ignore.case = T)]
   if("MSet.normIllumina" %in% objects.preprocess) {object4Swan <- "MSet.normIllumina"} else {
       if ("MSet.raw" %in% objects.preprocess) {object4Swan <- "MSet.raw"} else {
         object4Swan <- NULL}
         }
     
  ###  Preprocessing (normalization) before QC (generating MethylSet Object from RGChannelSet Object)
  cat("\nPreprocessing and normalising RGChannelSet object to MethylSet object(s): ", objects.preprocess, "\n")
  if("MSet.raw" %in% objects.preprocess) {
  MSet.raw <- preprocessRaw(RGset)   #  simply converting the Red and the Green channel into a Methy-lated and Unmethylated signal
  }
  # minfi also allows for background subtraction (also called background normalization) as
  # well as control normalization like in Genome studio. Both of these are optional and turning
  # both of them turned off is equivalent to raw preprocessing (preprocessRaw).
  # control normalization requires the selection of one array among the 12 arrays on a chip as a reference array. 
  # It is currently unclear how Genome Studio selects the reference. Array 1 is selected arbitrarily.
  if("MSet.normIllumina" %in% objects.preprocess) {
    MSet.normIllumina <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 1) # as in GenomStudio.
  }
  # Technical differences have been demonstrated to exist between the Type I and Type II assay designs within a single 450K array.
  # Using the SWAN method substantially reduces the technical variability between the assay designs whilst maintaining the 
  # important biological differences. 
  if("MSet.normSwan" %in% objects.preprocess && !is.null(object4Swan)) {
    MSet.normSwan <- preprocessSWAN(RGset, get(object4Swan))   
  }
  # also available: 
  # preprocessQuantile() implements stratified quantile normalization preprocessing
  # preprocessFunnorm() implements the functional normalization algorithm developed in http://biorxiv.org/content/early/2014/02/23/002956.
  
  
  #### Quality control using minfi
  if ("minfi" %in% QC_procedures) {
    cat("\nQuality control with minfi\n")
    
  # Overview of all internal control probes
    # the values of the control probes are stored in the initial RGChannelSet
    filename.qcReport <- file.path(projectfolder, paste0(projectname, "qcReport_noNorm.pdf"))
    cat("\n  Writing Quality control report to", filename.qcReport, "\n")
    qcReport(RGset, sampNames = pd[,sampleColumn], sampGroups = pd[,groupColumn], pdf = filename.qcReport)
    # Overview control features:
    # Bisulfite conversion I: 3 probes high, 9 low
    # Bisulfite conversion II: 4 probes high red
    # Extension: 2 high, 2 low
    # Hybridization: 1 high, 1 med, 1 low (all green)
    # Non-Polymorphic: 2 high, 2 low
    # Specificity I: 3 probes high, 3 low (all green)
    # Specificity II: 3 probes high (red)
    # Target Removal: 2 probes low (green)
  
 
  # Median_Intensity_plots and predicted sex plots
    cat("\n  Creating Median intensity plots and predicted sex plots.\n")
    for (i in objects.preprocess) {
      qcMeth <- minfiQC(get(i), fixOutliers = TRUE, verbose = FALSE) # needs class 'MethylSet' or 'GenomicMethylSet'
      png(filename=file.path(projectfolder, paste0(projectname, "QC_Median_Intensity_plot_", i, ".png")), width = 210, height = 210, units = "mm", res=600)
      plotQC(qcMeth$qc, badSampleCutoff = 10.5)
      dev.off()

      # To predict the gender, minfi separate the points by using a cutoff of -2 on log2 med(X) - log2 med(Y).
      png(filename=file.path(projectfolder, paste0(projectname, "QC_predicted_sex_plot_", i, ".png")), width = 210, height = 210, units = "mm", res=600)
      plotSex(qcMeth$qc, id=row.names(qcMeth$qc))
      dev.off()
      
      gender.columns <- grep("(gender)|(sex)|(geschlecht)", names(pData(qcMeth$object)), ignore.case = T, value=T) # includes new column "predictedSex"
      if(length(gender.columns) >1) { # save tables with given and predicted gender if columns with given gender are found in phenotype data.
        write.table(pData(qcMeth$object)[, gender.columns, drop=F], row.names = T, quote = F, sep="\t",
                    file = file.path(projectfolder, paste0(projectname, "QC_predicted_sex_table_", i, ".txt")))
        }
      }
 
  # Multi-dimensional scaling (MDS) plots 
    cat("\n  Creating Multi-dimensional scaling (MDS) plots.\n")
    for (i in objects.preprocess) {
      pdf(file.path(projectfolder, paste0(projectname, "MDSplot_", i, ".pdf")), width = 10, height = 10) 
      minfi::mdsPlot(get(i), numPositions = numPositionsMDSplot, sampNames = pd[,sampleColumn], sampGroups = pd[,groupColumn],
                     main=paste0("MDS plot for Methylation data (", i,")"))
      dev.off()  ##### mdsPlot() my be masked by base
    }
      
  # Comparison of regular Illumina normalisiation with additional SWAN normalisation
      if("MSet.normSwan" %in% objects.preprocess && !is.null(object4Swan)) {
      cat("\n  Prepare comparison of", object4Swan, "with additional SWAN normalisation.\n")
    png(filename=file.path(projectfolder, paste0(projectname, "SWAN_normalisation.png")), width = 297, height = 210, units = "mm", res=600)
    par(mfrow=c(1,2))
    plotBetasByType(MSet.normIllumina[,1], main = object4Swan)
    plotBetasByType(MSet.normSwan[,1], main = "SWAN")
    par(mfrow=c(1,1))
    dev.off()
    }
  } # end minfi QC
  
  

  
  
  ###### QC with COHCAP
  if ("COHCAP" %in% QC_procedures) {
    cat("\nQuality control with COHCAP\n")

  if (!file.exists(file.path(projectfolder, "temp_data"))) {dir.create(file.path(projectfolder, "temp_data"))}
  
  # Sample file und beta-file have to be stored temporarily. 
  # The beta-file needs to be annotated with COHCAP.annotate().
  sampleSheetMethFile <- file.path(projectfolder, "temp_data", paste0(projectname, "tempSampleSheetMeth_COHCAP.txt")) # sampleSheet-file for CHOCAP
  # if (matchvarMeth=="none") {colsTargets <- c(sampleColumn,groupColumn)} else {colsTargets <- c(sampleColumn,groupColumn,matchvarMeth)}
  colsTargets <- c(sampleColumn, groupColumn)
  
  # REMARK: when run on UNIX system, COHCAP cannot load the sampleSheetMethFile 
  # unless saved with windows-like carriage return: eol = "\r\n"
  endofline <- if(Sys.info()["sysname"] == "Windows") {"\n"} else {"\r\n"}
  
  write.table(pData(RGset)[,colsTargets], file=sampleSheetMethFile, sep="\t", quote=F, row.names=F, col.names = F, eol = endofline)
  

  for (i in objects.preprocess) {
    
    betaForCOCAP.table <- data.frame(SiteID=rownames(get(i)), getBeta(get(i)))
    beta.file <- file.path(projectfolder, "temp_data",  paste0(projectname, "betaForCOCAP_", i, ".txt"))
    write.table(betaForCOCAP.table, file=beta.file, sep="\t", quote=F, row.names=F)
    # creates annotated beta table in "Raw_Data" folder
    beta.table <- COHCAP.annotate(beta.file, project.name= paste0(projectname,i), project.folder= projectfolder, platform=methPlatform, 
                                      annotation.file= NULL, output.format="txt")
    COHCAP.qc(sample.file=sampleSheetMethFile, beta.table=beta.table, 
              project.name = paste0(projectname,"COHCAP_", i), 
              project.folder = projectfolder) # COHCAP uses QC folder in projectfolder
      }

  } # end COHCAP QC
  

  
  
  #### Quality control with WGCNA
  if ("WGCNA" %in% QC_procedures) {
  
    # Use beta values for sample clustering:
    cat("\nClustering with WGCNA\n")

  enableWGCNAThreads()
  

       for (i in objects.preprocess) {
         
        MSet <- get(i)
        cat("\n  Creating Dendrogram for", i, "\n")
        netdatMT <- as.data.frame(t(getBeta(MSet)))

        # Implementing an additional color row within the dendrogram for outlier detection 
        adjMatrix <- adjacency(t(netdatMT), type="distance")   # matrix re-transponed
        # this calculates the whole network connectivity
        k=as.numeric(apply(adjMatrix,2,sum))-1
        # standardized connectivity
        Z.k=scale(k)
        # Designate samples as outlying if their Z.k value is below the threshold
        # threshold.Z.k=-5 # often -2.5    # defined in function call
        # the color vector indicates outlyingness (red)  
        outlierColor=ifelse(Z.k<threshold.Z.k,"red","black")   # Outlier highlighted in red
        
      
        # Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
        # outliers. We use the function flashClust that provides faster hierarchical clustering than the standard function hclust.
        # sampleTreeMT = flashClust(dist(netdatMT), method = flashClustMethod)  # calculated with dist()
        # Calculated from Adjacency Matrix (normalised by max(dist))
        sampleTreeMT_adj = flashClust(as.dist(1-adjMatrix), method = "average") 
 
        ### Loading Clinical trait data
        traitDataMT = pData(MSet.raw) # pData is identical for MSet.raw, MSet.normIllumina and MSet.normSwan
        datTraitsMT = as.data.frame(traitDataMT[,phDendro])  # choose phenotypes for Dendrogram.
        names(datTraitsMT) <- phDendro

        traitColorsMT = labels2colors(datTraitsMT) 
        datColors=data.frame(outlier=outlierColor,traitColorsMT)  # add line for Outlier detection
        
        
        # Plot the sample dendrogram and the colors underneath.
        png(filename=file.path(projectfolder, paste0(projectname, "SampleDendrogram_Betavalues_", i, ".png")), width = 297, height = 210, units = "mm", res=600)
        
        plotDendroAndColors(sampleTreeMT_adj, colors=datColors, 
                            groupLabels = c("outlier", names(datTraitsMT)),
                            main = paste("Sample dendrogram and trait heatmap -", i), 
                            cex.dendroLabels = cex.dendroLabels)
        dev.off()
        
    } # end of loop

  disableWGCNAThreads()
  
  } # end WGCNA QC
  
  par(orig_par) # resetting graphical parameter to original values
  
  # Detaching libraries not needed any more
  detach_package(unique(pks2detach))

}



