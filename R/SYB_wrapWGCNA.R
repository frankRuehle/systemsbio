
#' Wrapper for weighted gene co-expression network analysis.
#' 
#' This function makes use of the \code{WGCNA}-package from Steve Horvath and Peter Langfelder to construct 
#' weighted gene co-expression networks and correlates detected gene modules with phenotypes.
#' 
#' 
#' Before starting network construction, an appropriate \code{softThresholdPower} must be selected for correlation coefficients.  
#' If no value is given in \code{softThresholdPower}, the function analyses scale free topology for multiple soft 
#' thresholding powers to help choosing the appropriate value for obtaining an approximately scale free network topology.
#' For each power, the scale free topology fit index is calculated and returned along with other information on connectivity.
#' If \code{softThresholdPower} is set to 'auto' and the function determines an appropriate value and directly starts network construction.
#' Network construction is performed in block-wise manner with respect to \code{maxBlockSize}. Genes are clustered 
#' using average linkage hierarchical clustering and coexpressed gene modules are identified in the resulting dendrogram by the 
#' Dynamic Hybrid tree cut. Modules whose module eigengenes (MEs) are highly correlated are merged. 
#' The function calculates the following parameter:
#' \itemize{
#'   \item kME: INTRAmodular connectivity for finding intramodular hubs. Also known as module membership measure (MM). 
#'        Correlation of the gene with the corresponding module eigengene. kME close to 1 means that the gene is a hub gene.
#'   \item GS: gene significance: correlation of the gene with a phenotype.
#'   \item Module-trait relationship: correlation of a module eigengene with a phenotype.
#'   }
#' Phenotypes are taken from phenotype data of \code{data} as specified in \code{phModule}. Furthermore, membership of samples in groups
#' which are defined in \code{groupsets} are also used as phenotypes (e.g. two groups from a differential gene expression experiment). 
#' When correlation with group membership is calculated, only those samples are included which belong to the denoted groupset 
#' (mind that gene modules were calculated using expression data from all samples).
#' All correlation coefficients are calculated using Pearson correlation. Categorical variables with only two levels are 
#' coded numerically.
#' 
#' 
#' @param data ExpressionSet, SummarizedExperiment, DESeqDataSet or MethylSet. If \code{data} is character containing a filepath, the functions assumes
#'            previously stored network object to be loaded from this path. If \code{data} is "load_default", network 
#'            object is loaded from default directory \code{"file.path(projectfolder, "TOM", "networkConstruction-auto.RData")"}.
#' @param projectfolder character with directory for output files (will be generated if not existing).
#' @param softThresholdPower soft-thresholding power for network construction. If "auto", function selects 
#'                     soft-thresholding power automatically. If Null, network construction is omitted.
#' @param corType character string specifying the correlation to be used. Allowed values are "pearson" and "bicor", corresponding 
#'          to Pearson and biweight midcorrelation, respectively. Missing values are handled using the pairwise.complete.obs option.
#' @param networkType character with network type. Allowed values are "unsigned", "signed", "signed hybrid".
#'              "unsigned" means negative correlation of genes are treated the same as positive correlation.
#'              In an "signed" network, negatively correlated genes will not be put into one module, but will be treated as not correlated.
#' @param TOMType character with one of "none", "unsigned", "signed". If "none", adjacency will be used for clustering. 
#'           If "unsigned", the standard TOM will be used (more generally, TOM function will receive the adjacency as input). 
#'           If "signed", TOM will keep track of the sign of correlations between neighbours.
#' @param maxBlockSize integer giving maximum block size for module detection. If the number of genes in \code{data} exceeds \code{maxBlockSize}, 
#'               genes will be pre-clustered into blocks whose size should not exceed \code{maxBlockSize} (\code{maxBlockSize} must not exceed 46340).
#'               It's intended to use as big block sizes as possible, but mind that big blocksizes will heavily impact memory usage.
#' @param TOMplot boolean. If TRUE make Topological Overlap Matrix (TOM) plot (also known as connectivity plot) of the network connections.
#'          Light color represents low topological overlap and progressively darker red color represents higher overlap.
#'          Modules correspond to red squares along the diagonal.
#' @param MDSplot boolean. If TRUE make Multidimensional scaling plot (MDS) to visualize pairwise relationships specified 
#'          by a dissimilarity matrix. Each row of the dissimilarity matrix is visualized by a point in a Euclidean space.
#'          Each dot (gene) is colored by the module assignment.
#' @param phDendro character vector with phenotypes of \code{data} object to be displayed in sample dendrogram. 
#' @param phModule character vector with phenotypes to correlate module eigengenes with in heatmap.
#' @param sampleColumn character with column name of Sample names in pheno data of \code{data}.
#' @param groupColumn character with column name of group names in pheno data of \code{data}. 
#' @param groupsets character vector with names of group sets in format "groupA-groupB". Groups summarized in parentheses
#'            "(groupA-groupB)" are coded as ONE group. They are used for correlation of module eigengenes with 
#'            corresponding samples of selected groupsets. Mind that eigengenes are calculated using all samples,
#'            while correlation is calculated for samples of denoted groupsets only.
#'            Group names must match names in \code{groupColumn}. Omitted if NULL. 
#' @param symbolColumn character with name of feature identifier in feature data of \code{data}.
#' @param flashClustMethod character with agglomeration method used for hierarchical clustering in \code{flashClust}-package. 
#'                   Either "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". 
#' @param moduleBoxplotsPerFigure numeric. Number of module boxplots to be included in a single figure.
#' @param figure.res numeric resolution for png.
#' @param dendroRowText boolean. If TRUE, phenotype names are plotted beneath the sample dendrogram.
#' @param autoColorHeight boolean. If TRUE, the height of the color area below the dendrogram is adjusted 
#'                  automatically for the number of phenotypes.
#' @param colorHeight numeric specifying the height of the color area under dendrogram as a fraction of the height of the dendrogram area. 
#'              Only effective when autoColorHeight above is FALSE.
#' @param cex.labels numeric with character expansion factor for dendrogram and heatmap labels.
#' @param ... further arguments to be passed to the \code{blockwiseModules}-function of the \code{WGCNA}-package.
#' 
#' 
#' @return The returned value depends on parameter \code{softThresholdPower}. If a \code{softThresholdPower} is given or 
#' could be chosen automatically, value is a list with the following components:
#' \itemize{
#'   \item colors: a vector of color or numeric module labels for all genes
#'   \item unmergedColors: a vector of color or numeric module labels for all genes before module merging
#'   \item MEs: a data frame containing module eigengenes of the found modules (given by colors).
#'   \item goodSamples: numeric vector giving indices of good samples, that is samples that do not have too many missing entries.
#'   \item goodGenes: numeric vector giving indices of good genes, that is genes that do not have too many missing entries.
#'   \item dendrograms: a list whose components contain hierarchical clustering dendrograms of genes in each block.
#'   \item TOMFiles: character vector (one string per block), giving the file names in which blockwise topological overlaps were saved.
#'   \item blockGenes: a list whose components give the indices of genes in each block.
#'   \item blocks: a vector of length equal number of genes giving the block label for each gene. 
#'          Note that block labels are not necessarily sorted in the order in which the blocks were processed 
#'   \item MEsOK: logical indicating whether the module eigengenes were calculated without errors. 
#'   } 
#' 
#' If \code{softThresholdPower} is NULL or could not be chosen automatically, value is a list with the following components:
#' \itemize{
#'   \item powerEstimate: estimate of an appropriate soft-thresholding power: the lowest power for which the 
#'                 scale free topology fit R^2 exceeds RsquaredCut. If R^2 is below RsquaredCut for all powers, NA is returned.
#'   \item fitIndices: data frame containing the fit indices for scale free topology. The columns contain the 
#'              soft-thresholding power, adjusted R^2 for the linear fit, the linear coefficient, adjusted R^2 for a 
#'              more complicated fit models, mean connectivity, median connectivity and maximum connectivity. 
#'          }
#' \strong{Side-effects}: 
#' Diagrams for SoftThreshold power, gene and sample dendrograms generated by hierarchical clustering with phenotypes 
#' given in \code{phDendro} or \code{phModule} printed underneath as well as correlation heatmaps are plotted into the projectfolder. 
#' Additionally, tables with module eigengenes and correlation results of eigengenes with phenotypes and groupsets are generated.
#' Scatterplots are generated with module membership and gene significance for each phenotype/groupset and the 8 top associated modules.
#'
#'   
#' @references \code{https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/}
#' 
#' @author Frank Ruehle
#' 
#' @export
#' 
#' @note The procedure is divided in several steps:
#' \enumerate{
#'   \item Selection of an appropriate softThresholdPower for network construction
#'   \item Automatic network construction and module detection
#'   \item Plot sample dendrogram and gene dendrogramm with phenotype information.
#'         Calculate gene significance for traits: \code{GS.datTraits(i) = |cor(gene,Trait)|} and 
#'         \code{GSPvalue[i] = corPvalueStudent(GS.datTraits[i], nSamples)}.
#'   \item Correlation of modules with phenotypes (traits)
#'      \itemize{
#'         \item \code{moduleTraitCor =  cor(MEs, Trait)} 
#'         \item \code{moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)}
#'         \item \code{moduleGroupsetCor = cor(MEs, groupsetMat)}. Make \code{groupsetMat} as phenotype matrix 
#'         from group memberships for groups denoted in \code{groupsets}. 
#'         \item  \code{moduleGroupsetPvalue = corPvalueStudent(moduleGroupsetCor, nSamplesInGroups)}
#'         \item  Heatmaps are generated for correlation results of phenotypes and groupsets.
#'         }
#'   \item Intramodular analysis - Find hub genes in modules
#'      \itemize{
#'         \item datKME: INTRAmodular connectivity for finding intramodular hubs. Also known as module membership measure (MM).
#'         \item \code{MMPvalue = corPvalueStudent(datKME, nSamples)}  
#'         \item Calculate gene significance for group memberships: \code{geneGroupsetCor = cor(gene,groupsetMat)} and 
#'              \code{geneGroupsetPvalue = corPvalueStudent(geneGroupsetCor, nSamplesInGroups)}
#'         \item Scatterplots are generated with module membership and gene significance for each phenotype/groupset 
#'           and the 8 top associated modules.
#'           }
#'   \item Generate output tables
#'      \itemize{
#'         \item \code{networkDatOutput = data.frame(featuredata, moduleColors, GS.datTraits, GSPvalue, geneGroupsetCor, geneGroupsetPvalue)}
#'         \item \code{networkDatOutput_incl_MM = data.frame(networkDatOutput, datKME[,modOrder], MMPvalue[,modOrder])}
#'         }
#'   \item Visualization of networks within R (TOMplot, MDSplot).
#'         }


wrapWGCNA <- function(data, 
                  projectfolder= "GEX/WGCNA",
                  softThresholdPower="auto", 
                  corType="bicor", 
                  networkType = "signed", 
                  TOMType = "signed", 
                  maxBlockSize = 45000,
                  TOMplot=FALSE, MDSplot=FALSE,
                  phDendro=NULL, 
                  phModule=NULL, 
                  sampleColumn = "Sample_Name", 
                  groupColumn  = "Sample_Group",  
                  groupsets=NULL,  
                  symbolColumn = NULL, 
                  flashClustMethod = "average",
                  moduleBoxplotsPerFigure = 16, 
                  figure.res = 300,
                  dendroRowText=F, autoColorHeight = FALSE, 
                  colorHeight=0.1, cex.labels = 0.6,
                  ...
) {
  
  
  

  # load required libraries
  pkg.cran <- c("WGCNA", "ICC", "flashClust")
  pkg.bioc <- c("Biobase", "limma", "DESeq2", "SummarizedExperiment")
  pks2detach <- attach_package(pkg.cran, pkg.bioc)
  
  orig_par <- par(no.readonly=T)      # make a copy of current settings
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads()

  
  # Creating output directories if not yet existing
    if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
    if (!file.exists(file.path(projectfolder, "networkConstruction"))) {dir.create(file.path(projectfolder, "networkConstruction")) }
    if (!file.exists(file.path(projectfolder, "Intramodular_analysis_Traits"))) {dir.create(file.path(projectfolder, "Intramodular_analysis_Traits")) }
    #if (!file.exists(file.path(projectfolder, "VisANT"))) {dir.create(file.path(projectfolder, "VisANT")) }
    #if (!file.exists(file.path(projectfolder, "Cytoscape"))) {dir.create(file.path(projectfolder, "Cytoscape")) }
    
    if (!is.null(groupsets)) {
      if (!file.exists(file.path(projectfolder, "Intramodular_analysis_Groupsets"))) {dir.create(file.path(projectfolder, "Intramodular_analysis_Groupsets")) }
    }
    
  
  
    
    ### Check if input object 'data' is an ExpressionSet or MethySet object or character containing filepath to 
    # previously strored networl object.
    if(is.character(data)) {
      if (data=="load_default") {data <- file.path(projectfolder, "TOM", "networkConstruction-auto.RData")}
      load(file = data) # file.path(projectfolder, "TOM", "networkConstruction-auto.RData")
      } else {
    # Get expression/methylation data, and feature data from input object 
    # for class ExpressionSet
    if (grepl("Expression", class(data), ignore.case=T ) ) {
      cat("\nClass", class(data), "detected\n")
          if("Status" %in% colnames(fData(data))) { # define regular probes if control probes still present in expression set
            ids_regular <- rownames(fData(data)[fData(data)$Status=="regular",])
            } else {ids_regular <- 1:nrow(fData(data))}
      
      featuredata <- fData(data)[ids_regular,]
      featuredata <- data.frame(rowfeatures= rownames(featuredata), featuredata)
      assaydata <- exprs(data)[ids_regular,]
      plot.label="normalised log2(expression) data"
      plot.label.pruned="GEX"
      traitData = pData(data) # phenotype data 
      
           } else {
             
             if(class(data) %in% c("SummarizedExperiment", "DESeqDataSet", "DESeqTransform")) {
               cat("\nClass", class(data), "detected\n")
               featuredata <- as.data.frame(rowData(data,use.names=TRUE))
               featuredata <- data.frame(rowfeatures= rownames(featuredata), featuredata)
               assaydata <- assay(data)
               plot.label="normalised log2(count) data"
               plot.label.pruned="counts"
               traitData = colData(data) # phenotype data 
               
             } else {
        
        # for class MethylSet      
        if(grepl("Methyl", class(data), ignore.case=T ) ) {
          cat("\nClass", class(data), "detected\n")
          attach_package(pkg.bioc="minfi")
          featuredata <- mcols(data, use.names=T)
          featuredata <- data.frame(rowfeatures= rownames(featuredata), featuredata)
          assaydata <- getBeta(data)
          plot.label="normalised beta values"
          plot.label.pruned="MT"
          traitData = pData(data) # phenotype data 
          
        } else {cat("wrong object class!");return()}
      }
    }
  
  ## We transpose the expression data for further analysis.
  objectdat <- as.data.frame(t(assaydata))
  }
  
  # Define numbers of genes and samples
  nGenes = ncol(objectdat)
  nSamples = nrow(objectdat)
  
 
  
 
 
#### Step1: Selection of softThresholdPower for network construction
  if (is.null(softThresholdPower) || softThresholdPower=="auto") { # if no softThresholdPower selected yet

    # Horvath: CONSIDER ONLY THOSE PARAMETER VALUES IN THE ADJACENCY FUNCTION THAT RESULT IN APPROXIMATE SCALE FREE TOPOLOGY, 
    # i.e. high scale free topology fitting index R^2
    # we choose the lowest beta (softThresholdPower) that results in approximate scale free topology as measured by 
    # the scale free topology fitting index. In practice, we use the lowest value where the curve starts to "saturate"
    # Rationale: 
    #  - Empirical finding: Many co-expression networks based on expression data from a single tissue exhibit scale free topology
    #  - Many other networks e.g. protein-protein interaction networks have been found to exhibit scale free topology (SFT)
    # Caveat: When the data contains few very large modules, then the criterion may not apply. In this case, use the default choices.
    
    # Remark: "hard threshold" means unweighted network  
    
    # choose RsquaredCut
    RsquaredCut = 0.85
    
    # Call the network topology analysis function
    powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
    sft = WGCNA::pickSoftThreshold(objectdat, dataIsExpr = TRUE, RsquaredCut = RsquaredCut, 
                            powerVector = powers,
                            networkType = networkType, verbose =5)
  
    # Plot the results:
    cat("\nPlot SoftThreshold results to", file.path(projectfolder, "SoftThreshold.png"), "\n")
    png(file.path(projectfolder, "SoftThreshold.png"), width = 296, height = 148, units = "mm", res=figure.res)
    #pdf(file.path(projectfolder, "SoftThreshold.pdf"), width = 14, height = 7) 
    # sizeGrWindow(9, 5)
    par(mfrow = c(1,2))
    cex1 = 0.9
    
    # Scale-free topology fit index (SFT) as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1, col="red");
    
    # this line corresponds to using an R^2 cut-off of h   
    abline(h=RsquaredCut, col="red")
    
    # Gene connectivity = row sum of the adjacency matrix
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    
    dev.off()
    par(mfrow = c(1,1))
    
  
    # abbort function if no softThresholdPower is selected
    if (is.null(softThresholdPower)) {
      cat("\nSelect soft-thresholding power from plots and re-run function!\n")
      return(sft)}
    
    if (is.na(sft$powerEstimate)) {
      cat("\nNo softThresholdPower found with R2 >", RsquaredCut, ".")
      cat("\nSelect soft-thresholding power from plots and re-run function!\n")
      return(sft)}
    
    # Otherwise use softThresholdPower from sft object
    softThresholdPower <- sft$powerEstimate
    cat("\nSoftThresholdPower", softThresholdPower, "selected with R2 >", RsquaredCut, "\n")    

  } #  end selection of softThresholdPower




  
##### Step2: Automatic network construction and module detection (if soft-thresholding power is given)
  if(!is.character(data)) { # check if 'net' has already been loaded from file.
    cat(paste("\nConstructing the gene network and identifying modules; softThresholdPower =", softThresholdPower, "\n"))
  
    net <- WGCNA::blockwiseModules(objectdat, power = softThresholdPower,
                           networkType=networkType, TOMType = TOMType, 
                           corType=corType, 
                           numericLabels = FALSE, # modules labeled by colors (FALSE), or by numbers (TRUE)
                           pamStage = TRUE,  
                           saveTOMs = FALSE, saveTOMFileBase = file.path(projectfolder, "networkConstruction", "TOM"),
                           verbose = 1, maxBlockSize = maxBlockSize, ...
                           # deepSplit = 2, # default value
                           # detectCutHeight = 0.995, # default value
                           # minModuleSize = 20 # default: min(20, ncol(datExpr)/2 )
                           # mergeCutHeight = 0.15) # default value
                           ) 
  }
  
  # function blockwiseModules implements 3 steps:
  #   1) Variant of k-means to cluster variables into blocks
  #   2) constructs a correlation network, Hierarchical clustering and branch cutting in each block (Dynamic Hybrid tree cut).
  #   3) Merge modules across blocks (based on correlations between module eigengenes)
  
  # TOM = Topological overlap matrix 
  # a measure proportional to the number of neighbors that a pair of nodes share in common. 
  # The topological overlap measure can be interpreted as a measure of agreement between the m=1 step neighborhoods of 2 nodes.
  # Horvath generalized the topological overlap matrix to >=2 step neighborhoods
  
  # The argument networkType determines correlation (type one of "unsigned", "signed", "signed hybrid"): 
  # for type = "unsigned", adjacency = |cor|^power; 
  # for type = "signed", adjacency = (0.5 * (1+cor) )^power; 
  # for type = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise; 
  
  # number of blocks genes have been distributed to due to maxBlockSize
  blockCount <- length(na.omit(unique(net$blocks)))
  
  # Module Definition: branch of a cluster tree (see dendrogram). We use average linkage hierarchical clustering coupled 
  # with the topological overlap dissimilarity measure
  
  # Question: How does one summarize the expression profiles in a module?
  # Math answer: module eigengene (ME) = first principal component
  # Network answer: the most highly connected intramodular hub gene. Both turn out to be equivalent
  
  # module eigengenes
  MEs <- WGCNA::orderMEs(net$MEs)
  
  # To see how many modules were identified and what the module sizes are:
  cat("\n\nNetwork Modules\n")
  moduleColors = net$colors
  module.table <- as.data.frame(table(moduleColors))
  module.table <- module.table[match(gsub("ME", "", names(MEs)), module.table[,1]),] # same order of eigengenes as in MEs (needed for heatmap labels) 
  print(module.table)
  write.table(module.table, file.path(projectfolder, "Module_Sizes.txt"), quote=F, row.names=F, sep="\t")

  MEcount <- length(colnames(MEs))
  rownames(MEs) <- rownames(traitData)
  write.table(MEs, file.path(projectfolder, "ModuleEigengenes_colorLabel.txt"), row.names=F, quote=F, sep="\t")
  modNames = substring(names(MEs), 3)  # remove "ME" at the beginning of module eigengene names

  
  ## boxplots of modules vs. sample groups (16 module boxplots per figure)
  # prepare as many figures as necessary with respect to number of modules
  cat(paste("\nPreparing boxplots for each module with", moduleBoxplotsPerFigure, "plots per figure.\n"))
  countfigures <- ceiling(MEcount / moduleBoxplotsPerFigure)
  plotslastfig <- MEcount %% moduleBoxplotsPerFigure
  
  for (fig in 1:countfigures) { 
      
      png(filename=file.path(projectfolder, paste0("module_boxplot_", fig, "_of_", countfigures, ".png")), width = 210 , height = 210, units = "mm", res=figure.res)
      
      MEsfig <- MEs[, ((fig-1)*moduleBoxplotsPerFigure+1) : min(MEcount, (fig*moduleBoxplotsPerFigure))] # subset MEs for each figure
    
      layout(matrix(c(1:moduleBoxplotsPerFigure), ncol= floor(sqrt(moduleBoxplotsPerFigure)), byrow=TRUE))  #, heights=rep(3, times=length(colnames(MEs))))
      for (color in gsub("ME", "", colnames(MEsfig))){
        MEs_color <- MEsfig[, which(colnames(MEsfig) == paste0("ME", color)), drop = FALSE]
        traitsinfo <- cbind(MEs_color, factor(traitData[, groupColumn])) # dataframe with 1 ME and group assignments of all samples
        bp <- graphics::boxplot(traitsinfo[,1] ~ traitsinfo[,2], ylab = "module eigengene values", 
                                ylim = c(min(MEs_color[, 1]) - 0.2, max(MEs_color[, 1]) + 0.2), 
                                main = paste0(color, " module"), col = color, xaxt="n")
        
        axis(1, at=seq(1, length(bp$names), by=1), labels = FALSE)
        text(x=seq(1, length(bp$names), by=1), y=par("usr")[3] - 0.2, labels = bp$names, srt = 30, pos=1, xpd = TRUE) 
      }
      dev.off()
  }  
      
  
  
  
 
  # We now save the module assignment and module eigengene information necessary for subsequent analysis.
  save(net, objectdat, traitData, MEs, moduleColors, file = file.path(projectfolder, "networkConstruction", "networkConstruction-auto.RData"))
  
  
  
  
  
######## Step3: Make SAMPLE DENDROGRAM and GENE DENDROGRAM with correlatied phenotypes underneath.
  # Samples are clustered using 'flashClust' that provides faster hierarchical clustering than the standard function 'hclust'.
  datTraits.dendro = as.data.frame(traitData[,phDendro, drop=F])  # select phenotypes for Dendrogram
  sampleTree = flashClust(dist(objectdat), method = flashClustMethod)
  traitColors = labels2colors(datTraits.dendro)
  if (dendroRowText==F) {rowText=NULL} else {rowText=datTraits.dendro}
  
  ### plot sample dendrogram with phenotypes given in 'phDendro'
  png(file.path(projectfolder, "Sample_Dendrogram.png"), width = 297 , height = 210, units = "mm", res=figure.res)
  # pdf(file.path(projectfolder, "Sample_Dendrogram.pdf"), width = 12, height = 9) 
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  WGCNA::plotDendroAndColors(sampleTree, colors=traitColors, 
                      rowText=rowText, rowTextAlignment="left", cex.rowText= cex.labels/2, # cex.rowText throws error if labels are too large for plotting
                      autoColorHeight = autoColorHeight, colorHeight=colorHeight, addGuide=T, guideAll=T,
                      groupLabels = names(datTraits.dendro), 
                      cex.dendroLabels = cex.labels, cex.colorLabels = cex.labels, 
                      main = paste("Sample Dendrogram -", plot.label))
  dev.off()

  ##### GENE DENDROGRAM with correlation of genes with modules and traits 
  # Equivalent definitions of a gene significance measure
  #   GS.ClinicalTrait(i) = |cor(x(i),ClinicalTrait)| where x(i) is the gene expression profile of the i-th gene
  #   GS(i)=|T-test(i)| of differential expression between groups defined by the trait
  #   GS(i)=-log(p-value) GS.datTraits <- sapply(datTraits, cor, y=objectdat, use="p")
  # 
  # phenotypes need to be numeric to be correlated by Pearson. If factor variables have only 2 levels, they are 
  # converted to numeric by factor level. Other factor variables are correlated by Intraclass Correlation Coefficient (ICC)
  
  # select phenotypes
  datTraits.module = as.data.frame(traitData[,phModule, drop=F])  # select phenotypes for Dendrogram.
  names(datTraits.module) <- phModule

  ## Calculate correlation of genes with phenotypes
  # initialisation of data.frames
  datTraits <- datTraits.module
  GS.datTraits <- data.frame(row.names=colnames(objectdat))
  GSPvalue <- data.frame(row.names=colnames(objectdat))

  # loop(i) for every phenotype to calculate correlation depending on variable class
  for(i in 1:length(phModule)) {
      
      # check if any factor variable has 2 levels and can be converted to numeric (0/1 coded) 
       if(is.numeric(datTraits.module[,i])) {
          datTraits[,i] <- datTraits.module[,i]} else { 
            # if(length(levels(factor(datTraits.module[,i])))<=2) { ########### all categorical traits converted to numerici until ICC replaced by classifier!
              datTraits[,i] <- (as.numeric(factor(datTraits.module[,i])))-1
              cat("\n", colnames(datTraits)[i], "converted to numeric by factor level!\n")
              # } else { ########### all categorical traits converted to numerici until ICC replaced by classifier!
              #     datTraits[,i] <- factor(datTraits.module[,i]) ########### all categorical traits converted to numerici until ICC replaced by classifier!
              #    } ########### all categorical traits converted to numerici until ICC replaced by classifier!
          }
       
      # if phenotype is numeric variable, correlation with expression by pearson correlation
           # if(is.numeric(datTraits[,i])) { ##################### replace ICC by classifier
            GS.datTraits[i] <- as.data.frame(cor(objectdat, datTraits[i], use="pairwise.complete.obs"))
            names(GS.datTraits)[i] <- paste("cor",names(datTraits)[i], sep=".")
          # } else {
          #   # if phenotype is factor variable, correlation with expression by Intraclass Correlation Coefficient (ICC).
          #   # (Calculation is not vectorized, i.e. every value is concatenated by rbind)
          #   ICCgene <- numeric()
          #   for (j in 1:length(objectdat)) {
          #     ICCgene <- rbind(ICCgene, ICCest(datTraits[,i], objectdat[,j])$ICC)
          #                   }
          #   GS.datTraits[i] <- ICCgene
          #   names(GS.datTraits)[i] <- paste("icc",names(datTraits)[i], sep=".")
          #        }
        
      # Calculation of correlation p-value independantly of type of correlation
      GSPvalue[i] <- as.data.frame(corPvalueStudent(as.matrix(GS.datTraits[i]), nSamples))  # not used for Dendrogram
      colnames(GSPvalue)[i] <- paste("p", names(GS.datTraits)[i], sep=".")
  
  } # end of i-loop

  cat("\nPhenotypes are correlated with gene expression by Pearson correlation.", 
  # factor variables are correlated by Intraclass Correlation Coefficient (ICC). 
  "Factor variables with only two levels are converted to numeric.\n\n")
  print(data.frame(trait=phModule, class=sapply(datTraits[,,drop=F],class)))

  ###### plot dendrogram for each block:  
  # colors for traits plottedunderneath Dendrogram
  GS.TraitColor <- numbers2colors(GS.datTraits, signed=T)
  colnames(GS.TraitColor) <- colnames(GS.datTraits)
  # 'moduleColors' (first row underneath dendrogram) are module names.
  datColors <- data.frame(moduleColors, GS.TraitColor)  
  
  for (b in 1:blockCount) {
    png(filename=file.path(projectfolder, paste0("ModuleDendrogram_Block_", b, "_of_", blockCount,".png")), width = 297 , height = 210, units = "mm", res=figure.res)
    # pdf(file.path(projectfolder, paste0("ModuleDendrogram_Block_", b, "_of_", blockCount,".pdf")), width = 14, height = 10) 
    # Plot the dendrogram and the module colors underneath
    WGCNA::plotDendroAndColors(net$dendrograms[[b]], colors=datColors[net$blockGenes[[b]],],
                      groupLabels=colnames(datColors),
                      main=paste("Cluster Dendrogram - Block",b,"of",blockCount),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      cex.rowText = cex.labels, cex.dendroLabels = cex.labels, cex.colorLabels = cex.labels)
    dev.off()
    }

  

  


  
######### Step4: Correlation of modules with phenotypes and groupsets
# calculate moduleTraitCor as correlation of MEs with numeric traits by pearson correlation, 
# and for factor traits by Intraclass Correlation Coefficient (ICC).
# Calculate moduleTraitPvalue from moduleTraitCor (independently of type of correlation).

  ### Correlation of Eigengenes with traits (separated by trait type)
  # initialisation of data.frames
  moduleTraitCor <- data.frame(row.names=colnames(MEs))
  moduleTraitPvalue <- data.frame(row.names=colnames(MEs))
  # loop(i) for every phenotype to calculate correlation depending on variable class
  for(i in 1:length(phModule)) {
     
    # if phenotype is numeric variable, correlation with expression by pearson correlation
       # if(is.numeric(datTraits[,i])) { ##################### replace ICC by classifier
        moduleTraitCor[i] <- as.data.frame(cor(MEs, datTraits[i], use="pairwise.complete.obs"))
        names(moduleTraitCor)[i] <- paste("cor",names(datTraits)[i], sep=".")
      #   } else {
      #     # if phenotype is factor variable, correlation with expression by Intraclass Correlation Coefficient (ICC).
      #     # (Calculation is not vectorized, i.e. every value is concatenated by rbind)
      #     ICCmodule <- numeric()
      #     for (j in 1:length(MEs)) {
      #       ICCmodule <- rbind(ICCmodule, ICCest(factor(datTraits[,i]), MEs[,j])$ICC)
      #     }
      #     moduleTraitCor[i] <- ICCmodule
      #     names(moduleTraitCor)[i] <- paste("icc",names(datTraits)[i], sep=".")
      # }

      # Calculation of correlation p-value independently of type of correlation
    moduleTraitPvalue[i] <- corPvalueStudent(as.matrix(moduleTraitCor[i]), nSamples)  # not used for Dendrogram
    colnames(moduleTraitPvalue)[i] <- paste("p", colnames(moduleTraitCor)[i], sep=".")
  } # end of i-loop
  
  moduleTraitCor <- as.matrix(moduleTraitCor)

  #### HEATMAP Module-trait relationships
  # Will display correlations and their p-values
  cat("\nGenerating Heatmap_Module-trait_relationship.pdf\n")
  textMatrix = paste(signif(moduleTraitCor, 2), ", p=",
                     signif(as.matrix(moduleTraitPvalue), 1), sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  
  png(file.path(projectfolder, "Heatmap_Module-trait_relationship.png"), width = 210 , height = 297, units = "mm", res=figure.res)
  par(mar = c(6, 10, 4, 4));
  # Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(moduleTraitCor),
                 yLabels = names(MEs),
                 ySymbols = paste0(module.table[,1], " (", module.table[,2], ")"),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8 * cex.labels,
                 cex.lab = cex.labels,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships -", plot.label))
  dev.off()
  
  ### Correlation of Eigengenes with groupsets (separated by trait type)
  if(!is.null(groupsets)) {
    
    # initialisation of data.frames
    groupsetMat <- data.frame(row.names=rownames(traitData))
    samplecount <- numeric()
    
    # Matrix 'groupsetMat' generated with column for each groupset. NA assigned to samples not included in groupset.
    # Samples included in groupset coded by 0,1,2... in order of appearance in 'groupsets'.
    # Groups joined in parentheses are coded as ONE group.
    # When correlation calculated later with paramerter 'use = "pairwise.complete.obs"', all samples with NA are ignored. 
    for (i in 1:length(groupsets)) { # build groupsetMat (matrix with corresponding samples indicated)
      names(groupsets)[i] <- if(is.null(names(groupsets)[i]) || is.na(names(groupsets)[i]) || names(groupsets)[i] =="") {
                                 groupsets[i]} else {names(groupsets)[i]}
        
      groups2fuse <- gsub("[\\(\\)]", "", regmatches(groupsets[i], gregexpr("\\(.*?\\)", groupsets[i]))[[1]])
    
      if(length(groups2fuse) > 0) {
      for(j in 1:length(groups2fuse)) {groupsets[i] <- gsub(groups2fuse[j],gsub("-", "&", groups2fuse[j]), groupsets[i])} # replace "-" within parentheses by temp separator "&"
      groupsets[i] <- gsub("[\\(\\)]", "", groupsets[i]) # remove parentheses
      }
      
      groups2plot <- as.list(unlist(strsplit(groupsets[i], "-"))) # to have the strsplit elments as single list elements
      groups2plot <- sapply(groups2plot, strsplit, "&") # split the list elements in character vector if groups to fuse.
      
      sampleTable <- traitData[,c(sampleColumn,groupColumn)] # temporary table with sample and group names
      sampleTable$sets <- NA # initialise new column with NAs
      
      for(j in 1:length(groups2plot)) {
        sampleTable[sampleTable[,groupColumn] %in% groups2plot[[j]], "sets"] <- (j-1) # advise number to samples of involved groups
        }
      
      groupsetMat[,i] <- data.frame(sampleTable$sets) # build matrix with corresponding samples indicated and give names to columns
      names(groupsetMat)[i] <- names(groupsets)[i]
      
      samplecount[i] <- length(na.omit(groupsetMat[,i])) # necessary for moduleGroupsetPvalue
      
    } # end i-loop 
    
    cat("\nPhenotype matrix made from group sets: \n(stored in", file.path(projectfolder, "Intramodular_analysis_Groupsets", "Phenotype_matrix_Groupsets.txt"), ")\n")
    print(groupsetMat)
    write.table(groupsetMat, file=file.path(projectfolder, "Intramodular_analysis_Groupsets", "Phenotype_matrix_Groupsets.txt"), quote=F, sep="\t" )
    
    #### Correlation coefficient (type of correlation depending on number of groups)  
    # Pearson correlation ME - groupset (only samples without NA used)
      moduleGroupsetCor = cor(MEs, groupsetMat, use = "pairwise.complete.obs")  
      colnames(moduleGroupsetCor) <- names(groupsetMat)

    # # loop(i) for every phenotype to calculate ICC correlation if more then two groups
    # ##################### replace ICC by classifier
    # for(i in 1: ncol(groupsetMat)) {
    #   # if more then two groups, correlation with expression by Intraclass Correlation Coefficient (ICC).
    #   if(length(na.omit(unique(groupsetMat[,i]))) > 2) {
    #       samplesNA <- is.na(groupsetMat[,i]) # samples not to include in correlation (ICCest has't option "pairwise.complete.obs")
    #       ICCmodule <- numeric()
    #       for (j in 1:length(MEs)) {  # (Calculation is not vectorized, i.e. every value is concatenated by rbind)
    #         ICCmodule <- rbind(ICCmodule, ICCest(factor(groupsetMat[!samplesNA,i]), MEs[!samplesNA,j])$ICC)
    #         }
    #       moduleGroupsetCor[,i] <- ICCmodule
    #       colnames(moduleGroupsetCor)[i] <- paste("icc", colnames(moduleGroupsetCor)[i], sep=".")
    #     }
    # } # end of i-loop

    # prepare matrix with sample counts for corPvalueStudent
    samplecount.df <- data.frame(t(samplecount))
    samplecount.df <- samplecount.df[rep(1,each=nrow(moduleGroupsetCor)),] # copy row with sample counts for each module
    samplecount.df <- as.matrix(samplecount.df)
    
    moduleGroupsetPvalue = corPvalueStudent(moduleGroupsetCor, nSamples=samplecount.df)
    colnames(moduleGroupsetPvalue) <- colnames(moduleGroupsetCor)
    rownames(moduleGroupsetPvalue) <- rownames(moduleGroupsetCor)
    

    ### HEATMAP Module-GROUPSET relationships
    cat("\nGenerating Heatmap_Module-Groupset_relationship.pdf\n")
    textMatrix = paste(signif(moduleGroupsetCor, 2), ", p=",
                       signif(moduleGroupsetPvalue, 1), sep = "")
    dim(textMatrix) = dim(moduleGroupsetCor)
    
    png(file.path(projectfolder, "Heatmap_Module-Groupset_relationship.png"), width = 210 , height = 297, units="mm", res=figure.res)
    par(mar = c(6, 10, 4, 4));
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleGroupsetCor,
                   xLabels = gsub("-", " vs.\n", colnames(moduleGroupsetCor)),
                   yLabels = names(MEs),
                   ySymbols = paste0(module.table[,1], " (", module.table[,2], ")"),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.8 * cex.labels,
                   cex.lab = cex.labels,
                   zlim = c(-1,1),
                   main = paste("Module - Groupset relationships -", plot.label))
    dev.off()
 } # end of if (!is.null(groupsets)) 
    



############ Step5: Intramodular analysis: Find hub genes in modules
# Define two alternative measures of INTRAmodular connectivity for finding intramodular hubs:
# - Intramodular connectivity kIN: Row sum across genes inside a given module based on adjacency matrix
#   Disadvantage: strongly depends on module size
# - kME: Module eigengene based connectivity, also known as module membership measure: kME(i) = cor(x(i),ME)
#   kME(i) is simply the correlation between the i-th gene expression profile and the module eigengene.
#   kME close to 1 means that the gene is a hub gene

# Single network analysis: Intramodular hubs in biologically interesting modules are often very interesting
# Differential network analysis: Genes that are intramodular hubs in one condition but not in another are often very interesting

  # calculate the module membership values (aka. module eigengene based connectivity kME):
  datKME <- signedKME(objectdat, MEs)  # equals geneModuleMembership
  colnames(datKME) <- sub("kME", "MM.", colnames(datKME))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(datKME), nSamples));
  colnames(MMPvalue) <- sub("MM.", "p.MM.", colnames(MMPvalue))
  
  ######## Intramodular analysis: Phenotypes
  # identifying genes with high GS (Gene Significance) and MM (Module Membership)
  # (Gene significance ('GS.datTraits' and 'GSPvalue') has already been calculated in step3: gene dendrogram)
  # Using the GS and MM measures, we can identify genes that have a high significance for a phenotype
  # as well as high module membership in interesting modules
  cat("\nGenerating plots for Intramodular analysis (phenotypes)\n")
  colorOfColumn <- substring(names(datKME),4)

  for (trait in phModule) {
    modTraitName <- grep(paste0(trait,"$"), names(moduleTraitPvalue), value=T) # with prefix "p.cor" or "p.icc"
    # select 8 top associated modules for each trait given in 'phModule': 
      cat(paste("\ntrait:", trait, "\n"))
      selectModules <- rownames(moduleTraitPvalue[order(moduleTraitPvalue[,modTraitName], decreasing=F), , drop=F])[1:min(8, ncol(MEs))]
      selectModules <- substring(selectModules,3) # remove substring "ME"

      # plot 8 scatter plots for each trait
      png(file.path(projectfolder, "Intramodular_analysis_Traits", paste0("Intramodular_analysis_", trait, ".png")), width = 210 , height = 297, units ="mm", res=figure.res)
      par(mfrow=c(length(selectModules)/2,2))
      par(mar=c(6, 8, 4, 4) + 0.1)
        
        for (module in selectModules) {
          column <- match(module,colorOfColumn)
          restModule <- moduleColors==module
          WGCNA::verboseScatterplot(datKME[restModule,column],GS.datTraits[restModule, grep(paste0(trait,"$"),names(GS.datTraits), value=T)],
                       xlab=paste("Module Membership:",module,"module"),ylab=paste("Gene significance:", trait),
                       main=paste(plot.label.pruned, "Module membership vs. gene significance\n"), 
                       cex = cex.labels, cex.axis = cex.labels, cex.lab = cex.labels, cex.main = 1.3 * cex.labels,
                       pch=21, col="black", bg=module)
        }
      dev.off()
  } # end of trait-loop
  
  par(orig_par) # resetting graphical parameter to original values

  ######## Intramodular analysis: groupsets
  if(!is.null(groupsets)) {
  # geneGroupsetCor: Correlate genes in objectdat with groupsetMat (analogous to GS.datTraits above)
  # geneGroupsetPvalue: Calculate p-vale from geneGroupsetCor (analogous toGSPvalue)
  # groupsetMat and samplecount have been computed earlier in step4
  
  #### Correlation coefficient (type of correlation depending on number of groups)  
  # Pearson correlation gene - groupset (only samples without NA used)
  geneGroupsetCor <- cor(objectdat, groupsetMat, use="pairwise.complete.obs")   
  colnames(geneGroupsetCor) <- paste0("cor.",colnames(geneGroupsetCor))
  
  # # If necessary, pearson correlation coefficients are overwritten by ICC
  # ##################### replace ICC by classifier
  # # loop(i) for every groupset to calculate ICC correlation if more then two groups
  # for(i in 1: ncol(groupsetMat)) {
  #   # if more then two groups, correlation with expression by Intraclass Correlation Coefficient (ICC).
  #   if(length(na.omit(unique(groupsetMat[,i]))) > 2) {
  #     samplesNA <- is.na(groupsetMat[,i]) # samples not to include in correlation (ICCest has't option "pairwise.complete.obs")
  #     ICCgene <- numeric()
  #     for (j in 1:length(objectdat)) {  # (Calculation is not vectorized, i.e. every value is concatenated by rbind)
  #       ICCgene <- rbind(ICCgene, ICCest(factor(groupsetMat[!samplesNA,i]), objectdat[!samplesNA,j])$ICC)
  #     }
  #     geneGroupsetCor[,i] <- ICCgene
  #     colnames(geneGroupsetCor)[i] <- paste("icc", colnames(geneGroupsetCor)[i], sep=".")
  #   }
  # } # end of i-loop

  
  # prepare matrix with sample counts for corPvalueStudent(). Calculation is independent from correlation type
  samplecount.df.group <- data.frame(t(samplecount))
  samplecount.df.group <- samplecount.df.group[rep(1,each=nrow(geneGroupsetCor)),] # copy row with sample counts for each module
  samplecount.df.group <- as.matrix(samplecount.df.group)
  
  geneGroupsetPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneGroupsetCor), nSamples=samplecount.df.group))  
  colnames(geneGroupsetPvalue) <- sub("cor.", "p.cor.", colnames(geneGroupsetCor))
  rownames(geneGroupsetPvalue) <- rownames(geneGroupsetCor)

  # plot scatter plots for 8 top-associated modules
  cat("\nGenerating plots for Intramodular analysis (group sets)\n")
  colorOfColumn <- substring(names(datKME),4)
  
  for (gset in names(groupsets)) {
    # select 8 top associated modules for each Groupset given in 'groupsets': 
    cat(paste("\nGroupset:", gset, "\n"))
    moduleGroupsetPvalue <- as.data.frame(moduleGroupsetPvalue)
    selectModules <- rownames(moduleGroupsetPvalue[order(moduleGroupsetPvalue[,gset], decreasing=F), , drop=F])[1:min(8, ncol(MEs))] # drop=F to avoid loss of dimensions when just 1 Groupset
    selectModules <- substring(selectModules,3) # remove substring "ME"
      
    # plot 8 scatter plots for each group Groupset
    png(file.path(projectfolder, "Intramodular_analysis_Groupsets", paste0("Intramodular_analysis_", gset, ".png")), 
         width = 210 , height = 297, units="mm", res=figure.res)
    par(mfrow=c(length(selectModules)/2,2))
    par(mar=c(6, 8, 4, 4) + 0.1)
    
    for (module in selectModules) {
      column <- match(module,colorOfColumn)
      restModule <- moduleColors==module  # logical vector with length of number of genes. TRUE for genes of selected module
      WGCNA::verboseScatterplot(datKME[restModule,column],geneGroupsetCor[restModule, grep(gset,colnames(geneGroupsetCor), value=T)],
                         xlab=paste("Module Membership:",module,"module"),ylab=paste("Gene significance:", gset),
                         main=paste(plot.label.pruned, "Module membership vs. gene significance\n"), 
                         cex = cex.labels, cex.axis = cex.labels, cex.lab = cex.labels, cex.main = 1.3 * cex.labels,
                         pch=21, col="black", bg=module)
    }
    dev.off()
  } # end of gset-loop
  par(orig_par) # resetting graphical parameter to original values
  
  } # end of if (!is.null(groupsets)) 
  ####################### end of intramodular analysis





  ################## Step6: Create output tables
  # 'networkDatOutput' contains feature data of all genes, module assignment, Gene correlation with traits
  # and (if selected) gene correlation with groupsets
  #
  # featuredata = feature annotations
  # moduleColors = Modul assignment
  # GS.datTraits equals geneTraitSignificance. prefix: "cor" or "icc"
  # GSPvalue = p_value(geneTraitSignificance). prefix: "p.cor" or "p.icc"
  # datKME equals geneModuleMembership (eigengene-based connectivity, also known as module membership, "MM")
  # MMPvalue = p_value(geneModuleMembership). "p.MM"
  
  cat("\nGenerating Network Output Tables\n")
  networkDatOutput0 <- data.frame(featuredata, moduleColors, GS.datTraits, GSPvalue)
  
  if(!is.null(groupsets)) {
    networkDatOutput0 <- data.frame(networkDatOutput0, geneGroupsetCor, geneGroupsetPvalue, check.names = F)
  } 
  
  # sorting rows and columns if groupColumn exists
  groupColumn.name <- grep(paste0("^p.*", groupColumn),names(networkDatOutput0), value=T) # get name of pvalue-column from cor with groups
  if(length(groupColumn.name)==1) {
    # sort sequence of module-membership-columns with respect to the correlation of modules with groupColumn before adding to networkDatOutput
    # (if correlation between modules and groupColumn is already done in 'moduleTraitCor')
       modOrder = order(-abs(moduleTraitPvalue[,groupColumn.name])) 
        # Add module membership information in 'datKME' and MMPvalue in the chosen order
       networkDatOutput_MM <- data.frame(networkDatOutput0, datKME[,modOrder], MMPvalue[,modOrder], check.names = F)
       
       # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance (if column with pvalue from cor with groupColumn exists)
       geneOrder = order(networkDatOutput0$moduleColors, -abs(networkDatOutput0[,groupColumn.name]))
       networkDatOutput0 = networkDatOutput0[geneOrder, ]
       networkDatOutput_MM  = networkDatOutput_MM[geneOrder, ]
       }
  
  write.table(networkDatOutput0,file.path(projectfolder, "networkDatOutput.txt"), row.names=F, quote=F, sep="\t")
  write.table(networkDatOutput_MM,file.path(projectfolder, "networkDatOutput_incl_MM.txt"), row.names=F, quote=F, sep="\t")

  ### saving pruned result tables: the 8 most significant modules for each trait
  for (trait in phModule) {
      # select 8 top associated modules for each trait given in 'phModule': 
      modTraitName <- grep(paste0(trait,"$"), names(moduleTraitPvalue), value=T) # with prefix "p.cor" or "p.icc"
      selectModules <- rownames(moduleTraitPvalue[order(moduleTraitPvalue[,modTraitName], decreasing=F), , drop=F])[1:min(8, ncol(MEs))]
      selectModules <- substring(selectModules,3) # remove substring "ME"
      i <- 0      
       for (module in selectModules) {
         i <- i+1
         restModule <- networkDatOutput_MM$moduleColors==module
         defColnames <- c(colnames(featuredata), "moduleColors", 
                          grep(paste0(trait,"$"), names(GS.datTraits), value=T), 
                          grep(paste0(trait,"$"), names(GSPvalue), value=T))  
         moduleColnames <- grep(paste0("MM.",module, "$"), names(networkDatOutput_MM), value=T)
         
         Output.pruned <- networkDatOutput_MM[restModule, c(defColnames, moduleColnames)]
         write.table(Output.pruned,file.path(projectfolder, "Intramodular_analysis_Traits", paste0(paste("out",trait,"no",i,"module",module, sep="_"),".txt")), 
                     row.names=F, quote=F, sep="\t")
       } # end of module-loop
  } # end of trait-loop 

  ### saving pruned result tables: the 8 most significant modules for each groupset
  if(!is.null(groupsets)) {
    for (gsets in names(groupsets)) {
    # select 8 top associated modules for each trait given in 'phModule': 
    selectModules <- rownames(moduleGroupsetPvalue[order(moduleGroupsetPvalue[,gsets], decreasing=F), , drop=F])[1:min(8, ncol(MEs))]
    selectModules <- substring(selectModules,3) # remove substring "ME"
    
    i <- 0      
    for (module in selectModules) {
      i <- i+1
      restModule <- networkDatOutput_MM$moduleColors==module
      defColnames <- c(colnames(featuredata), "moduleColors", 
                       grep(gsets, colnames(geneGroupsetCor), value=T), 
                       grep(gsets, colnames(geneGroupsetPvalue), value=T))  
      moduleColnames <- grep(paste0("MM.",module, "$"), names(networkDatOutput_MM), value=T)
      
      Output.pruned <- networkDatOutput_MM[restModule, c(defColnames, moduleColnames)]
      write.table(Output.pruned,file.path(projectfolder, "Intramodular_analysis_Groupsets", paste0(paste("out",gsets,"no",i,"module",module, sep="_"),".txt")), 
                  row.names=F, quote=F, sep="\t")
     } # end of module-loop
    } # end of gsets-loop 
  } # end of if(!is.null(groupsets))



  
  
  
  
##########################################
# Systems genetic approach for integrating gene co-expression network analysis with
# genetic data  (Ghazalpour A, et al 2006,Presson A, et al 2008).

# A module quantitative trait locus (modulQTL) is a chromosomal location (e.g. a SNP)
# which correlates with many gene expression profiles of a given module.

# The following steps summarize our overall approach:
#   (1) Construct a gene co-expression network from gene expression data 
#   (2) Study the functional enrichment (gene ontology etc) of network modules (e.g. with DAVID)
#   (3) Relate modules to the clinical traits 
#   (4) Identify chromosomal locations (QTLs) regulating modules and traits 
#       (use a correlation test if the SNP is numerically (additive) coded)
#   (5) Use QTLs as causal anchors in network edge orienting analysis to find causal drivers underlying the clinical traits.
#       we use genetic markers to calculate Local Edge Orienting (LEO) scores for 
#         a) determining whether the module eigengene has a causal effect on the phenotype, and 
#         b) for identifying genes inside the trait related modules that have a causal effect on the phenotype.




####### Step7: Visualization of networks within R                   

### Topological Overlap Matrix (TOM) plot (also known as connectivity plot) of the network connections.
# Light color represents low topological overlap and progressively darker red color represents higher overlap.
# Modules correspond to red squares along the diagonal.
if (TOMplot | MDSplot) {

    cat("\nGenerating Topological Overlap Matrix (TOM) Plot\n")

    # Calculate topological overlap new: this could be done more efficiently by saving the TOM
    # calculated during module detection, but let us do it again here.
    # TOM <- 1-TOMsimilarityFromExpr(objectdat, power = softThresholdPower, corType=corType,
    #                                   networkType=networkType, TOMType = TOMType)
  TOMnew <- 1-TOMsimilarityFromExpr(objectdat, power = 5, corType="bicor", networkType="signed", TOMType = "signed")
    
  ###### load TOM from network
     for (tomfile in net$TOMFiles) {load(tomfile)}
    # wie kann ich TOMs aus verschiedenen Bl?cken mergen??

    # Distance Matrix
    dissTOM <- 1-TOM
}


if(TOMplot) {
    # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
    plotTOM <- dissTOM^7;
    # Set diagonal to NA for a nicer plot
    # diag(plotTOM) <- NA; # Error: only matrix diagonals can be replaced
    # Call the plot function: Topological Overlap Matrix
    geneTree <- flashClust(as.dist(dissTOM), method="average")
    png(file.path(projectfolder, paste("Network_heatmap_plot.png")), width = 210 , height = 297, units="mm", res=figure.res)
    #pdf(file.path(projectfolder, paste("Network_heatmap_plot.pdf")), width = 10, height = 14) 
    WGCNA::TOMplot(dissim=plotTOM, dendro=geneTree, Colors=moduleColors, main = "Network heatmap plot, all genes")                 
    dev.off()
}


### Multidimensional scaling (MDS) 
# visualizing pairwise relationships specified by a dissimilarity matrix.
# Each row of the dissimilarity matrix is visualized by a point in a Euclidean space.
# Each dot (gene) is colored by the module assignment.
if(MDSplot) {
  cat("\nGenerating Multidimensional scaling (MDS) Plot\n")

  cmd1 <- cmdscale(as.dist(dissTOM), k=2)  # classical MDS plot using 2 scaling dimensions.
  png(file.path(projectfolder, paste("WGCNA_MDS_plot.png")), width = 210 , height = 297, units="mm", res=figure.res)
  #pdf(file.path(projectfolder, paste("WGCNA_MDS_plot.pdf")), width = 10, height = 14) 
  plot(cmd1, col=moduleColors, main="WGCNA MDS plot", xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
  dev.off()
}




if(FALSE) {
###### VisANT 

# software framework for visualizing, mining, analyzing and modeling multi-scale biological networks. 
cat("\nGenerating input files for VisANT\n")

# for WGCNA::exportNetworkToVisANT. Plot Symbol in network if available, else plot probename
probeToGene <- data.frame(probes=rownames(featuredata), 
                          Symbols= ifelse(featuredata[,symbolColumn] != "" & !is.na(featuredata[,symbolColumn]), 
                                          featuredata[,symbolColumn], rownames(featuredata)) )


for (trait in phModule) {
  # select 8 top associated modules for each trait given in 'phModule': 
  selectModules <- rownames(moduleTraitPvalue[order(moduleTraitPvalue[,trait], decreasing=F),])[1:8]
  selectModules <- substring(selectModules,3) # remove substring "ME"
    
  index <- 0
  for (module in selectModules) {
    index <- index+1
    # 8 network exports for each trait
    
    probes = names(objectdat)
    inModule <- moduleColors==module
    modProbes = probes[inModule]
    
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule]
    # Because the module is rather large, we focus on the 30 top intramodular hub genes
    nTopHubs = 30
    
    # intramodular connectivity
    kIN = softConnectivity(objectdat[, modProbes])
    selectHubs = (rank (-kIN) <= nTopHubs)
    vis = exportNetworkToVisANT(modTOM[selectHubs,selectHubs],
                file=file.path(projectfolder, "VisANT", paste0("VisANT_input_", paste(trait, index, module, sep="_"), ".txt")),
                weighted=TRUE, threshold = 0, 
                probeToGene=probeToGene)
        } # end of module loop
   } # end of trait-loop




###### Cytoscape  
cat("\nGenerating input files for Cytoscape\n")

for (trait in phModule) {
  # select 8 top associated modules for each trait given in 'phModule': 
  selectModules <- rownames(moduleTraitPvalue[order(moduleTraitPvalue[,trait], decreasing=F),])[1:8]
  selectModules <- substring(selectModules,3) # remove substring "ME"
  
  index <- 0
  for (module in selectModules) {
    index <- index+1
    # 8 network exports for each trait
    
    probes = names(objectdat)
    
    # Select module probes
    inModule <-is.finite(match(moduleColors,module))
    modProbes <- probes[inModule]
    
    match1=match(modProbes,probeToGene$probes)
    modGenes=probeToGene$Symbols[match1]
    
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule]
    dimnames(modTOM) = list(modProbes, modProbes)
    
    # Export the network into edge and node list files for Cytoscape
    cyt = exportNetworkToCytoscape(modTOM,
                 edgeFile=file.path(projectfolder, "Cytoscape", paste0("CytoEdge_",paste(trait, index, module, sep="_"),".txt")),
                 nodeFile=file.path(projectfolder, "Cytoscape", paste0("CytoNode_",paste(trait, index, module, sep="_"),".txt")),
                 weighted = TRUE, threshold = 0.02, nodeNames=modProbes,
                 altNodeNames = modGenes, nodeAttr = moduleColors[inModule])
    
      } # end of module loop
    } # end of trait-loop


}







# resetting options to default
options(stringsAsFactors = TRUE)
disableWGCNAThreads()
par(orig_par) # resetting graphical parameter to original values

# Detaching libraries not needed any more
detach_package(unique(pks2detach))

return(net)

}
