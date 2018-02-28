
#' Principal component analysis for gene expression data 
#' 
#' Wrapper function for PCA routines from \code{pcaGoPromoter}-package incl. PCA-plots and 
#' enrichment analysis of PC loadings.
#' 
#' The pcaGoPromoter::pca function uses prcomp to do the principal component analysis. 
#' The input data is scaled and centered, so constant variables (sd = 0) will be removed to avoid divison by zero.
#' 2-dim and 3-dim PCA plots are generated for desired samples in the given ExpressionSet \code{expca}. 
#' Tables of PC-associated probes and transcription factor binding sites and GO terms enriched in top correlated probes 
#' are generated for any number of principal components in positive and negative orientation.
#' All output data is stored in supplied \code{projectfolder}.
#' 
#' 
#' @param expca \code{ExpressionSet} object or a table with expression data with variables (probes) in rows and 
#'              observations in columns (samples). In latter case, rows of data matrix must be named after probe 
#'              identifiers selected in \code{inputType}.
#' @param groupsoi character vector with sample groups of interest to be included in PCA (if \code{expca} is an \code{ExpressionSet}). 
#'           Respective samples are taken from \code{phenoData} of \code{expca}.
#'           groupnames must match entries in column given in \code{groupby}. 
#' @param groupby character with column name of phenoData of \code{expca} used for group names if \code{expca} is 
#' an \code{ExpressionSet}. Otherwise, \code{groupby} must be vector of group assignments in the same order as samples in the data matrix.
#' @param sample.name.column Character with column name of phenoData of \code{expca} used for sample names
#' @param samples2exclude Character vector for optionally exclusion of individual samples. Used as 
#'                  regular expression for lookup of samples. \code{Null} if no sample to exlude.
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param figure.res numeric resolution for png.
#' @param inputType Character vector with description of the input type. Must be Affymetrix chip type, 
#'            "geneSymbol" or "entrezID".
#' @param print.sample.names boolean indicating whether sample names shall be plotted in PCA plots 
#'            (for pcainfoplot they are plotted anyway).  
#' @param print.symbol.colors boolean indicating whether the symbols should be plotted with colors.
#' @param org a character vector specifying the organism. Either "Hs" (homo sapiens), "Mm" (mus musculus) 
#' or "Rn" (rattus norwegicus).
#' @param annotation.packages character with bioconductor annotation packages to load. 
#' E.g. c("pcaGoPromoter.Hs.hg19", "org.Hs.eg.db") for human or c("pcaGoPromoter.Mm.mm9", "org.Mm.eg.db") for mouse. 
#' @param PCs4table numeric or numeric vector. Indicates number of PCs (numeric) or distinct PCs (numeric vector) 
#'            for which result tables of enriched transcription factor binding sites and GO-terms are calculated.
#' @param PCs2plot numeric or numeric vector. Indicates number of PCs (numeric) or distinct PCs (numeric vector) 
#'           to use in 2-dim and 3-dim PCA plots. For 2-dim PCA plots all possible pairs of PCs are plotted.
#'           Additionally, a 3D plot is generated with the first 3 PCs in PCs2plot.
#'           Note that pca informative plot (containing TFBS and GO annotation on the axes) 
#'           is restricted to first two PCs only! 
#' @param probes2enrich numeric. Number of top PC-associated probes to look for enriched TFBS and GO terms. 
#'              A value \code{<= 1} is interpreted as fraction of total number of probes.
#'              
#'              
#' @return  Several plots and files are generated as side-effects and stored are in the designated projectfolder.
#' The returned value is a list of 4 objects.
#' \itemize{
#'    \item PCA: Principal component matrix
#'    \item loadsperPC: Top associated probes for every PC in pos and neg direction
#'    \item TFtables: dataframes containing enriched TFBS for every PC in pos and neg direction (over- and underrepresented)
#'    \item GOtreeOutput: dataframes containing enriched GO terms for every PC in pos and neg direction
#'   }    
#' 
#' @author Frank Ruehle
#' 
#' @export

  
wrapPCAgoprom <- function(expca, 
                          groupsoi = NULL, 
                          groupby ="Sample_Group", 
                          sample.name.column = "Sample_Name",
                          samples2exclude = NULL,
                          projectfolder= file.path("pcaGoPromoter"), 
                          projectname=NULL, 
                          figure.res = 300,
                          inputType="geneSymbol", 
                          print.sample.names = TRUE,
                          print.symbol.colors = TRUE,
                          org = "Hs", 
                          annotation.packages = c("pcaGoPromoter.Hs.hg19", "org.Hs.eg.db"),
                          PCs4table = 2,  
                          PCs2plot = c(1,2,3), 
                          probes2enrich = 0.025
) {
  

  # load required libraries
  pkg.bioc <- c("pcaGoPromoter", "GO.db", annotation.packages)
  pkg.cran <- c("rgl")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  ### Create result directory if not yet existing 
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
  if (!file.exists(file.path(projectfolder, "Gene_Ontology"))) {dir.create(file.path(projectfolder, "Gene_Ontology")) }
  if (!file.exists(file.path(projectfolder, "Transcription_factors"))) {dir.create(file.path(projectfolder, "Transcription_factors")) }
  if (!file.exists(file.path(projectfolder, "PC_loadings"))) {dir.create(file.path(projectfolder, "PC_loadings")) }
  
  if(!is.null(projectname)) {projectname <- paste0(projectname, "_")}
  

  
  
  
 if(is.matrix(expca)) { # make expressionSet if expca is matrix
   cat("\nRows of data matrix must be named after probe identifiers selected in inputType 
                    (Affymetrix probe set ID, geneSymbol or entrezID)\n")
   cat("groupby must be vector of group assignments in the same order as samples in data matrix")
   expca <- ExpressionSet(assayData = expca,
                        phenoData = data.frame(Sample_Name=colnames(expca), Sample_Group = factor(groupby)), 
                        experimentData = MIAME(),
                        annotation = character())
   fData(expca) <- data.frame(features= rownames(expca))
   }
  
 
  
  if(class(expca) %in% c("ExpressionSetIllumina", "ExpressionSet")) {
 
       if(!is.null(groupsoi)) {
         groups.found <- groupsoi %in%  pData(expca)[,groupby] 
         groups.included <- groupsoi[groups.found]
         groups.included.samples <- pData(expca)[,groupby] %in% groups.included # boolean. Samples of included groups
         expca <- expca[, groups.included.samples] # subset expression set for samples of included groups
       }
         
        if(!is.null(samples2exclude)) { # exclude samples matching reg expression in samples2exclude
           samples2exclude.found <- grepl(samples2exclude, pData(expca)[,sample.name.column])
           expca <- expca[, !samples2exclude.found] # subset expression set for samples2exclude
         }
   
      groupvector <- factor(pData(expca)[,groupby]) 
      samplevector <- pData(expca)[,sample.name.column]
      expressionMatrix <- exprs(expca)
      featureData <- fData(expca)
         
   } else {
    
      if(class(expca) %in% c("SummarizedExperiment", "DESeqDataSet")) {
 
        if(!is.null(groupsoi)) {
          groups.found <- groupsoi %in%  colData(expca)[,groupby] 
          groups.included <- groupsoi[groups.found]
          groups.included.samples <- colData(expca)[,groupby] %in% groups.included # boolean. Samples of included groups
          expca <- expca[, groups.included.samples] # subset expression set for samples of included groups
        }
        
        if(!is.null(samples2exclude)) { # exclude samples matching reg expression in samples2exclude
          samples2exclude.found <- grepl(samples2exclude, colData(expca)[,sample.name.column])
          expca <- expca[, !samples2exclude.found] # subset expression set for samples2exclude
        }

        groupvector <- factor(colData(expca)[,groupby]) 
        samplevector <- colData(expca)[,sample.name.column]
        expressionMatrix <- assay(expca)
        featureData <- as.data.frame(rowData(expca,use.names=TRUE))
        rownames(featureData) <- rownames(assay(expca))
        } 
  }
  
  cat("\nSample list included: ", paste(samplevector, collapse=", "))

  
  
 
  
  
   
  
  
  
  ### Row names should be probe identifiers given in inputType (Affymetrix probe set ID, "geneSymbol" or "entrezID")
  # for ExpressionSet objects, the annotation can be obtained from the feature data table

        if (inputType=="geneSymbol") {
            if("SYMBOLREANNOTATED" %in% names(featureData)) {
              probe.identifier <- "SYMBOLREANNOTATED"
                } else { # grep first column containing 'symbol'
                  probe.identifier <- grep("symbol", names(featureData), value=T, ignore.case = TRUE)[1]
                }
            if(!is.na(probe.identifier)) {
              cat("\nUsing column", probe.identifier, "as gene Symbol.")
              rownames(expressionMatrix) <- featureData[,probe.identifier]
            } else {
              cat("\nno gene symbols found in feature data. Using rownames instead") # default anyway
              }
      }
      
      if (inputType=="entrezID") {
          if("ENTREZREANNOTATED" %in% names(featureData)) {
            probe.identifier <- "ENTREZREANNOTATED"
            } else { # grep first column containing 'entrez'
              probe.identifier <- grep("entrez", names(featureData), value=T, ignore.case = TRUE)[1]
            }
          if(!is.na(probe.identifier)) {
            cat("\nUsing column", probe.identifier, "as EntrezID.")
            rownames(expressionMatrix) <- featureData[,probe.identifier]
          } else {
            cat("no entrez IDs found in feature data. Using rownames instead")
            }
      }
  
  
  
  ### Number of PC-associated probes to look for TFBS. Either fraction of total probe count or total number otherwise.
  probesInDataset <- nrow(featureData)
  noProbes <- if(probes2enrich <= 1) {round(probes2enrich*probesInDataset)} else {round(probes2enrich)}
  cat("\n", noProbes, "probes used for enrichment analysis.\n")  
  
  
  
  
  
  ### principal component analysis (pca)
  pcaOutput <- pcaGoPromoter::pca(expressionMatrix)
  

  ### Make PCA informative plot. Restricted to PC1 and PC2 only!
  filename.pcaInfoPlot <- file.path(projectfolder, paste0(projectname, "pcainfoplot_PC1_2.png"))
  cat("\n\nSave pcaInfoPlot for PC1 and PC2 to", filename.pcaInfoPlot, "\n")
  png(filename= filename.pcaInfoPlot, width = 150, height = 150, units = "mm", res=figure.res) 
  #tiff(filename= filename.pcaInfoPlot, width = 5500 , height = 5500, res=600, compression = "lzw") # width = 7016, height = 4960
    pcaInfoPlot(expressionMatrix,inputType=inputType, org = org, 
              printNames = TRUE,
              groups= groupvector, 
              noProbes = noProbes, GOtermsAnnotation = TRUE, primoAnnotation = TRUE)
  dev.off()
  

  ### PCA plot not restricted to PC1 and PC2 (but still 2D)
  # plot.pca makes same pca-plot as pcaInfoPlot() but without TFBS annotation
  if(length(PCs2plot)==1) {PCs2plot <- 1:PCs2plot}
  
  PC.combinations <- combn(PCs2plot,2) # get all possible PC pairs from PCs2plot
  cat("Generate 2D PCA plots for", ncol(PC.combinations), "PC pair(s).\n")
  if(ncol(PC.combinations)>10) {
    cat("\nWarning: ", ncol(PC.combinations), "combinations of PCs possible. Plots are restricted to first 10 PC pairs.\n")
    PC.combinations <- PC.combinations[,1:10]}
  
  for (i in 1:ncol(PC.combinations)) {  
    PCs <- PC.combinations[,i]
    
    filename.pcaplot <- file.path(projectfolder, paste0(projectname, "pcaplot_PC", PCs[1], "_", PCs[2], ".png"))
    cat("\nSave pcaplot to", filename.pcaplot, "\n")
    png(filename= filename.pcaplot, width = 150, height = 150, units = "mm", res=figure.res) 
    #tiff(filename=filename.pcaplot, width = 5500 , height = 5500, res=600, compression = "lzw") # width = 7016, height = 4960
    pcaGoPromoter::plot.pca(pcaOutput, groups= groupvector, PCs = PCs, 
             printNames = print.sample.names, symbolColors = print.symbol.colors, plotCI = TRUE)
    dev.off()
  }   
  

  
  # 3D scatter plot with rgl package
  if (length(PCs2plot)>=3) {
    PCs3Dplot <- PCs2plot[1:3]
    # resize window
    par3d(windowRect = c(100, 100, 1000, 1000))
    filename.3dplot <- file.path(projectfolder, paste0(projectname, "pcaplot3d_PC", paste(PCs3Dplot, collapse="_"), ".png"))
    cat("\n3D-plot generated for PCs", paste(PCs3Dplot, collapse="_"), "and stored at", filename.3dplot, ".\n")
    # plot groups
    rgl::plot3d(pcaOutput$scores[,PCs3Dplot], size=2, type = "s", col=as.numeric(groupvector))
    grid3d("x")
    grid3d("y")
    grid3d("z")
    rgl::legend3d("topleft", legend=unique(as.character(groupvector)), pch=16, cex=1.5, inset=c(0.02),
             col=unique(as.numeric(groupvector)))
    # capture snapshot
    rgl::snapshot3d(filename = filename.3dplot, fmt = 'png')
  }
  

  
  
  
  ### Calculate significant TFBS and GO-Terms
  cat("\nCalculate enriched transcription factor binding sites and GO-Terms for first", noProbes, "probes of each PC.")
  cat("\nSave GOtreeOutput plots and tables to", file.path(projectfolder, "Gene_Ontology"), ".")
  cat("\nSave over-representated TF tables to", file.path(projectfolder, "Transcription_factors"), ".")
  
  if(length(PCs4table)==1) {PCs4table <- 1:PCs4table}
  
  # initialisation for loop
  loadsperPC <- list()
  GOtreeOutput <- list()
  TFtables <- list()
  for (p in PCs4table) {
    
    for (d in c("pos", "neg")){
      cat("\nProcess PC", p, d, "direction\n")
      
      # Get loadings (probes) from PCs. We would like to have the top probe ids from every PC in pos and neg direction.
      loadsperPC[[paste0("PC",p,d)]] <- getRankedProbeIds(pcaOutput, pc=p, decreasing=(d=="pos") )[1:noProbes]
      # Hansen 2012: "The function getRankedProbeIds generates a ranked list of the probe set IDs that mostly 
      # contribute for placing experiments along the chosen principal component"
      
      # save loadings
      write.table(loadsperPC[[paste0("PC",p,d)]], row.names=F, quote=F, sep="\t", col.names = F,
                  file=file.path(projectfolder, "PC_loadings", paste0(projectname, "Loading_", "PC", p, d, ".txt")))
      
      
      # Create Gene Ontology tree from loadings
      GOtreeOutput[[paste0("PC",p,d)]] <- GOtree(input = loadsperPC[[paste0("PC",p,d)]], 
                                                 inputType = inputType, org=org, 
                                                 statisticalTest = "binom", binomAlpha = NA,
                                                 p.adjust.method = "fdr")
      filename.GOtreeOutput <- file.path(projectfolder, "Gene_Ontology", paste0(projectname, "GOtreeOutput_", "PC", p, d))
     
         
#       # plot GO term trees # currently not working: Error in addNode("n0", g) : argument "to" is missing, with no default  
#       GO2plot <- GOtreeOutput[[paste0("PC",p,d)]]
#       if(nrow(GO2plot[['sigGOs']]) >=2) { # plot only if at least 2 GO-Terms
#           tiff(filename=paste0(filename.GOtreeOutput, ".tiff"), width = 4960 , height = 7016, res=600, compression = "lzw")
#           plot.GOtree(GO2plot, legendPosition = NULL, boxes = min(20, nrow(GO2plot[['sigGOs']])))
#           dev.off() # for unknown reason plot is not colored according to legend. Therefore legend omitted.
#         }
        
      # save significant GO-Terms as table
        sigGOTerms <- GOtreeOutput[[paste0("PC",p,d)]][["sigGOs"]]
        write.table(sigGOTerms[sigGOTerms$pValue < 0.05,], row.names=F, quote=F, sep="\t",
                    file=paste0(filename.GOtreeOutput, ".txt"))
        
      
      
      # Get list of enriched transcription factors
      TFtables[[paste0("PC",p,d)]] <- primo(loadsperPC[[paste0("PC",p,d)]], inputType = inputType, org=org,
                                            PvalueCutOff = 0.05, cutOff = 0.9,
                                            p.adjust.method = "fdr", printIgnored = FALSE , primoData = NULL)
      write.table(TFtables[[paste0("PC",p,d)]][["overRepresented"]], row.names=F, quote=F, sep="\t",
                  file= file.path(projectfolder, "Transcription_factors", paste0(projectname, "overrepresentedTFs_", "PC", p, d,".txt")))
      
      write.table(TFtables[[paste0("PC",p,d)]][["underRepresented"]], row.names=F, quote=F, sep="\t",
            file= file.path(projectfolder, "Transcription_factors", paste0(projectname, "underrepresentedTFs_", "PC", p, d,".txt")))


    } # end of d-loop
  } # end of p-loop
  
  
  # Detaching libraries not needed any more
   detach_package(c(pkg.cran, pkg.bioc))
  
  
  return(list(PCA=pcaOutput, loadsperPC=loadsperPC, TFtables=TFtables, GOtreeOutput=GOtreeOutput))
  
} # end of function definition



