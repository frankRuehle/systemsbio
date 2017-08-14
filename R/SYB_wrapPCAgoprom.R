
#' Principal component analysis for gene expression data 
#' 
#' Wrapper function for PCA routines from \code{pcaGoPromoter}-package incl. PCA-plots and 
#' enrichment analysis of PC loadings.
#' 
#' 
#' 2-dim and 3-dim PCA plots are generated for desired samples in the given ExpressionSet \code{expca}. 
#' Tables of PC-associated probes and transcription factor binding sites and GO terms enriched in top associated probes 
#' are generated for any number of principal components in positive and negative orientation.
#' All output data is stored in supplied \code{projectfolder}.
#' 
#' 
#' @param expca ExpressionSet object or a table with expression data with variables (probes) in rows and 
#'              observations in columns (samples).
#' @param groupsoi character vector with sample groups of interest to be included in PCA. 
#'           Respective samples are taken from \code{phenoData} of \code{expca}.
#'           groupnames must match entries in column given in \code{groupby}. 
#' @param groupby Column name of phenoData of \code{expca} used for group names.
#' @param sample.name.column Column name of phenoData of \code{expca} used for sample names
#' @param samples2exclude Character vector for optionally exclusion of individual samples. Used as 
#'                  regular expression for lookup of samples. \code{Null} if no sample to exlude.
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param inputType a character vector description of the input type. Must be Affymetrix chip type, 
#'            "geneSymbol" or "entrezID".
#' @param print.sample.names boolean indicating whether sample names shall be plotted in PCA plots 
#'                     (for pcainfoplot they are plotted anyway).  
#' @param org a character vector specifying the organism. Either "Hs" (homo sapiens), "Mm" (mus musculus) or "Rn" (rattus norwegicus).
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
#' The returned value is a list of 3 lists.
#' \itemize{
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
                          projectfolder= file.path("GEX", "pcaGoPromoter"), 
                          projectname=NULL, 
                          inputType="geneSymbol", 
                          print.sample.names = TRUE, 
                          org = "Hs", 
                          PCs4table = 2,  
                          PCs2plot = c(1,2,3), 
                          probes2enrich = 0.025 
) {
  

  # load required libraries
  pkg.bioc <- c("pcaGoPromoter", "pcaGoPromoter.Hs.hg19")
  pkg.cran <- c("rgl")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  ### Create result directory if not yet existing 
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
  if (!file.exists(file.path(projectfolder, "Gene_Ontology"))) {dir.create(file.path(projectfolder, "Gene_Ontology")) }
  if (!file.exists(file.path(projectfolder, "Transcription_factors"))) {dir.create(file.path(projectfolder, "Transcription_factors")) }
  if (!file.exists(file.path(projectfolder, "PC_loadings"))) {dir.create(file.path(projectfolder, "PC_loadings")) }
  
  
  
  ### Select Samples from groups of interest
  if(!is.null(groupsoi)) {
    
    groups.found <- groupsoi %in%  pData(expca)[,groupby]
    groups.included <- groupsoi[groups.found]
    
    cat("\n\nGroups included: ", paste(groups.included, collapse=", "))
    if(!all(groups.found)) {cat("\nGroup(s)", groupsoi[!groups.found], "not found in pheno data of expression set!")}
    
    groups.included.samples <- pData(expca)[,groupby] %in% groups.included # boolean. Samples of included groups
    expca <- expca[, groups.included.samples] # subset expression set for samples of included groups
       
  } else { # include all groups
    groups.included <- as.character(unique(pData(expca)[,groupby]))
    cat("\n\nAll groups are included:", paste(groups.included, collapse=", "))
  } # end of group selection
  
  
  
  if(!is.null(samples2exclude)) { # exclude samples matching reg expression in samples2exclude
    samples2exclude.found <- grepl(samples2exclude, pData(expca)[,sample.name.column])
    cat("\nSamples removed: ", paste(pData(expca)[samples2exclude.found,sample.name.column], collapse=", "))
    expca <- expca[, !samples2exclude.found] # subset expression set for samples2exclude
  }  
  
  cat("\nFinal sample list included: ", paste(pData(expca)[,sample.name.column], collapse=", "))
  # end of sample selection
  
  
  
  ### Row names should be probe identifiers given in inputType (Affymetrix probe set ID, "geneSymbol" or "entrezID")
  if (inputType=="geneSymbol") {
    if("SYMBOLREANNOTATED" %in% names(fData(expca))) {
      probe.identifier <- "SYMBOLREANNOTATED"
    } else { # grep first column containing 'symbol'
      probe.identifier <- grep("symbol", names(fData(expca)), value=T, ignore.case = TRUE)[1]
    }
    if(is.na(probe.identifier)) {stop("no gene symbols found in feature data!")}
    rownames(exprs(expca)) <- fData(expca)[,probe.identifier]
  }
  
  if (inputType=="entrezID") {
    if("ENTREZREANNOTATED" %in% names(fData(expca))) {
      probe.identifier <- "ENTREZREANNOTATED"
    } else { # grep first column containing 'entrez'
      probe.identifier <- grep("entrez", names(fData(expca)), value=T, ignore.case = TRUE)[1]
    }
    if(is.na(probe.identifier)) {stop("no entrez IDs found in feature data!")}
    rownames(exprs(expca)) <- fData(expca)[,probe.identifier]
  }
  
  
  
  ### Number of PC-associated probes to look for TFBS. Either fraction of total probe count or total number otherwise.
  noProbes <- if(probes2enrich <= 1) {round(probes2enrich*nrow(fData(expca)))} else {round(probes2enrich)}
  cat("\n", noProbes, "probes used for enrichment analysis.\n")  
  
  ### principal component analysis (pca)
  pcaOutput <- pcaGoPromoter::pca(exprs(expca))
  
  
  ### Make PCA informative plot. Restricted to PC1 and PC2 only!
  filename.pcaInfoPlot <- file.path(projectfolder, paste(projectname, "pcainfoplot_PC1_2.tiff", sep="_"))
  cat("\n\nSave pcaInfoPlot for PC1 and PC2 to", filename.pcaInfoPlot, "\n")
  tiff(filename= filename.pcaInfoPlot, width = 5500 , height = 5500, res=600, compression = "lzw") # width = 7016, height = 4960
  pcaInfoPlot(exprs(expca),inputType=inputType, org = org, 
              printNames = TRUE,
              groups=factor(pData(expca)[,groupby]), 
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
    
    filename.pcaplot <- file.path(projectfolder, paste0(projectname, "_pcaplot_PC", PCs[1], "_", PCs[2], ".tiff"))
    cat("\nSave pcaplot to", filename.pcaplot, "\n")
    tiff(filename=filename.pcaplot, width = 5500 , height = 5500, res=600, compression = "lzw") # width = 7016, height = 4960
    plot.pca(pcaOutput, groups=factor(pData(expca)[,groupby]), PCs = PCs, 
             printNames = print.sample.names, symbolColors = TRUE, plotCI = TRUE)
    dev.off()
  }   
  
  
  # 3D scatter plot with rgl package
  if (length(PCs2plot)>=3) {
    PCs3Dplot <- PCs2plot[1:3]
    # resize window
    par3d(windowRect = c(100, 100, 1000, 1000))
    filename.3dplot <- file.path(projectfolder, paste0(projectname, "_pcaplot3d_PC", paste(PCs3Dplot, collapse="_"), ".png"))
    cat("\n3D-plot generated for PCs", paste(PCs3Dplot, collapse="_"), "and stored at", filename.3dplot, ".\n")
    # plot groups
    plot3d(pcaOutput$scores[,PCs3Dplot], size=2, type = "s", col=as.numeric(factor(pData(expca)[,groupby])))
    grid3d("x")
    grid3d("y")
    grid3d("z")
    legend3d("topleft", legend=unique(as.character(factor(pData(expca)[,groupby]))), pch=16, cex=1.5, inset=c(0.02),
             col=unique(as.numeric(factor(pData(expca)[,groupby]))))
    # capture snapshot
    snapshot3d(filename = filename.3dplot, fmt = 'png')
  }
  
  
  
  
  
  ### Calculate singnificant TFBS and GO-Terms
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
                  file=file.path(projectfolder, "PC_loadings", paste0(projectname, "_Loading_", "PC", p, d, ".txt")))
      
      
      # Create Gene Ontology tree from loadings
      GOtreeOutput[[paste0("PC",p,d)]] <- GOtree(input = loadsperPC[[paste0("PC",p,d)]], 
                                                 inputType = inputType, org=org, 
                                                 statisticalTest = "binom", binomAlpha = NA,
                                                 p.adjust.method = "fdr")
      filename.GOtreeOutput <- file.path(projectfolder, "Gene_Ontology", paste0(projectname, "_GOtreeOutput_", "PC", p, d))
     
         
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
                  file= file.path(projectfolder, "Transcription_factors", paste0(projectname, "_overrepresentedTFs_", "PC", p, d,".txt")))
      
      write.table(TFtables[[paste0("PC",p,d)]][["underRepresented"]], row.names=F, quote=F, sep="\t",
            file= file.path(projectfolder, "Transcription_factors", paste0(projectname, "_underrepresentedTFs_", "PC", p, d,".txt")))


    } # end of d-loop
  } # end of p-loop
  
  
  # Detaching libraries not needed any more
   detach_package(c(pkg.cran, pkg.bioc))
  
  
  return(list(loadsperPC=loadsperPC, TFtables=TFtables, GOtreeOutput=GOtreeOutput))
  
} # end of function definition



