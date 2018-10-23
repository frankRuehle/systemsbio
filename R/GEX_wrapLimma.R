
#' Differential expression analysis with limma
#' 
#' Differential Gene expression analysis for multiple group comparisons incl. output of 
#' Venn diagrams, volcano plots and heatmaps.
#' 
#' Function generates design matrix, fit and differential expression results for dedicated 
#' group comparisons of the input object. Analysis is performed either for one group comparison 
#' in paired design or for any number of unpaired group comparisons. The group designations in 
#' the \code{groupColumn} must match the contrasts given in \code{comparisons} (\code{"groupA-groupB"}).
#' P-value and foldchange thresholds may be applied to filter result data. 
#' Volcano plots are generated for each group comparison given in \code{comparisons} with highligted 
#' significance thresholds.
#' Heatmaps with sample signal intensities are generated for top differentially expressed genes for each 
#' group comparison. Additionally heatmaps indicating foldchanges can be generated 
#' for each set of group comparisons given in \code{FC.heatmap.comparisons}. Selection of probes with respect 
#' to these group comparisons is characterised in \code{FC.heatmap.geneselection}. If the resulting number 
#' of probes exceeds \code{maxHM}, probes are prioritized by either F-Test p-value or by minimum p-value of
#' selected group comparisons.
#' 
#' 
#' @param data ExpressionSet, SummarizedExperiment, DESeqDataSet or MethylSet.
#' @param comparisons character vector with group comparisons in format \code{"groupA-groupB"} or 
#'              nested comparisons in format \code{"(groupA-groupB)-(groupC-groupD)"}.
#' @param p.value.threshold numeric p-value threshold 
#' @param adjust.method adjustment method for multiple testing (\code{"none", "BH", "BY" and "holm"})
#' @param FC.threshold numeric foldchange threshold. If data is log2-transformed, mind to transform threshold, too.
#' 
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' 
#' @param matchvar (character) if paired design, name of column to be used to group samples. \code{NULL} for unpaired design.
#' @param Symbol.column character with column name of Gene Symbols in \code{data} or \code{NULL}.
#' @param sampleColumn character with column name of Sample names in \code{data}
#' @param groupColumn character with column name of group names in \code{data}. Group names must match \code{comparisons}!
#' @param useWeights (boolean) if \code{TRUE}, arrays are weighted for their quality.
#' @param geneAnno2keep (character vector) annotation columns in feature data to keep in output
#' @param venn_comparisons character vector or (named) list of character vectors with group comparisons to be included 
#'                   in Venn Diagram (max 5 comparisons per Venn diagramm). Must be subset of \code{comparisons} or \code{NULL}.
#'
#' @param maxHM (numeric) max number of top diff. regulated elements printed in Heatmap
#' @param scale character indicating if the values should be centered and scaled in either the row direction 
#' or the column direction, or none (default). Either \code{"none"}, \code{"row"} or \code{"column"}.
#' @param HMcexRow labelsize for rownames in Heatmap
#' @param HMcexCol labelsize for colnames in Heatmap
#' @param figure.res numeric resolution for png.
#' @param HMincludeRelevantSamplesOnly (boolean) if \code{TRUE} include only Samples in heatmap, which belong to the groups of the 
#' respective group comparison. If \code{FALSE}, all samples are plotted in the heatmaps. If a custom selection of sample groups 
#' is required for the heatmap of each group comparison, \code{HMincludeRelevantSamplesOnly} can be a named list of the form 
#' \code{list("groupA-groupB" = c("groupA", "groupB", "groupX"))}
#' @param color.palette select color palette for heatmaps                     
#' @param FC.heatmap.comparisons character vector or (named) list of character vectors with group comparsions to be 
#'                         included in foldchange heatmaps. One FC heatmap is generated for each character vector.
#'                         If \code{NULL} no FC heatmaps are generated.
#' @param FC.heatmap.geneselection character vector containing elements of \code{c("all", "intersect", "union")}. 
#'                           Which probes shall be used for foldchange heatmaps? \code{"allprobes"}: no restriction,  
#'                           \code{"intersect"}: only probes are used which are differentially expressed in all respective group comparisons.
#'                           \code{"union"}: only probes are used which are differentially expressed in any of the respective group comparisons.
#'                           Probes are prioritised by F-Test p-value calculated for all group comparisons denoted in \code{comparisons}
#'                           as well as by minimum adjusted p-value of respective group comparisons plotted in this heatmap.
#'                           Obsolete if \code{FC.heatmap.comparisons} is \code{NULL}.
#'                           
#'
#' @return List of 2 elements 
#' \itemize{
#'   \item DEgenes list of dataframes with differential expressed genes according to significance 
#'   thresholds for each comparison (sorted by p-value). List elements are named by the respective group comparison. 
#' \item DEgenes.unfilt list of dataframes with unfiltered differential expression data for each comparison (sorted by p-value).
#' }
#' An F-Test and a Venn Diagramm of all comparisons as well as filtered and unfiltered result tables, heatmaps and volcano plots
#' for each group comparison are stored as side-effects in the project folder.
#' 
#' @author Frank Ruehle
#' 
#' @export 



wrapLimma <- function(data, 
                      comparisons,
                      p.value.threshold = 0.05, 
                      adjust.method="BH", 
                      FC.threshold = log2(1.5),
                      
                      projectfolder = file.path("GEX/limma"),
                      projectname = "", 
                      
                      matchvar=NULL, 
                      Symbol.column = "SYMBOL",
                      sampleColumn = "Sample_Name",   
                      groupColumn= "Sample_Group", 
                      useWeights=TRUE, 
                      geneAnno2keep=NULL, 
                      venn_comparisons= NULL, 

                      # Heatmap parameter:
                      maxHM=50, 
                      scale = c("none"),
                      HMcexRow=1.1, 
                      HMcexCol=1.1, 
                      figure.res = 300,
                      HMincludeRelevantSamplesOnly=TRUE,
                      color.palette="heat.colors",
                      FC.heatmap.comparisons=NULL,
                      FC.heatmap.geneselection = c("allprobes", "intersect", "union")
                      ) {


                        
  
 
  
  
  # load required packages. Packages are not detached afterwards
  pkg.cran <- c("gplots", "VennDiagram")
  pkg.bioc <- c("limma", "Biobase")
  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  projectname <- if (!is.null(projectname)  && projectname!="" && !grepl("_$", projectname)) {paste0(projectname, "_")} else {""}
  
  
  # Create output directories if not yet existing 
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
  if (!file.exists(file.path(projectfolder, "limma_unfiltered"))) {dir.create(file.path(projectfolder, "limma_unfiltered")) }
  if (!file.exists(file.path(projectfolder, "limma_filtered"))) {dir.create(file.path(projectfolder, "limma_filtered")) }
  if (!file.exists(file.path(projectfolder, "Heatmaps"))) {dir.create(file.path(projectfolder, "Heatmaps")) }
  if (!file.exists(file.path(projectfolder, "Volcano_plots"))) {dir.create(file.path(projectfolder, "Volcano_plots")) }
  
  
  
  #### data preparation for expression or methylation data input

  
  
  if(class(data) %in% c("SummarizedExperiment", "DESeqDataSet", "DESeqTransform")) {
    cat("\nClass", class(data), "detected\n")
    Exp <- "GEX"
    phenodata = colData(data) # phenotype data 
    featuredata <- as.data.frame(rowData(data,use.names=TRUE))
    featuredata <- data.frame(rowfeatures= rownames(featuredata), featuredata)
    assaydata <- assay(data)
    heatmap_legend_xaxis = "read counts"
    
  } else {

    if (grepl("expression", class(data), ignore.case=T) ) {
      cat("\nClass", class(data), "detected\n")
      Exp <- "GEX"
      phenodata <- pData(data)
      if("Status" %in% colnames(fData(data))) { # define regular probes if control probes still present in expression set
          ids_regular <- rownames(fData(data)[fData(data)$Status=="regular",])
          } else {ids_regular <- 1:nrow(exprs(data))}
      featuredata <- fData(data)[ids_regular,] # remove control probes from expression set
      assaydata <- exprs(data)[ids_regular,]
      heatmap_legend_xaxis="log2(expression)" # for heatmap
      # pvalue <- "P.Value" # character with column name of p values
        }  else {  
        
        if(grepl("methyl", class(data), ignore.case=T) ) {
          cat("\nClass", class(data), "detected\n")
          Exp <- "MT"
          phenodata <- pData(data)
          featuredata <- mcols(data, use.names=T)
          assaydata <- getBeta(data)
          heatmap_legend_xaxis="beta values" # for heatmap
          # pvalue <- "P.Value" # character with column name of p values
        } else {return("wrong object class!")}
      }
  }
  
  
  # modify plot labels (heatmap legend)
  heatmap_legend_xaxis <- "log2(expression)"
  if (grepl("row", scale)) {heatmap_legend_xaxis <- "row z-score"}
  if (grepl("col", scale)) {heatmap_legend_xaxis <- "column z-score"}
  

  # if no thresholds given for p-value and foldchange, values are applied to omit filtering
  if(is.null(p.value.threshold)) {p.value.threshold <- 1}
  if(is.null(FC.threshold)) {FC.threshold <- 0}
  
  

# Remark paired design from limma vignette:
# The approach used for paired samples can be applied in any situation where there are batch
# effects or where the experiment has been conducted in blocks. The treatments can be adjusted for
# differences between the blocks by using a model formula of the form:
#  design <- model.matrix(~Block+Treatment)
# In this type of analysis, the treatments are compared only within each block.


  ### Prepare Design Matrix:

  # no contrasts generated from a 0/1 phenotype (levels must by syntactically valid names in R). Therefore renamed.
  if(is.numeric(phenodata[,groupColumn]) & all(phenodata[,groupColumn] %in% c(0,1))) {
      rna <- factor(make.names(phenodata[,groupColumn])) # c("0","1") -> c("X0","X1")
           } else {
        rna <- factor(phenodata[,groupColumn])
        }
  
  
  if (!is.null(matchvar)) {
        cat("\nPrepare design matrix for paired design\n")
        matched <- phenodata[,matchvar] # wenn paired design
        design <- model.matrix(~matched+rna)
  } else {  
    cat("\nPrepare design matrix\n")
    design <- model.matrix(~0+rna) # default: no paired design
    colnames(design) <- levels(rna)
        }
  
  # array weights: arrays weighted for their quality
  if(useWeights) {
    cat("\nArray quality weighting\n")
    aw <- arrayWeights(assaydata, design)  
      } else { aw <- NULL}





### Prepare fit: "MArrayLM"-Object ready for diff analysis.
# input for lmfit: A matrix-like data object containing log-ratios or log-expression values for a series of arrays
  if (!is.null(matchvar)) {fit <- lmFit(assaydata, design, weights=aw) # if paired design  
                         fit <- eBayes(fit)
    } else { # if no paired design
      contrasts <- makeContrasts(contrasts=comparisons, levels=design) 
      fit <- lmFit(assaydata, design, weights=aw)  
      fit <- eBayes(contrasts.fit(fit, contrasts))
      }


# Venn-Diagramm for group comparisons:
  if(!is.null(venn_comparisons)) {
    if(!is.list(venn_comparisons)) {venn_comparisons <- list(venn_comparisons)}
    for(v in 1:length(venn_comparisons)) {
      
      if(is.null(names(venn_comparisons)[v])) {names(venn_comparisons)[v] <- paste0("Vennset",v)}
      filename.Venn <- file.path(projectfolder, paste0("Venn_Diagram_", projectname, names(venn_comparisons)[v], ".png"))
      cat("\nWrite Venn diagramm to", filename.Venn, "\n")
  
     resultDiffGenes <- decideTests(fit, p.value=p.value.threshold, lfc=FC.threshold, adjust.method=adjust.method) # default adjust.method="BH"
    if (length(venn_comparisons[[v]]) <= 5) {
      png(filename.Venn, width = 300, height = 300, units = "mm", res=figure.res)
      # pdf(filename.Venn, width = 7, height = 7) 
      vennDiagram(resultDiffGenes[,colnames(resultDiffGenes) %in% venn_comparisons[[v]]], names=NULL)  #optional: select desired group comparisons for Venn diagramm (max: 5)!
      dev.off()
      }
  } # end v-loop
}
  
# Diff expressed genes caculated by F-Test via all groups (contr.fit$F.p.value)
# topTable with coef=NULL is the same as topTableF, unless the fitted model fit has only one column.
# DiffAllGroups_Ftest contains entries with at least one absolute log-fold-changes greater than lfc. 
# Filtering for F-Test p-values is done when writing table. 
# Columns with names of group comparisons represent log foldchanges
  DiffAllGroups_Ftest <- topTable(fit, coef=NULL, confint=TRUE, p.value=1, lfc=FC.threshold, 
                                  adjust.method=adjust.method, number=Inf, 
                                  genelist=featuredata[,names(featuredata) %in% geneAnno2keep, drop=F])
  if(Exp=="MT") {names(DiffAllGroups_Ftest)[names(DiffAllGroups_Ftest)=="AveExpr"] <- "AveBeta"}
  
  DiffAllGroups_Ftest <- data.frame(PROBE_ID=rownames(DiffAllGroups_Ftest), DiffAllGroups_Ftest)
  # in columns names of comparisons hyphend are replaced by dots.
  # This is reversed here and "logFC_" added to color names of comparisons
  ftest.columns <- sub("\\.", "-", names(DiffAllGroups_Ftest))
  columns2replace <- ftest.columns %in% comparisons
  ftest.columns <- paste0("logFC_", ftest.columns)
  names(DiffAllGroups_Ftest)[columns2replace] <- ftest.columns[columns2replace]
  
  
  cat("\nWrite F-Test result to", file.path(projectfolder, "limma_filtered", paste0(projectname, "all_Groups_Ftest.txt")), "\n")
  DiffAllGroups_Ftest.sign <- DiffAllGroups_Ftest[DiffAllGroups_Ftest$adj.P.Val<p.value.threshold,]
  write.table(DiffAllGroups_Ftest.sign, 
              file= file.path(projectfolder, "limma_filtered", paste0(projectname, "all_Groups_Ftest.txt")), sep="\t", quote=F, row.names=F)
  
  
  
      ### Heatmap with all samples. Genes selected by F-Test
      if(!is.null(DiffAllGroups_Ftest.sign)) {
        if(nrow(DiffAllGroups_Ftest.sign) >=2) {
          cat("\nWrite F-Test Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "all_Samples_Ftest.png", sep="" )), "\n")  
        
        if (nrow(DiffAllGroups_Ftest.sign) > maxHM) {DEgenesHM.Ftest <- DiffAllGroups_Ftest.sign[1:maxHM,]} else {DEgenesHM.Ftest <- DiffAllGroups_Ftest.sign}
        plotmatrix <- assaydata[rownames(assaydata) %in% rownames(DEgenesHM.Ftest), ,drop=F]
        
        # If Symbols are available, rows are annotated with symbols instead of probe IDs   
        indexRownames <- match(rownames(plotmatrix), rownames(featuredata))
        
        if(!is.null(Symbol.column)) {
          rownames(plotmatrix) <- ifelse(featuredata[indexRownames,Symbol.column] != "" & !is.na(featuredata[indexRownames,Symbol.column]), 
                                         as.character(featuredata[indexRownames,Symbol.column]), rownames(plotmatrix))
        }
        
        groupColorCode <- rainbow(length(unique(phenodata[,groupColumn])))[as.numeric(factor(phenodata[,groupColumn]))] # for ColSideColors in heatmaps
        
        png(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "all_Samples_Ftest.png", sep="" )), width = 210, height = 297, units = "mm", res=figure.res) 
        #pdf(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "all_Samples_Ftest.pdf", sep="" )), width = 10, height = 14) 
        heatmap.2(plotmatrix, main=paste0(projectname, " all_Samples_Ftest"),  
                  margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both", scale = scale, 
                  ColSideColors= groupColorCode,
                  cexRow=HMcexRow, cexCol=HMcexCol, trace="none", density.info="histogram",  
                  key.xlab= heatmap_legend_xaxis, key.ylab=NA, key.title="Color Key", col=color.palette)  
        dev.off()
      }
  }

## Diff expressed genes caculated by ANOVA for all groups 
      
      aov.result <- apply(assaydata, 1, function(x) {
        unlist(summary(aov(x ~ phenodata[,groupColumn]))[[1]][1, "Pr(>F)"])})  # ANOVA for each gene results in p-value vector
      
      aov.result.sign <- aov.result[aov.result< p.value.threshold]  # filtering and sorting ANOVA result
      aov.result.sign <- sort(aov.result.sign)
      cat("\nWrite ANOVA result to", file.path(projectfolder, "limma_filtered", paste0(projectname, "all_Groups_ANOVA.txt")), "\n")
      # write.table(as.data.frame(aov.result.sign), row.names=T, col.names = "Pr(>F)",
      #             file= file.path(projectfolder, "limma_filtered", paste0(projectname, "all_Groups_ANOVA.txt")), sep="\t", quote=F)
      write.table(data.frame(feature=names(aov.result.sign), p=aov.result.sign), row.names=F, 
                  file= file.path(projectfolder, "limma_filtered", paste0(projectname, "all_Groups_ANOVA.txt")), sep="\t", quote=F)
      
      
      
      ### Heatmap with all samples. Genes selected by ANOVA
      if(length(aov.result.sign) >=2) {
        cat("\nWrite ANOVA Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "all_Samples_ANOVA.png", sep="" )), "\n")  
     
         maxgenesANOVA <- maxHM # 1000
      
         if (length(aov.result.sign) > maxgenesANOVA) {aov.result.sign <- aov.result.sign[1:maxgenesANOVA]} 
        plotmatrix.aov <- assaydata[rownames(assaydata) %in% names(aov.result.sign), ,drop=F]

      # If Symbols are available, rows are annotated with symbols instead of probe IDs   
      indexRownames <- match(rownames(plotmatrix.aov), rownames(featuredata))
      
      if(!is.null(Symbol.column)) {
        rownames(plotmatrix.aov) <- ifelse(featuredata[indexRownames,Symbol.column] != "" & !is.na(featuredata[indexRownames,Symbol.column]), 
                                       as.character(featuredata[indexRownames,Symbol.column]), rownames(plotmatrix.aov))
      }
      
      groupColorCode <- rainbow(length(unique(phenodata[,groupColumn])))[as.numeric(factor(phenodata[,groupColumn]))] # for ColSideColors in heatmaps
      
      png(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "all_Samples_ANOVA.png", sep="" )), width = 210, height = 297, units = "mm", res=figure.res) 
      #pdf(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "all_Samples_ANOVA.pdf", sep="" )), width = 10, height = 14) 
      heatmap.2(plotmatrix.aov, main=paste0(projectname, " all_Samples_ANOVA"),  
                margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both", # labRow=F, 
                ColSideColors= groupColorCode,
                cexRow=HMcexRow, cexCol=HMcexCol, scale = scale, trace="none", density.info="histogram",  
                key.xlab= heatmap_legend_xaxis, key.ylab="", key.title="Color Key", col=color.palette)  
      dev.off()  
      
      }
      
   
      
         

######### differential group comparisons
DEgenes.unfilt <- list()
DEgenes <- list()

# loop for all dedicated group comparisons
for (i in 1:length(comparisons)) {
  cat("\n\nProcessing comparison", i, ":", comparisons[i])
  
  # Calculate differential expressed genes for each group comparison (filtered and unfiltered)
  # coef is column of investigation from design matrix. In matched design the first two columns are "(Intercept)" and "matched".
  if (is.null(matchvar)) {coef = comparisons[i]} else {coef = 2+i}  
  
  DEgenes.unfilt[[comparisons[i]]] <- topTable(fit, coef=coef, confint=TRUE, number=Inf, sort.by="p",
                                               adjust.method=adjust.method, 
                                               genelist=featuredata[,names(featuredata) %in% geneAnno2keep, drop=F])
  
  if(Exp=="MT") {names(DEgenes.unfilt[[i]])[names(DEgenes.unfilt[[i]])=="AveExpr"] <- "AveBeta"}
  
  
  ### Filtering gene list
  DEgenes[[comparisons[i]]] <- filterGeneLists(DEgenes.unfilt[[comparisons[i]]],
                                               newheader=NULL,
                                               filtercat1 = "adj.P.Val",
                                               filtercat1.decreasing = FALSE,
                                               filtercat1.function = identity,
                                               filtercat1.threshold=p.value.threshold,
                                               filtercat2 = "logFC",
                                               filtercat2.decreasing = TRUE,
                                               filtercat2.function = abs,
                                               filtercat2.threshold = FC.threshold)
  
  

  # Writing result gene tables
  helptable.unfilt <- data.frame(PROBE_ID=rownames(DEgenes.unfilt[[comparisons[i]]]), DEgenes.unfilt[[comparisons[i]]])
  helptable.filt   <- data.frame(PROBE_ID=rownames(DEgenes[[comparisons[i]]]), DEgenes[[comparisons[i]]])
  
  write.table(helptable.unfilt, sep="\t", quote=F, row.names=F,
              file= file.path(projectfolder, paste0("limma_unfiltered"), paste0(projectname, comparisons[i], "_unfilt.txt")))
  
  write.table(helptable.filt, sep="\t", quote=F, row.names=F,
              file= file.path(projectfolder, paste0("limma_filtered"), paste0(projectname, comparisons[i], ".txt")))
  
  cat(paste(nrow(helptable.filt),"differentially regulated elements for comparison:",comparisons[i]))
  cat(paste("\nWrite gene tables to", file.path(projectfolder, paste0("limma_unfiltered"), paste0(projectname, comparisons[i], "_unfilt.txt")),
            "and", file.path(projectfolder, paste0("limma_filtered"), paste0(projectname, comparisons[i], ".txt"))))
  
  
 
  ######## Volcano plots per group comparison 
  
  volcanoplot_filename <- file.path(projectfolder, "Volcano_plots", paste("Volcano_", projectname, comparisons[i], ".png", sep="" ))
  
  cat("\nWrite Volcano plot to", volcanoplot_filename)  
  
  png(file=volcanoplot_filename, width = 150, height = 150, units = "mm", res=figure.res) 
  plot(DEgenes.unfilt[[comparisons[i]]]$logFC, -log10(DEgenes.unfilt[[comparisons[i]]]$adj.P.Val), main=comparisons[i],
       col="gray", pch=16, cex=0.5, xlab="fold change", ylab="-log10 p-value") 
  if(!is.null(p.value.threshold)) {abline(h=-log10(p.value.threshold),lty=c(1))}
  if(!is.null(FC.threshold)) {abline(v=c(FC.threshold, -FC.threshold),lty=2)}
  
  if(nrow(DEgenes[[comparisons[i]]]) >=1) {  # no highlighting if no diff expressed genes
    points(DEgenes[[comparisons[i]]]$logFC, -log10(DEgenes[[comparisons[i]]]$adj.P.Val),pch=16, col="black", cex=0.8) 
  } 
  dev.off()
  
  
  
  
   # add adjusted p-value to 'DiffAllGroups_Ftest' for later foldchange heatmaps
  DiffAllGroups_Ftest[,paste0("p_", comparisons[i])] <- helptable.unfilt[match(rownames(DiffAllGroups_Ftest), rownames(helptable.unfilt)),"adj.P.Val"]
  
  
  
  ######## Heatmaps per group comparison with signal intensities
      if(nrow(DEgenes[[i]]) >1) {  # no Heatmap if just one diff expressed gene
        
        cat("\nWrite Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_", comparisons[i], ".png", sep="" )), "\n")  
        
        if (nrow(DEgenes[[i]]) > maxHM) {DEgenesHM <- DEgenes[[i]][1:maxHM,]} else {DEgenesHM <- DEgenes[[i]]}
        
        plotmatrix <- assaydata[rownames(assaydata) %in% rownames(DEgenesHM), ,drop=F]
        
        # If Symbols are available, rows are annotated with symbols instead of probe IDs   
        indexRownames <- match(rownames(plotmatrix), rownames(featuredata))

        if(!is.null(Symbol.column)) {
          rownames(plotmatrix) <- ifelse(featuredata[indexRownames,Symbol.column] != "" & !is.na(featuredata[indexRownames,Symbol.column]), 
                                         as.character(featuredata[indexRownames,Symbol.column]), rownames(plotmatrix))
          }
        
        groupColorCode <- rainbow(length(unique(phenodata[,groupColumn])))[as.numeric(factor(phenodata[,groupColumn]))] # for ColSideColors in heatmaps
        

          #  if HMincludeRelevantSamplesOnly=F, all samples are plotted in every group comparison
          if (isTRUE(HMincludeRelevantSamplesOnly) | is.list(HMincludeRelevantSamplesOnly)) {
            
            if (is.list(HMincludeRelevantSamplesOnly)) { 
              if(is.null(HMincludeRelevantSamplesOnly[[comparisons[i]]])) {stop(paste("\nname", comparisons[i], "not found as list element in HMincludeRelevantSamplesOnly"))} 
              groups2plot <- HMincludeRelevantSamplesOnly[[comparisons[i]]] # plot samples of groups given as character in HMincludeRelevantSamplesOnly
            } else {
              groups2plot <- unlist(strsplit(comparisons[i], "-") )  # plot only samples considered belonging to the respective group comparison 
              groups2plot <- unique(sub("[(|)]", "", groups2plot))  # for comparison of comparison
            }
           
          sampleTable <- phenodata[,c(sampleColumn,groupColumn)]
          samples2plot <- character()
          
          for (j in 1:length(groups2plot)) {
            samples2plot <- c(samples2plot, as.character(sampleTable[sampleTable[,groupColumn]==groups2plot[j], sampleColumn])) }
          
          cols2plot <- colnames(plotmatrix) %in% samples2plot
          plotmatrix <- plotmatrix[, cols2plot]
       
          groupColorCode <- groupColorCode[phenodata[,groupColumn] %in% groups2plot] # adjust ColSideColors in heatmaps to relevant groups
         }
        
        png(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, comparisons[i], ".png", sep="" )), width = 210, height = 297, units = "mm", res=figure.res) 
        #pdf(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, comparisons[i], ".pdf", sep="" )), width = 10, height = 14) 
        heatmap.2(plotmatrix, main=comparisons[i], margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both",  
                  ColSideColors= groupColorCode, 
                  cexRow=HMcexRow, cexCol=HMcexCol, scale = scale, trace="none", density.info="histogram",  
                  key.xlab=heatmap_legend_xaxis, key.ylab="", key.title="Color Key", col=color.palette)  
        
        # Settings for plotting color key instead of dendrogram on top of heatmap:  
        # lmat=rbind(c(0,3),c(0,4), c(2,1)), lhei = c(0.3,0.5,3.8), lwid = c(0.5,4), key.par=list(mar=c(4,2,2,13)), Rowv=TRUE, Colv=FALSE, dendrogram="row", keysize=0.8 
        dev.off()
      } # end if nrow(DEgenes[[i]]) >1
   

  
} # end i-loop for group comparisons
     

  
  ###### Heatmaps with selected group Foldchanges 
  # (use 'DiffAllGroups_Ftest' and 'resultDiffGenes' from above)
  # DiffAllGroups_Ftest contains all probes with at least one absolute log-fold-changes greater than lfc. 
  if(!is.null(FC.heatmap.comparisons)) {

    # modify plot labels
    heatmap_legend_xaxis <- "log2(FC)"
    if (grepl("row", scale)) {heatmap_legend_xaxis <- "row z-score (FC)"}
    if (grepl("col", scale)) {heatmap_legend_xaxis <- "column z-score (FC)"}
    
    
    if(!is.list(FC.heatmap.comparisons)) {FC.heatmap.comparisons <- list(FC.heatmap.comparisons)}
    
    for(fc in 1:length(FC.heatmap.comparisons)) {

      for(g in FC.heatmap.geneselection) { 
          # define probes which are intersections or unions for respective group comparisons
          probes2use <- switch(g, allprobes=rownames(resultDiffGenes),
                         intersect=rownames(resultDiffGenes[apply(resultDiffGenes[,colnames(resultDiffGenes) %in% FC.heatmap.comparisons[[fc]]], 1, function(x) {all(x!=0)}),,drop=F]),
                         union=rownames(resultDiffGenes[apply(resultDiffGenes[,colnames(resultDiffGenes) %in% FC.heatmap.comparisons[[fc]]], 1, function(x) {any(x!=0)}),,drop=F])
                         )
         
          # top genes are filtered for F-Test p-value or minimum p-value from selected group comparison
          for(f in c("Ftest", "minPvalue")) {
            column2sort <- switch(f, Ftest="adj.P.Val", minPvalue="minimumP")
            DiffAllGroups_Ftest$minimumP <- apply(DiffAllGroups_Ftest[,names(DiffAllGroups_Ftest) %in% paste0("p_", FC.heatmap.comparisons[[fc]]),drop=F], 1, min)
            DiffAllGroups_Ftest <- DiffAllGroups_Ftest[order(DiffAllGroups_Ftest[,column2sort], decreasing=F),]
            names(DiffAllGroups_Ftest) <- sub("logFC_", "", names(DiffAllGroups_Ftest))
                
              # make plotmatrix with relevant probes and relevant group comparisons
              plotmatrix <- as.matrix(DiffAllGroups_Ftest[rownames(DiffAllGroups_Ftest) %in% probes2use, 
                                                          names(DiffAllGroups_Ftest) %in% FC.heatmap.comparisons[[fc]], drop=F])
            
            # If Symbols are available, rows are annotated with symbols instead of probe IDs
            indexRownames <- match(rownames(plotmatrix), rownames(featuredata))
            if(!is.null(Symbol.column)) {
              rownames(plotmatrix) <- ifelse(featuredata[indexRownames,Symbol.column] != "" & !is.na(featuredata[indexRownames,Symbol.column]),
                                             as.character(featuredata[indexRownames,Symbol.column]), rownames(plotmatrix))
            }

            if(nrow(plotmatrix) >1) {  # no Heatmap if just one diff expressed gene
              if (nrow(plotmatrix) > maxHM) {plotmatrix <- plotmatrix[1:maxHM,]} 
              if(is.null(names(FC.heatmap.comparisons)[fc])) {names(FC.heatmap.comparisons)[fc] <- paste0("set",fc)}
              filename.FCHM <- file.path(projectfolder, "Heatmaps", paste0("Heatmap_Foldchanges_", projectname, names(FC.heatmap.comparisons)[fc], "_", g, "_prior_by_", f, ".png"))
              cat("\nWrite Heatmap to", filename.FCHM, "\n")  
              
              png(filename.FCHM, width = 210, height = 297, units = "mm", res=figure.res) 
              #pdf(filename.FCHM, width = 10, height = 14) 
              heatmap.2(plotmatrix, main=paste0("Heatmap Foldchanges ", projectname, names(FC.heatmap.comparisons)[fc], " ", g, "\nprobes prioritised by ", f), 
                        margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both",  
                        cexRow=HMcexRow, cexCol=HMcexCol, scale = scale, trace="none", density.info="histogram",  
                        key.xlab= heatmap_legend_xaxis, key.ylab="", key.title="Color Key", col=color.palette)  
              # Settings for plotting color key instead of dendrogram on top of heatmap:  
              # lmat=rbind(c(0,3),c(0,4), c(2,1)), lhei = c(0.3,0.5,3.8), lwid = c(0.5,4), key.par=list(mar=c(4,2,2,13)), Rowv=TRUE, Colv=FALSE, dendrogram="row", keysize=0.8 
              dev.off()
            
        } # end if nrow(plotmatrix) >1
      } # end f-loop
    } # end g-loop
  } # end fc-loop 
} # end if(!is.null(FC.heatmap.comparisons))
 ############ end foldchange heatmaps 
  
 

# return unfiltered genes as list of dataframes for each comparison (sorted by p-value)
return(list(DEgenes=DEgenes, DEgenes.unfilt=DEgenes.unfilt))


}


