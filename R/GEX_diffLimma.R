
## Description
# Differential expression analysis with limma

## Usage 
diffLimma <- function(GEXMTSet, 
                      comparisons,
                      p.value.threshold = 0.05, 
                      adjust.method="BH", 
                      FC.threshold = log2(1.5),
                      
                      projectfolder = file.path("GEX/Diff_limma"),
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
                      HMcexRow=1.1, 
                      HMcexCol=1.1, 
                      HMincludeRelevantSamplesOnly=TRUE,
                      color.palette="heat.colors",
                      FC.heatmap.comparisons=NULL,
                      FC.heatmap.geneselection = c("allprobes", "intersect", "union")
                      ) {


                        
  # ## Arguments
  # GEXMTSet: ExpressionSet or MethylSet
  # comparisons: character vector with group comparisons in format "groupA-groupB" or 
  #              nested comparisons in format "(groupA-groupB)-(groupC-groupD)".
  # p.value.threshold: numeric p-value threshold 
  # adjust.method: adjustment method for multiple testing ("none", "BH", "BY" and "holm")
  # FC.threshold: numeric foldchange threshold. If data is log2-transformed, mind to transform threshold too.
  # 
  # projectfolder: character with directory for output files (will be generated if not exisiting).
  # projectname: optional character prefix for output file names.
  # 
  # matchvar: (character) if paired design, name of column to be used to group samples. NULL for unpaired design.
  # Symbol.column: character with column name of Gene Symbols in 'DEgenes.unfilt' and 'GEXMTSet' or NULL.
  # sampleColumn: character with column name of Sample names in pData(GEXMTSet)
  # groupColumn: character with column name of group names in pData(GEXMTSet). Group names must match comparisons!
  # useWeights: (boolean) if TRUE, arrays are weighted for their quality.
  # geneAnno2keep: (character vector) annotation columns in feature data to keep in output
  # venn_comparisons: character vector or (named) list of character vectors with group comparisons to be included 
  #                   in Venn Diagram (max 5 comparisons per Venn diagramm)
  #
  # # Heatmap parameter:
  # maxHM: (numeric) max number of top diff. regulated elements printed in Heatmap
  # HMcexRow: labelsize for rownames in Heatmap
  # HMcexCol: labelsize for colnames in Heatmap
  # HMincludeRelevantSamplesOnly: (boolean) if TRUE include only Samples in heatmap, 
  #                                which are in groups of the respective group comparison.  
  # color.palette: select color palette for heatmaps                     
  # FC.heatmap.comparisons: character vector or (named) list of character vectors with group comparsions to be 
  #                         included in foldchange heatmaps. One FC heatmap is generated for each character vector.
  #                         If NULL no FC heatmaps are generated.
  # FC.heatmap.geneselection: character vector containing elements of c("all", "intersect", "union"). 
  #                           Which probes shall be used for foldchange heatmaps? "allprobes": no restriction,  
  #                           "intersect": only probes are used which are differentially expressed in all respective group comparisons.
  #                           "union": only probes are used which are differentially expressed in any of the respective group comparisons.
  #                           Probes are prioritised by F-Test p-value calculated for all group comparisons 
  #                           as well as by minimum adjusted p-value of respective group comparisons
  #                           Obsolete if 'FC.heatmap.comparisons' is NULL.
  
  
  ## Details
  # Function generates design matrix, fit and differential expression results for dedicated group comparisons of input object.
  # Analysis is performed either for one group comparison in paired design or for any number of unpaired group comparisons.
  # The group designations in the groupColumn must match the contrasts given in comparisons ("groupA-groupB").
  # P-value and foldchange thresholds may be applied to filter result data. 
  # Heatmaps with sample signal intensities are generated for top differentially expressed genes for each group comparison 
  # given in 'comparisons'. Additionally heatmaps indicating foldchanges are generated for each set of group comparisons
  # given in 'FC.heatmap.comparisons'. Selection of probes with respect to these group comparisons is characterised in 
  # 'FC.heatmap.geneselection'. If the resulting number of probes exceeds 'maxHM', probes are prioritized by either 
  # F-Test p-value or by minimum p-value of selected group comparisons.
  
  
  ## Value 
  # List of dataframes with unfiltered differential expression data for each comparison (sorted by p-value).
  # List elements are named by the respective group comparison. 
  # An F-Test and a Venn Diagramm of all comparisons as well as filtered and unfiltered result tables and heatmaps 
  # for each group comparison are stored as side-effects in the project folder.
  
  
  ## Author(s) 
  # Frank Rühle 
  
  
  
  # load required packages. Packages are not detached afterwards
  pkg.cran <- c("gplots")
  pkg.bioc <- c("limma", "Biobase")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  # Create output directories if not yet existing 
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
  if (!file.exists(file.path(projectfolder, "Diff_limma_unfiltered"))) {dir.create(file.path(projectfolder, "Diff_limma_unfiltered")) }
  if (!file.exists(file.path(projectfolder, "Heatmaps"))) {dir.create(file.path(projectfolder, "Heatmaps")) }
  
  
  
  #### data preparation for expression or methylation data input
  phenoGEXMTSet <- pData(GEXMTSet)
  
  if (grepl("expression", class(GEXMTSet), ignore.case=T) ) {
    Exp <- "GEX"
    if("Status" %in% colnames(fData(GEXMTSet))) { # define regular probes if control probes still present in expression set
        ids_regular <- rownames(fData(GEXMTSet)[fData(GEXMTSet)$Status=="regular",])
        } else {ids_regular <- 1:nrow(exprs(GEXMTSet))}
    featureGEXMTSet <- fData(GEXMTSet)[ids_regular,] # remove control probes from expression set
    dataGEXMTSet <- exprs(GEXMTSet)[ids_regular,]
    plot.label="log2(expression)" # for heatmap
    # pvalue <- "P.Value" # character with column name of p values
      }  else {  
      
      if(grepl("methyl", class(GEXMTSet), ignore.case=T) ) {
        Exp <- "MT"
        featureGEXMTSet <- mcols(GEXMTSet, use.names=T)
        dataGEXMTSet <- getBeta(GEXMTSet)
        plot.label="beta values" # for heatmap
        # pvalue <- "P.Value" # character with column name of p values
      } else {return("wrong object class!")}
    }
  
   
  
  
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

  # no contrasts generated from a 0/1 phenotype. Therefor recoded to "y"/"n"
  if(is.numeric(phenoGEXMTSet[,groupColumn]) & all(phenoGEXMTSet[,groupColumn] %in% c(0,1))) {
      rna <- factor(ifelse(phenoGEXMTSet[,groupColumn]==0, "n", "y"))
           } else {
        rna <- factor(phenoGEXMTSet[,groupColumn])
        }
  
  
  if (!is.null(matchvar)) {
        cat("\nPrepare design matrix for paired design\n")
        matched <- phenoGEXMTSet[,matchvar] # wenn paired design
        design <- model.matrix(~matched+rna)
  } else {  
    cat("\nPrepare design matrix\n")
    design <- model.matrix(~0+rna) # Normalfall ohne paired design
    colnames(design) <- levels(rna)
        }
  
  # array weights: arrays weighted for their quality
  if(useWeights) {
    cat("\nArray quality weighting\n")
    aw <- arrayWeights(dataGEXMTSet, design)  
      } else { aw <- NULL}





### Prepare fit: "MArrayLM"-Object ready for diff analysis.
if (!is.null(matchvar)) {fit <- lmFit(dataGEXMTSet, design, weights=aw) # if paired design  
                         fit <- eBayes(fit)
    } else { # if no paired design
      contrasts <- makeContrasts(contrasts=comparisons, levels=design) 
      fit <- lmFit(dataGEXMTSet, design, weights=aw)  
      fit <- eBayes(contrasts.fit(fit, contrasts))
      }


# Venn-Diagramm for group comparisons:
  if(!is.null(venn_comparisons)) {
    if(!is.list(venn_comparisons)) {venn_comparisons <- list(venn_comparisons)}
    for(v in 1:length(venn_comparisons)) {
      
      if(is.null(names(venn_comparisons)[v])) {names(venn_comparisons)[v] <- paste0("Vennset",v)}
      filename.Venn <- file.path(projectfolder, paste0("Venn_Diagram_", projectname, names(venn_comparisons)[v], ".pdf"))
      cat("\nWrite Venn diagramm to", filename.Venn, "\n")
  
     resultDiffGenes <- decideTests(fit, p.value=p.value.threshold, lfc=FC.threshold, adjust.method=adjust.method) # default adjust.method="BH"
    if (length(venn_comparisons[[v]]) <= 5) {
      pdf(filename.Venn, width = 14, height = 7) 
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
                                  genelist=featureGEXMTSet[,names(featureGEXMTSet) %in% geneAnno2keep, drop=F])
  if(Exp=="MT") {names(DiffAllGroups_Ftest)[names(DiffAllGroups_Ftest)=="AveExpr"] <- "AveBeta"}
  
  DiffAllGroups_Ftest <- data.frame(PROBE_ID=rownames(DiffAllGroups_Ftest), DiffAllGroups_Ftest)
  # in columns names of comparisons hyphend are replaced by dots.
  # This is reversed here and "logFC_" added to color names of comparisons
  ftest.columns <- sub("\\.", "-", names(DiffAllGroups_Ftest))
  columns2replace <- ftest.columns %in% comparisons
  ftest.columns <- paste0("logFC_", ftest.columns)
  names(DiffAllGroups_Ftest)[columns2replace] <- ftest.columns[columns2replace]
  
  
  cat("\nWrite F-Test result to", file.path(projectfolder, paste0(projectname, "all_Groups_Ftest.txt")), "\n")
  DiffAllGroups_Ftest.sign <- DiffAllGroups_Ftest[DiffAllGroups_Ftest$adj.P.Val<p.value.threshold,]
  write.table(DiffAllGroups_Ftest.sign, 
              file= file.path(projectfolder, paste0(projectname, "all_Groups_Ftest.txt")), sep="\t", quote=F, row.names=F)
  
  
  
      ### Heatmap with all samples. Genes selected by F-Test
      if(!is.null(DiffAllGroups_Ftest.sign)) {
      cat("\nWrite Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_all_Samples_Ftest.pdf", sep="" )), "\n")  
      
      if (nrow(DiffAllGroups_Ftest.sign) > maxHM) {DEgenesHM.Ftest <- DiffAllGroups_Ftest.sign[1:maxHM,]} else {DEgenesHM.Ftest <- DiffAllGroups_Ftest.sign}
      plotmatrix <- dataGEXMTSet[rownames(dataGEXMTSet) %in% rownames(DEgenesHM.Ftest), ,drop=F]
      
      # If Symbols are available, rows are annotated with symbols instead of probe IDs   
      indexRownames <- match(rownames(plotmatrix), rownames(featureGEXMTSet))
      
      if(!is.null(Symbol.column)) {
        rownames(plotmatrix) <- ifelse(featureGEXMTSet[indexRownames,Symbol.column] != "" & !is.na(featureGEXMTSet[indexRownames,Symbol.column]), 
                                       as.character(featureGEXMTSet[indexRownames,Symbol.column]), rownames(plotmatrix))
      }
      
      pdf(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_all_Samples_Ftest.pdf", sep="" )), width = 10, height = 14) 
      heatmap.2(plotmatrix, main=paste0(projectname, " all_Samples_Ftest"),  
                margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both",  
                cexRow=HMcexRow, cexCol=HMcexCol, trace="none", density.info="histogram",  
                key.xlab="log2(expression)", key.ylab="", key.title="Color Key", col=color.palette)  
      dev.off()
      }

## Diff expressed genes caculated by ANOVA for all groups 
      
      aov.result <- apply(dataGEXMTSet, 1, function(x) {
        unlist(summary(aov(x ~ phenoGEXMTSet[,groupColumn]))[[1]][1, "Pr(>F)"])})  # ANOVA for each gene results in p-value vector
      
      aov.result.sign <- aov.result[aov.result< p.value.threshold]  # filtering and sorting ANOVA result
      aov.result.sign <- sort(aov.result.sign)
      
      cat("\nWrite ANOVA result to", file.path(projectfolder, paste0(projectname, "all_Groups_ANOVA.txt")), "\n")
      write.table(as.data.frame(aov.result.sign), row.names=T, col.names = "Pr(>F)",
                  file= file.path(projectfolder, paste0(projectname, "all_Groups_ANOVA.txt")), sep="\t", quote=F)
      
      
      
      ### Heatmap with all samples. Genes selected by ANOVA
      cat("\nWrite Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_all_Samples_ANOVA.pdf", sep="" )), "\n")  
      maxgenesANOVA <- 1000
      if (length(aov.result.sign) > maxgenesANOVA) {aov.result.sign <- aov.result.sign[1:maxgenesANOVA]} 
      plotmatrix.aov <- dataGEXMTSet[rownames(dataGEXMTSet) %in% names(aov.result.sign), ,drop=F]

      pdf(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_all_Samples_ANOVA.pdf", sep="" )), width = 10, height = 14) 
      heatmap.2(plotmatrix.aov, main=paste0(projectname, " all_Samples_ANOVA"),  
                margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both", labRow=F, 
                cexCol=HMcexCol, trace="none", density.info="histogram",  
                key.xlab="log2(expression)", key.ylab="", key.title="Color Key", col=color.palette)  
      dev.off()  
      
      
      
   
      
         

######### differential group comparisons
DEgenes.unfilt <- list()
DEgenes <- list()

# loop for all dedicated group comparisons
for (i in 1:length(comparisons)) {
  
  # Calculate differential expressed genes for each group comparison (filtered and unfiltered)
  # coef is column of investigation from design matrix. In matched design the first two columns are "(Intercept)" and "matched".
  if (is.null(matchvar)) {coef = comparisons[i]} else {coef = 2+i}  
  
  DEgenes.unfilt[[comparisons[i]]] <- topTable(fit, coef=coef, confint=TRUE, number=Inf, sort.by="p",
                                               adjust.method=adjust.method, 
                                               genelist=featureGEXMTSet[,names(featureGEXMTSet) %in% geneAnno2keep, drop=F])
  
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
              file= file.path(projectfolder, paste0("Diff_limma_unfiltered"), paste0(projectname, comparisons[i], "_unfilt.txt")))
  
  write.table(helptable.filt, sep="\t", quote=F, row.names=F,
              file= file.path(projectfolder, paste0(projectname, comparisons[i], ".txt")))
  
  cat(paste("\n", nrow(helptable.filt),"differentially regulated elements for comparison:",comparisons[i]))
  cat(paste("\nWrite gene tables to", file.path(projectfolder, paste0("Diff_limma_unfiltered"), paste0(projectname, comparisons[i], "_unfilt.txt")),
            "and", file.path(projectfolder, paste0(projectname, comparisons[i], ".txt"))))
  
  
  # add adjusted p-value to 'DiffAllGroups_Ftest' for later foldchange heatmaps
  DiffAllGroups_Ftest[,paste0("p_", comparisons[i])] <- helptable.unfilt[match(rownames(DiffAllGroups_Ftest), rownames(helptable.unfilt)),"adj.P.Val"]
  
  
  
  ######## Heatmaps per group comparison with signal intensities
      if(nrow(DEgenes[[i]]) >1) {  # no Heatmap if just one diff expressed gene
        
        cat("\nWrite Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_", comparisons[i], ".pdf", sep="" )), "\n")  
        
        if (nrow(DEgenes[[i]]) > maxHM) {DEgenesHM <- DEgenes[[i]][1:maxHM,]} else {DEgenesHM <- DEgenes[[i]]}
        
        plotmatrix <- dataGEXMTSet[rownames(dataGEXMTSet) %in% rownames(DEgenesHM), ,drop=F]
        
        # If Symbols are available, rows are annotated with symbols instead of probe IDs   
        indexRownames <- match(rownames(plotmatrix), rownames(featureGEXMTSet))

        if(!is.null(Symbol.column)) {
          rownames(plotmatrix) <- ifelse(featureGEXMTSet[indexRownames,Symbol.column] != "" & !is.na(featureGEXMTSet[indexRownames,Symbol.column]), 
                                         as.character(featureGEXMTSet[indexRownames,Symbol.column]), rownames(plotmatrix))
          }
        
        
        #  if HMincludeRelevantSamplesOnly=T, only samples considered belonging to the respective group comparison 
        if (HMincludeRelevantSamplesOnly) {
          groups2plot <- unlist(strsplit(comparisons[i], "-") )  
          groups2plot <- unique(sub("[(|)]", "", groups2plot))  # for comparison of comparison
          sampleTable <- phenoGEXMTSet[,c(sampleColumn,groupColumn)]
          samples2plot <- character()
          
          for (j in 1:length(groups2plot)) {
            samples2plot <- c(samples2plot, as.character(sampleTable[sampleTable[,groupColumn]==groups2plot[j], sampleColumn])) }
          
          cols2plot <- colnames(plotmatrix) %in% samples2plot
          plotmatrix <- plotmatrix[, cols2plot]
        }
        
        pdf(file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, "_", comparisons[i], ".pdf", sep="" )), width = 10, height = 14) 
        heatmap.2(plotmatrix, main=comparisons[i], margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both",  
                  cexRow=HMcexRow, cexCol=HMcexCol, trace="none", density.info="histogram",  
                  key.xlab=plot.label, key.ylab="", key.title="Color Key", col=color.palette)  
        
        # Settings for plotting color key instead of dendrogram on top of heatmap:  
        # lmat=rbind(c(0,3),c(0,4), c(2,1)), lhei = c(0.3,0.5,3.8), lwid = c(0.5,4), key.par=list(mar=c(4,2,2,13)), Rowv=TRUE, Colv=FALSE, dendrogram="row", keysize=0.8 
        dev.off()
      } # end if nrow(DEgenes[[i]]) >1
   
} # end i-loop for group comparisons
     

  
  ###### Heatmaps with selected group Foldchanges 
  # (use 'DiffAllGroups_Ftest' and 'resultDiffGenes' from above)
  # DiffAllGroups_Ftest contains all probes with at least one absolute log-fold-changes greater than lfc. 
  if(!is.null(FC.heatmap.comparisons)) {

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
            indexRownames <- match(rownames(plotmatrix), rownames(featureGEXMTSet))
            if(!is.null(Symbol.column)) {
              rownames(plotmatrix) <- ifelse(featureGEXMTSet[indexRownames,Symbol.column] != "" & !is.na(featureGEXMTSet[indexRownames,Symbol.column]),
                                             as.character(featureGEXMTSet[indexRownames,Symbol.column]), rownames(plotmatrix))
            }

            if(nrow(plotmatrix) >1) {  # no Heatmap if just one diff expressed gene
              if (nrow(plotmatrix) > maxHM) {plotmatrix <- plotmatrix[1:maxHM,]} 
              if(is.null(names(FC.heatmap.comparisons)[fc])) {names(FC.heatmap.comparisons)[fc] <- paste0("set",fc)}
              filename.FCHM <- file.path(projectfolder, "Heatmaps", paste0("Heatmap_Foldchanges_", projectname, names(FC.heatmap.comparisons)[fc], "_", g, "_prior_by_", f, ".pdf"))
              cat("\nWrite Heatmap to", filename.FCHM, "\n")  
              
              pdf(filename.FCHM, width = 10, height = 14) 
              heatmap.2(plotmatrix, main=paste0("Heatmap Foldchanges ", projectname, names(FC.heatmap.comparisons)[fc], " ", g, "\nprobes prioritised by ", f), 
                        margins = c(15, 15), Rowv=TRUE, Colv=TRUE, dendrogram="both",  
                        cexRow=HMcexRow, cexCol=HMcexCol, trace="none", density.info="histogram",  
                        key.xlab="log2(FC)", key.ylab="", key.title="Color Key", col=color.palette)  
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
return(DEgenes.unfilt)


}


