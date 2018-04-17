#' Differential expression analysis with DESeq2
#' 
#' Differential Gene expression analysis for multiple group comparisons incl. output of 
#' Venn diagrams, volcano plots and heatmaps.
#' 
#' Function uses the \code{DESeq2} package for differential expression analysis of an 
#' un-normalized, un-transformed DESeqDataSet.
#' Analysis is performed for any number of group comparisons. The group designations in 
#' the \code{groupColumn} must match the contrasts given in \code{comparisons} (\code{"groupA-groupB"}).
#' P-value and foldchange thresholds may be applied to filter result data. 
#' Volcano plots are generated for each group comparison given in \code{comparisons} with highligted 
#' significance thresholds.
#' For generation of heatmaps, the count data is log2-transformed. Heatmaps are generated for top 
#' differentially expressed genes for each group comparison. 
#'  
#' @param dds DESeqDataSet
#' @param min_rowsum numeric. Minimum rowsum of countmatrix of \code{dds}. All rows with \code{rowSums < min_rowsum}
#' are removed from the count matrix. 
#' @param comparisons character vector with group comparisons in format \code{"groupA-groupB"} 
#' @param p.value.threshold numeric p-value threshold 
#' @param adjust.method adjustment method for multiple testing (\code{"none", "BH", "BY" and "holm"})
#' @param FC.threshold numeric foldchange threshold. 
#' 
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' 
#' @param Symbol.column character with column name of Gene Symbols in \code{dds} or NULL.
#' @param sampleColumn character with column name of Sample names in \code{dds}
#' @param groupColumn character with column name of group names in \code{dds}. Group names must match \code{comparisons}!
#' @param venn_comparisons character vector or (named) list of character vectors with group comparisons to be included 
#'                   in Venn Diagram (max 5 comparisons per Venn diagramm). Must be subset of \code{comparisons} or NULL.
#'
#' @param maxHM (numeric) max number of top diff. regulated elements printed in Heatmap
#' @param scale character indicating if the values should be centered and scaled in either the row direction 
#' or the column direction, or none (default). Either \code{"none"}, \code{"row"} or \code{"column"}.
#' @param HMcexRow labelsize for rownames in Heatmap
#' @param HMcexCol labelsize for colnames in Heatmap
#' @param figure.res numeric resolution for png.
#' @param HMincludeRelevantSamplesOnly (boolean) if TRUE include only Samples in heatmap, 
#'                                which are in groups of the respective group comparison.  
#' @param color.palette select color palette for heatmaps                     
#'                           
#'
#' @return List of elements 
#' \itemize{
#'   \item DEgenes list of dataframes with differential expressed genes according to significance 
#'   thresholds for each comparison (sorted by p-value). List elements are named by the respective group comparison. 
#' \item DEgenes.unfilt list of dataframes with unfiltered differential expression data for each comparison (sorted by p-value).
#' }
#' Filtered and unfiltered result tables, heatmaps, MA-plots and volcano plots for each group comparison 
#' as well as dispersion plot and Venn diagrams are stored as side-effects in the project folder.
#' 
#' @author Frank Ruehle
#' 
#' @export wrapDESeq



wrapDESeq <- function(dds, 
                  min_rowsum = 10,
                  comparisons,
                  p.value.threshold = 0.05, 
                  adjust.method="BH", 
                  FC.threshold = log2(1.5),
                      
                  projectfolder = file.path("GEX/deseq"),
                  projectname = "", 
                      
                  Symbol.column = NULL,
                  sampleColumn = "Sample_Name",   
                  groupColumn= "Sample_Group", 
                  venn_comparisons= NULL, 

                  # Heatmap parameter:
                  maxHM=50, 
                  scale = c("none"),
                  HMcexRow=1.1, 
                  HMcexCol=1.1, 
                  figure.res = 300,
                  HMincludeRelevantSamplesOnly=TRUE,
                  color.palette="heat.colors"
                  ) {


  
  # load required packages. Packages are not detached afterwards
  pkg.cran <- c("gplots", "VennDiagram")
  pkg.bioc <- c("DESeq2", "limma", "Biobase")
  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  projectname <- if (!is.null(projectname) && projectname!=""  && !grepl("_$", projectname)) {paste0(projectname, "_")} else {""}
  
  
  # Create output directories if not yet existing 
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
  if (!file.exists(file.path(projectfolder, "deseq_unfiltered"))) {dir.create(file.path(projectfolder, "deseq_unfiltered")) }
  if (!file.exists(file.path(projectfolder, "deseq_filtered"))) {dir.create(file.path(projectfolder, "deseq_filtered")) }
  if (!file.exists(file.path(projectfolder, "Heatmaps"))) {dir.create(file.path(projectfolder, "Heatmaps")) }
  if (!file.exists(file.path(projectfolder, "Volcano_plots"))) {dir.create(file.path(projectfolder, "Volcano_plots")) }
  if (!file.exists(file.path(projectfolder, "MA_plots"))) {dir.create(file.path(projectfolder, "MA_plots")) }
  
  
 
   # modify plot labels (heatmap legend)
  heatmap_legend_xaxis <- "log2(counts)"
  if (grepl("row", scale)) {heatmap_legend_xaxis <- "row z-score"}
  if (grepl("col", scale)) {heatmap_legend_xaxis <- "column z-score"}

  
  
  # if no thresholds given for p-value and foldchange, values are applied to omit filtering
  if(is.null(p.value.threshold)) {p.value.threshold <- 1}
  if(is.null(FC.threshold)) {FC.threshold <- 0}
  
  
  # Pre-filtering
  # not necessary, but reduces memory size. Note that more strict filtering to increase power is 
  # automatically applied via independent filtering on the mean of normalized counts within the results function.
  cat("\nApply pre-filtering for rowSums >=", min_rowsum)
  dds <- dds[rowSums(counts(dds)) >= min_rowsum,]
  
  
  # Differential expression
  # also generates sizefactors for library size normalisation used in results function.
  cat("\nDifferential expression analysis with DESeq2.")
  dds <- DESeq(dds)
  
  
  if (class(dds) == "DESeqDataSet") {
    cat("\nUsing variance Stabilizing Transformation for generating heatmaps")
    vsd <- DESeq2::varianceStabilizingTransformation(dds) # includes normalisation for library size
    # expmatrix <- DESeq2::rlog(dds, fitType="local") # class: DESeqTransform
    expmatrix <- assay(vsd)
    features <- mcols(dds,use.names=TRUE)
    # colData(dds) # sample phenotypes
  }
  
  
  

######### differential group comparisons
res.unfilt <- list()
table_unfilt <- list()
table_filt <- list()

# loop for all dedicated group comparisons
for (i in 1:length(comparisons)) {
  contrast_i <- c(groupColumn, sub("-.*$", "", comparisons[i]), sub("^.*-", "", comparisons[i]))
  
  # Calculate differential expressed genes (filtered and unfiltered)


  res.unfilt[[comparisons[i]]] <- DESeq2::results(dds, lfcThreshold = 0,  alpha = 0.05, 
                                                      pAdjustMethod=adjust.method, contrast=contrast_i)
  table_unfilt[[comparisons[i]]] <- as.data.frame(res.unfilt[[comparisons[i]]])
  # remove NA entries
  cat("\nRemove", sum(is.na(table_unfilt[[comparisons[i]]][,"padj"])), "entries with missing p-values")
  table_unfilt[[comparisons[i]]] <- table_unfilt[[comparisons[i]]][!is.na(table_unfilt[[comparisons[i]]][,"padj"]),]
  # add column with rownames id
  table_unfilt[[comparisons[i]]] <- data.frame(id.rownames= rownames(table_unfilt[[comparisons[i]]]), table_unfilt[[comparisons[i]]])
  # add Symbol column
  if(!is.null(Symbol.column)) {
    if(!(Symbol.column %in% names(mcols(dds)))) {stop("Symbolcolumn not found in feature data!")}
    table_unfilt[[comparisons[i]]][, Symbol.column] <- mcols(dds)[match(rownames(table_unfilt[[comparisons[i]]]), rownames(dds)), Symbol.column]
  }
  # add ENTREZID column
  if(any(grepl("entrezid", names(mcols(dds)), ignore.case = T))) {
    entrezIdColumnname <- grep("entrezid", names(mcols(dds)), ignore.case = T, value = T)[1] # get name of (first) EntrezID column 
        table_unfilt[[comparisons[i]]][, entrezIdColumnname] <- mcols(dds)[match(rownames(table_unfilt[[comparisons[i]]]), rownames(dds)), entrezIdColumnname]
  }
  
  
  # filter for significance
  table_filt[[comparisons[i]]]  <- filterGeneLists(table_unfilt[[comparisons[i]]],
                                               newheader=NULL,
                                               filtercat1 = "padj",
                                               filtercat1.decreasing = FALSE,
                                               filtercat1.function = identity,
                                               filtercat1.threshold=p.value.threshold,
                                               filtercat2 = "log2FoldChange",
                                               filtercat2.decreasing = TRUE,
                                               filtercat2.function = abs,
                                               filtercat2.threshold = FC.threshold)
  

  # comment bioconductor workflow on http://www.bioconductor.org/help/workflows/rnaseqGene/
  # The column log2FoldChange is the effect size estimate. It tells us how much the gene's expression seems to 
  # have changed due to treatment with dexamethasone in comparison to untreated samples. This value is reported 
  # on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene's expression 
  # is increased by a multiplicative factor of 2**1.5=2.82.
  
  
  # Writing result gene tables
  write.table(table_unfilt[[comparisons[i]]], sep="\t", quote=F, row.names=F,
              file= file.path(projectfolder, "deseq_unfiltered", paste0(projectname, comparisons[i], "_unfilt.txt")))
  
  write.table(table_filt[[comparisons[i]]], sep="\t", quote=F, row.names=F,
              file= file.path(projectfolder, "deseq_filtered", paste0(projectname, comparisons[i], ".txt")))
  
  cat(paste("\n", nrow(table_filt[[comparisons[i]]]),"differentially regulated elements for comparison:",comparisons[i]))
  cat(paste("\nWrite gene tables to", file.path(projectfolder, "deseq_unfiltered", paste0(projectname, comparisons[i], "_unfilt.txt")),
            "and", file.path(projectfolder, "deseq_filtered", paste0(projectname, comparisons[i], ".txt"))))
  
  
 
  ######## Volcano plots per group comparison 
  
  volcanoplot_filename <- file.path(projectfolder, "Volcano_plots", paste("Volcano_", projectname, comparisons[i], ".png", sep="" ))
  
  cat("\nWrite Volcano plot to", volcanoplot_filename, "\n")  
  
  png(file=volcanoplot_filename, width = 150, height = 150, units = "mm", res=figure.res) 
  plot(table_unfilt[[comparisons[i]]]$log2FoldChange, -log10(table_unfilt[[comparisons[i]]]$padj), main=comparisons[i],
       col="gray", pch=16, cex=0.5, xlab="fold change", ylab="-log10 p-value") 
  if(!is.null(p.value.threshold)) {abline(h=-log10(p.value.threshold),lty=c(1))}
  if(!is.null(FC.threshold)) {abline(v=c(FC.threshold, -FC.threshold),lty=2)}
  
  if(nrow(table_filt[[comparisons[i]]]) >=1) {  # no highlighting if no diff expressed genes
    points(table_filt[[comparisons[i]]]$log2FoldChange, -log10(table_filt[[comparisons[i]]]$padj),pch=16, col="black", cex=0.8) 
  } 
  dev.off()
  
  
  
  
 
  ######## Heatmaps per group comparison with signal intensities
      if(nrow(table_filt[[comparisons[i]]]) >1) {  # no Heatmap if just one diff expressed gene
        
        cat("\nWrite Heatmap to", file.path(projectfolder, "Heatmaps", paste("Heatmap_", projectname, comparisons[i], ".png", sep="" )), "\n")  
        
        if (nrow(table_filt[[comparisons[i]]]) > maxHM) {DEgenesHM <- table_filt[[comparisons[i]]][1:maxHM,]} else {DEgenesHM <- table_filt[[comparisons[i]]]}
 

        plotmatrix <- expmatrix[rownames(expmatrix) %in% rownames(DEgenesHM), ,drop=F]


        # If Symbols are available, rows are annotated with symbols instead of probe IDs   
        indexRownames <- match(rownames(plotmatrix), rownames(features))

        if(!is.null(Symbol.column)) {
          rownames(plotmatrix) <- ifelse(features[indexRownames,Symbol.column] != "" & !is.na(features[indexRownames,Symbol.column]), 
                                         as.character(features[indexRownames,Symbol.column]), rownames(plotmatrix))
          }
        
        groupColorCode <- rainbow(length(unique(colData(dds)[,groupColumn])))[as.numeric(factor(colData(dds)[,groupColumn]))] # for ColSideColors in heatmaps
        
        #  if HMincludeRelevantSamplesOnly=T, only samples considered belonging to the respective group comparison 
        if (HMincludeRelevantSamplesOnly) {
          groups2plot <- unlist(strsplit(comparisons[i], "-") )  
          groups2plot <- unique(sub("[(|)]", "", groups2plot))  # for comparison of comparison
          sampleTable <- colData(dds)[,c(sampleColumn,groupColumn)]
          samples2plot <- character()
          
          for (j in 1:length(groups2plot)) {
            samples2plot <- c(samples2plot, as.character(sampleTable[sampleTable[,groupColumn]==groups2plot[j], sampleColumn])) }
          
          cols2plot <- colnames(plotmatrix) %in% samples2plot
          plotmatrix <- plotmatrix[, cols2plot]
       
          groupColorCode <- groupColorCode[colData(dds)[,groupColumn] %in% groups2plot] # adjust ColSideColors in heatmaps to relevant groups
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
   

  # MA-plot
  filename.MA <- file.path(projectfolder, "MA_plots", paste0("MA_plot_", projectname, comparisons[i], ".png"))
  cat("\nWrite MA-plot to", filename.MA)
  png(filename= filename.MA, width = 150, height = 150, units = "mm", res=figure.res)
  DESeq2::plotMA(res.unfilt[[comparisons[i]]], ylim=c(-2,2), alpha = 0.05) # 
  dev.off()
  
  
  
  
  } # end i-loop for group comparisons
     

### Venn Diagram
if(is.null(venn_comparisons)) {venn_comparisons <- comparisons}
  if(!is.list(venn_comparisons)) {venn_comparisons <- list(venn_comparisons)}
  for(v in 1:length(venn_comparisons)) {
    
    if(is.null(names(venn_comparisons)[v])) {names(venn_comparisons)[v] <- paste0("Vennset",v)}
    filename.Venn <- file.path(projectfolder, paste0("Venn_Diagram_", projectname, names(venn_comparisons)[v], ".png"))
    cat("\nWrite Venn diagramm to", filename.Venn, "\n")
    
    if (length(venn_comparisons[[v]]) <= 5) {

      venndata <- table_filt[which(names(table_filt) %in% venn_comparisons[[v]])]
      category.names <- names(venndata) # insert \n in venn labels if character length >15
      category.names <- ifelse(nchar (category.names)>15, sub("-", "-\n", category.names), category.names)
      for (i in names(venndata)) {venndata[[i]] <- rownames(venndata[[i]])}
      png(filename= filename.Venn, width = 150, height = 150, units = "mm", res=figure.res)
      venn.plot <- venn.diagram(venndata, filename= NULL, category.names= category.names, cex = 1, cat.cex=1, alpha=0.3, margin=0.3)
      grid.draw(venn.plot)
      dev.off()
    }
  } # end v-loop


# Dispersion plot
filename.Disp <- file.path(projectfolder, paste0(projectname, "Dispersion_plot.png"))
cat("\nWrite Dispersion plot to", filename.Disp)
png(filename= filename.Disp, width = 150, height = 150, units = "mm", res=figure.res)
DESeq2::plotDispEsts(dds, main="Dispersion Estimates")
dev.off()


# pca plot
filename.pca <- file.path(projectfolder, paste0(projectname, "pca_plot.png"))
cat("\nWrite pca plot to", filename.pca)
png(filename= filename.pca, width = 150, height = 150, units = "mm", res=figure.res)
DESeq2::plotPCA(vsd, intgroup= groupColumn)
dev.off()



# # plot counts: examine the counts of reads for a single gene across the groups
# plotCounts(dds, gene=which.min(res$padj), intgroup=groupColumn)

# return unfiltered genes as list of dataframes for each comparison (sorted by p-value)
return(list(DEgenes=table_filt, DEgenes.unfilt=table_unfilt))

}


