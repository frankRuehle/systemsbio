
#' Wrapper for gene enrichment analysis using clusterProfiler
#' 
#' Applying overrepresentation analysis and gene set enrichment analysis to supplied gene list
#' 
#' Function uses genelist and optionally background gene list as input to perform enrichment analysis 
#' using the \code{clusterProfiler}-package. By default, overrepresentation analysis and gene set enrichment analysis 
#' (GSEA) is performed for all categories given in \code{enrichmentCat}.
#' The Gene list can be supplied as vector with ids given in \code{id.type},
#' dataframe with IDs of type \code{id.type} given in column \code{id.column} or character string with path of data file.
#' The IDs are converted to ENTREZ IDs (if necessary) prior to enrichment using the annotation package
#' for the species denoted in \code{org}.
#' Optionally, quanitative data can be included in \code{sortcolumn} for sorting and filtering (using \code{sortcolumn.threshold}) 
#' the data. If quantitative data are no fold changes, fold changes may be given additionally in \code{FCcolumn}.
#' They are used only for cnetplot of enrichment results only. If Null the data in \code{sortcolumn} is also used for 
#' the cnetplots. Be aware that the legend of the cnetplots will be "Fold Change" anyway!
#' If no quantitative data is provided, sorting and GSEA will be skipped.
#' Optionally, a background list can be provided for enrichment analysis.   
#'
#' 
#' @param genes vector with gene names to analyse or dataframe or character with path to dataframe. Supplied
#'        dataframe needs to include a column with identifiers specified in 'id.column'.
#' @param newheader optional character vector with new header information for 'genes' dataframe. Only relevant 
#'            if 'genes' is a dataframe (or character string with filepath to a table) with wrong or missing header. 
#'            NULL otherwise.
#' @param backgroundlist optional background list for enrichment analysis. Can either be a vector with gene IDs 
#'                 or dataframe or character with path to dataframe. If given, an 'id.column' needed as in 'genes'. 
#'                 If "genome", all ENTREZ IDs from the annotation package of the respective org
#'                 (denoted in 'org') are used as background.
#'                 If NULL, full ID list from 'genes' is used as background.
#' @param newheaderBackground optional character vector with new header information for 'backgroundlist' dataframe.
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param enrichmentCat character vector with categories to be enriched (GO: gene ontology (MF, BP, CC), 
#'                KEGG: KEGG pathways, Reactome: Reactome pathways, DO: Disease ontology).
#' @param maxInputGenes (numeric) max number of top diff regulated elements used for enrichment analysis.(or NULL).
#' @param id.type character with identifier type from annotation package (e.g. "ENTREZID", "SYMBOL", "UNIGENE")
#'          IDs Will be converted to EntrezIDs prior to enrichment analysis.
#' @param id.column character with column name for identifier variable in 'genes'. 
#' @param sortcolumn character with column name of quantitative data in 'genes' used for ordering.
#'             If Null, ranking of genes is omitted and GSEA not possible.
#' @param sortdecreasing (boolean) order parameter for hierarchy of values in 'sortcolumn'.
#'                 FALSE for increasing values (e.g. p-values), TRUE for decreasing values (e.g. fold changes).
#'                 If FALSE, values in sortcolumn will be transformed prior to GSEA.
#' @param sortcolumn.threshold numeric threshold for 'sortcolumn' to be included in overepresentation analysis.
#'                       If sortdecreasing=F, value < sortcolumn.threshold 
#'                       else value > sortcolumn.threshold
#' @param fun.transf.incr.vales function definition for transforming values in 'sortcolumn' prior to GSEA.
#'                        Is applied if sortdecreasing=F.
#' @param FCcolumn (character) optional column name of foldchanges in 'genes' if 'sortcolumn' is used elsewhere.
#'           Used only for cnetplot of enrichment results. Omitted if NULL
#' @param threshold_FC (numeric) Fold change threshold for filtering (threshold interpreted for log transformed foldchange values!)
#'               Only relevant for overrepresentation analysis if an unfiltered gene list is given in 'genes'
#'               to allow for parallel GSEA.
#' @param pAdjustMethod method for adjusting for multiple testing. 
#'                One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param org character with name of organism ("human", "mouse", "rat").
#' @param enrich.p.valueCutoff numeric p-value threshold for returned enrichment terms.
#' @param enrich.q.valueCutoff numeric q-value threshold for returned enrichment terms.
#' @param nPerm permutation numbers
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' 
#' 
#' @return List of enrichment-objects defined in \code{DOSE}-package (\code{enrichResult}-object for overrepresentation 
#'   analysis, \code{gseaResult}-objects for gene set enrichment analysis).                          
#'   Enrichment tables and plots are stored in the project folder as side effects.
#' 
#' @author Frank Ruehle
#' 
#' @export 

    
## Usage 
clusterprof <- function (genes, 
                         newheader = NULL, 
                         backgroundlist=NULL, 
                         newheaderBackground = NULL,
                         projectfolder= "GEX/clusterProfiler",
                         projectname="", 
                         enrichmentCat = c("GO", "KEGG", "Reactome", "DO"),
                         maxInputGenes = 100,  
                         id.type = "SYMBOL",
                         id.column = "SYMBOL",   
                         sortcolumn ="adj.P.Val",  
                         sortdecreasing = FALSE, 
                         sortcolumn.threshold = 0.05,
                         fun.transf.incr.vales = function(x) {-log10(x)},
                         FCcolumn = "logFC",
                         threshold_FC= log2(1.5),
                         org = "human",
                         pAdjustMethod = "BH", 
                         enrich.p.valueCutoff = 0.05, 
                         enrich.q.valueCutoff = 0.05,
                         nPerm        = 1000,
                         minGSSize    = 10) 
{
  
  
    
  ## create output directory
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T)} 
    
  
  # define annotation database for organism
  annotationdb <- switch(org, human = "org.Hs.eg.db",
                                  mouse = "org.Mm.eg.db",
                                  rat= "org.Rn.eg.db")
  
  
  
  ## install/load required packages from CRAN and Bioconductor
  pkg.bioc <- c("clusterProfiler", "DOSE", "ReactomePA", annotationdb)
  pkg.cran <- c("plyr")
  attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  if(!is.null(id.column)) {
  if(id.column!="ENTREZID" && "ENTREZID" %in% names(genes)) {
    warning("Column ENTREZIDs found. Consider these data for id.column")
    }
  }
  
  
  ### process data input
  # read file if 'genes' is character string with file path
  if(is.character(genes) & length(genes)==1) {
    cat("\n\nReading gene list:", genes, "\n")
    genes <- read.table(genes, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
  }
  
  if(is.data.frame(genes) & !is.null(newheader)) {# if 'newheader' is defined, it is used as names(DEgenes.unfilt)
    cat("\nNew header added to input file:", newheader, "\n")
    colnames(genes) <- newheader
  }
  
  columns2use <- c(id.column, sortcolumn, FCcolumn) 
  columnsfound <- columns2use[columns2use %in% names(genes)]
  if (!all(columns2use %in% names(genes))) {
    warning("column(s) ", columns2use[!(columns2use %in% names(genes))], " not found in data header.")
  }  
  genes <- genes[,columnsfound, drop=F]  # reduce dataframe to columns of interest
    
  if(is.vector(genes)) {
    if(!is.null(names(genes))) { # if named vector with IDs as element names and quantitative information as elements
      genes <- data.frame(names(genes), genes)
      names(genes) <- c(id.column, sortcolumn)        
        } else { # if unnamed vector with IDs but no quantitative information (no sortcolumn)
            genes <- data.frame(genes)
            names(genes) <- c(id.column)        
            }
        }
  
  # now 'genes' is a data.frame with 1, 2 or 3 columns (id.column, sortcolumn, FCcolumn)
      genes <- genes[!is.na(genes[,id.column]),] # remove rows without entry in id.column
      if (!is.null(sortcolumn)) {
      if (sortcolumn %in% names(genes)) {    
        cat(paste("\nSorting input dataframe for", sortcolumn, ", decreasing =", sortdecreasing, "\n"))
        genes <- genes[order(genes[,sortcolumn], decreasing=sortdecreasing),] # sorting genes for quantitative variable
        
        # if multiple probes per gene, only the first gene entry is used 
        genes <- genes[!duplicated(genes[,id.column]),] 
        }}
     
      
      
  ## convert IDs to ENTREZ IDs if necessary
  ### bitr noch ersetzen durch basic_anno!
    if(id.type!="ENTREZID") {
    cat(paste("\nConverting", id.type, "from input list to ENTREZIDs\n"))
    entrezids <- bitr(genes[,id.column], fromType=id.type, toType="ENTREZID", annoDb=annotationdb)
    cat(paste(nrow(entrezids), " of ",  length(unique(genes[,id.column])), "unique", id.type, "mapped to ENTREZIDs\n"))
    cat(paste("For dublicated", id.type, "only the first element is used (dataset was ordered for", sortcolumn, "decreasing =", sortdecreasing, ")\n")) 
    genes <- plyr::join(genes, entrezids, by=id.column, type="right", match="first") 
    # joined dataframe contains only entries with available ENTREZIDs, but keeps order of 'genes'
  } else {names(genes)[names(genes)==id.column] <- "ENTREZID"}
  
   # 'genes' now definitively provides a column 'ENTREZID'   
   # if multiple probes per gene, only the first gene entry is used 
   # (after sorting, if sortcolumn available).
   genes <- genes[!duplicated(genes[,"ENTREZID"]),] 
      
      
  
  ## Threshold applied to gene list: 'maxInputGenes' as well as 'sortcolumn.threshold' if applicable
  # filtered genes
  if (!is.null(sortcolumn) & !is.null(sortcolumn.threshold)) {
    if (sortcolumn %in% names(genes)) {   
        if(sortdecreasing==F) { # for increasing sortcolumn (e.g. p-values) value < threshold
          filtgenes <- genes[genes[,sortcolumn] < sortcolumn.threshold, ]
          cat(paste("\nGenelist filtered for", sortcolumn, "<", sortcolumn.threshold, "(", nrow(filtgenes), "genes remaining)\n"))
          
        } else { # for decreasing sortcolumn (e.g. fold changes) value > threshold
          filtgenes <- genes[genes[,sortcolumn] > sortcolumn.threshold, ]
          cat(paste("\nGenelist filtered for", sortcolumn, ">", sortcolumn.threshold, "(", nrow(filtgenes), "genes remaining)\n"))
        }
        if(!is.null(FCcolumn) & (FCcolumn!=sortcolumn) & !is.null(threshold_FC)) { # additional filtering for foldchange after filtering sortcolumn
          filtgenes <- filtgenes[abs(filtgenes[,FCcolumn]) > threshold_FC, ]
          cat(paste("\nGenelist filtered for abs(", FCcolumn, ") >", round(threshold_FC, digits=4), "(", nrow(filtgenes), "genes remaining)\n"))
        }
      
    } else {filtgenes <-  genes}
  } else {filtgenes <-  genes}
  
  cat(paste("\n", min(nrow(filtgenes), maxInputGenes), "unique genes for overrepresentation analysis.\n"))  
  filtgenes <- filtgenes[1:min(nrow(filtgenes), maxInputGenes),]
  
  
  ## Processing list for Gene Set Enrichment Analysis (GSEA). Ordered named Vector needed.
  if (!is.null(sortcolumn)) { # quantitative data available?
    gseagenes <- genes[,sortcolumn] 
    names(gseagenes) <- genes[,"ENTREZID"]
    # for increasing sortcolumn (e.g. p-values). Low values transformed to high values
    if(sortdecreasing==F) {gseagenes <- fun.transf.incr.vales(gseagenes)} 
    gseagenes <- sort(gseagenes, decreasing=T) # Vector in decreasing order.
    cat(paste("\n", length(gseagenes), "unique genes for gene set enrichment analysis.\n"))  
  }
  
  
  
  
  
  
  ###### read background list if available. Otherwise unfiltered ID list from 'genes' used as background Vector.
  if(!is.null(backgroundlist)) {
    
    if (all(backgroundlist=="genome")) {
      # all ENTREZIDs from the annotation db of the respective genome (given in org) used as background.
      backgroundlist <- keys(get(annotationdb), keytype="ENTREZID")
      cat(paste("\nUsing", length(backgroundlist), "ENTREZ IDs from",annotationdb, "as background.\n"))      
      } else {
    
    # read file if 'backgroundlist' is character string with file path
    if(is.character(backgroundlist) & length(backgroundlist)==1) {
      cat("\n\nReading background list:", backgroundlist, "\n")
      backgroundlist <- read.table(backgroundlist, header=is.null(newheaderBackground), sep="\t", na.strings = c("", " ", "NA")) # if no newheaderBackground, header must be in file
    }
    if(is.data.frame(backgroundlist) & !is.null(newheaderBackground)) {# if 'newheaderBackground' is defined, it is used as names(DEgenes.unfilt)
      cat("\nNew header added to background list:", newheaderBackground, "\n")
      colnames(backgroundlist) <- newheaderBackground
      }
    
    if(is.vector(backgroundlist)) {
      backgroundlist <- data.frame(unique(backgroundlist))
        } 
    
    # remove rows without entry in id.column and make vector with id.column     
    backgroundlist <- as.character(backgroundlist[!is.na(backgroundlist[,id.column]),id.column]) # remove rows without entry in id.column
    backgroundlist <- unique(backgroundlist)
  
    
    # convert IDs to ENTREZ IDs if necessary
    if(id.type!="ENTREZID") {
      cat(paste("\nConverting", id.type, "from background list to ENTREZIDs\n"))
      backgroundlist <- bitr(backgroundlist, fromType=id.type, toType="ENTREZID", annoDb=annotationdb)
      backgroundlist <- as.character(unique(backgroundlist[,"ENTREZID"]))
    } 
     
    cat(paste("\n", length(backgroundlist), "ENTREZ IDs in background list.\n"))
      
    
     } # end loading and processing of a supplied backgroundlist 
    
    } else { # if backgroundlist is NULL, all EntrezIDs from 'genes' used as background instead.
      backgroundlist <- as.character(unique(genes[,"ENTREZID"]))
      cat(paste("\n", length(backgroundlist), "ENTREZ IDs used from input object as background.\n"))
    }
  
  

  

  #################### Start Enrichment analysis
  
  # remove categories not supported by clusterProfile
  enrichmentCat.found <- enrichmentCat %in% c("GO", "KEGG", "Reactome", "DO")
  enrichmentCat <- enrichmentCat[enrichmentCat.found]
  cat(paste("\nCategories used for enrichment analysis with clusterProfiler:", paste(enrichmentCat, collapse=" "), "\n"))
  
  # initialise result list
  enrichresult <- list()
  
  ##### Gene Ontology Enrichment
  ## Over representation analysis

  if ("GO" %in% enrichmentCat) {

  for (go in c("MF", "BP", "CC")) {
    if(nrow(filtgenes)>1) {
      enrichresult[[paste0("Overrep_",go)]] <- enrichGO(gene = filtgenes[,"ENTREZID"],
                                    universe      = backgroundlist,
                                    organism      = org,
                                    ont           = go,
                                    minGSSize    = minGSSize,
                                    pAdjustMethod = pAdjustMethod, 
                                    pvalueCutoff  = enrich.p.valueCutoff,
                                    qvalueCutoff  = enrich.q.valueCutoff,
                                    readable      = TRUE) # EntrezIDs displayed as gene Symbols
      cat(paste("Overrepresentation analysis for", go, "GO-terms results in", nrow(summary(enrichresult[[paste0("Overrep_",go)]])), "enriched terms.\n"))
    }


  # GO Gene Set Enrichment Analysis
    if (!is.null(sortcolumn)) { # quantitative data available?
      cat(paste("GeneSet Enrichment analysis for", go, "GO-terms results in "))
      enrichresult[[paste0("GSEA_",go)]] <- gseGO(geneList = gseagenes,
                                                  organism     = org,
                                                  ont          = go,
                                                  nPerm        = nPerm,
                                                  minGSSize    = minGSSize,
                                                  pvalueCutoff = enrich.p.valueCutoff,
                                                  pAdjustMethod = pAdjustMethod, 
                                                  verbose      = FALSE)
  cat(paste(nrow(summary(enrichresult[[paste0("GSEA_",go)]])), "enriched terms.\n"))
  }
  
  } # end go-loop
} # end if GO


##### KEGG pathway enrichment
# KEGG over-representation test

  if ("KEGG" %in% enrichmentCat) {
    if(nrow(filtgenes)>1) {
      enrichresult[[paste0("Overrep_","KEGG")]] <- enrichKEGG(gene = filtgenes[,"ENTREZID"],
                     organism     = org,
                     universe    = backgroundlist, 
                     pvalueCutoff = enrich.p.valueCutoff, 
                     qvalueCutoff  = enrich.q.valueCutoff,
                     pAdjustMethod = pAdjustMethod, 
                     minGSSize    = minGSSize,
                     readable     = TRUE,
                     use_internal_data = FALSE)
    
    cat(paste("Overrepresentation analysis for KEGG pathways results in", nrow(summary(enrichresult[[paste0("Overrep_","KEGG")]])), "enriched pathways\n"))
    }

# KEGG Gene Set Enrichment Analysis
if (!is.null(sortcolumn)) { # quantitative data available?
  cat(paste("GeneSet Enrichment analysis for KEGG pathways results in "))
  enrichresult[[paste0("GSEA_","KEGG")]] <- gseKEGG(geneList  = gseagenes,
                 organism     = org,
                 nPerm        = nPerm,
                 minGSSize    = minGSSize,
                 pvalueCutoff = enrich.p.valueCutoff,
                 pAdjustMethod = pAdjustMethod, 
                 verbose      = FALSE,
                 use_internal_data = FALSE)
  
  cat(paste(nrow(summary(enrichresult[[paste0("GSEA_","KEGG")]])), "enriched pathways\n"))
  }
}


##### Reactome pathway enrichment
# Reactome over-representation test
  
if ("Reactome" %in% enrichmentCat) {
  if(nrow(filtgenes)>1) {
    enrichresult[[paste0("Overrep_","Reactome")]] <- enrichPathway(gene =  filtgenes[,"ENTREZID"],
                      organism     = org,
                      universe    = backgroundlist, 
                      pvalueCutoff = enrich.p.valueCutoff, 
                      qvalueCutoff  = enrich.q.valueCutoff,
                      pAdjustMethod = pAdjustMethod, 
                      minGSSize    = minGSSize,
                      readable     = TRUE)
  
  cat(paste("Overrepresentation analysis for Reactome pathways results in", nrow(summary(enrichresult[[paste0("Overrep_","Reactome")]])), "enriched pathways\n"))
  }



# Reactome Gene Set Enrichment Analysis
  if (!is.null(sortcolumn)) { # quantitative data available?
    cat(paste("GeneSet Enrichment analysis for Reactome pathways results in "))
    enrichresult[[paste0("GSEA_","Reactome")]] <- gsePathway(geneList  = gseagenes,
                      organism     = org,
                      nPerm        = nPerm,
                      minGSSize    = minGSSize,
                      pvalueCutoff = enrich.p.valueCutoff,
                      pAdjustMethod = pAdjustMethod, 
                      verbose      = FALSE)
  
  cat(paste(nrow(summary(enrichresult[[paste0("GSEA_","Reactome")]])), "enriched pathways\n"))
  }

}


  ##### Disease Ontology enrichment
  # DO over-representation test
  if ("DO" %in% enrichmentCat) {
    if(nrow(filtgenes)>1) {
      enrichresult[[paste0("Overrep_","DO")]] <- enrichDO(gene =  filtgenes[,"ENTREZID"],
                    ont="DO",
                    universe    = backgroundlist, 
                    pvalueCutoff = enrich.p.valueCutoff, 
                    qvalueCutoff  = enrich.q.valueCutoff,
                    pAdjustMethod = pAdjustMethod, 
                    minGSSize    = minGSSize,
                    readable     = TRUE)
    
    cat(paste("Overrepresentation analysis for Disease Ontology results in", nrow(summary(enrichresult[[paste0("Overrep_","DO")]])), "enriched terms\n"))
    } 


  # DO Gene Set Enrichment Analysis
  if (!is.null(sortcolumn)) { # quantitative data available?
    cat(paste("GeneSet Enrichment analysis for Disease Ontology results in "))
    enrichresult[[paste0("GSEA_","DO")]] <- gseAnalyzer(geneList  = gseagenes,
                          setType = "DO",    
                          organism     = org,
                          nPerm        = nPerm,
                          minGSSize    = minGSSize,
                          pvalueCutoff = enrich.p.valueCutoff,
                          pAdjustMethod = pAdjustMethod, 
                          verbose      = FALSE)
  
  cat(paste(nrow(summary(enrichresult[[paste0("GSEA_","DO")]])), "enriched terms.\n"))
  }


}



# Creating result tables and plot 
cat(paste("\nCreating result tables and plots in", file.path(projectfolder), "\n"))

if(!is.null(projectname)) { # add underline to projectname if given
  if(projectname!="") {
  projectname <- paste0(projectname, "_")
  }}



for(usedcat in names(enrichresult)) {

    if(nrow(summary(enrichresult[[usedcat]])) >=1) {
  
    cat(paste("\nProcessing", usedcat))
    # Result table
    write.table(summary(enrichresult[[usedcat]]), quote=F, row.names=F, sep="\t",
                file=file.path(projectfolder, paste0(projectname, usedcat, "_resulttable.txt")))
  
    # Enrichment map
    tiff(file.path(projectfolder, paste0(projectname, usedcat, "_enrichmentMap.tiff")), 
         width = 7016 , height = 4960, res=600, compression = "lzw")
    try(enrichMap(enrichresult[[usedcat]], n = 30), silent=T)
    dev.off()
    
    # cnet plot (Remark: legend will be 'Fold change' even if other quantitative data are included)
    tiff(file.path(projectfolder, paste0(projectname, usedcat, "_cnetPlot.tiff")), 
         width = 7016 , height = 4960, res=600, compression = "lzw")
    try(cnetplot(enrichresult[[usedcat]], showCategory = 5, categorySize="pvalue", 
             foldChange= if(!is.null(FCcolumn)) {filtgenes[,FCcolumn]} else {filtgenes[,sortcolumn]}), silent=T)
    dev.off()

    }
}

cat("\n")

  
  # Detaching libraries not needed any more
# detach_package(c(pkg.cran, pkg.bioc))


return(enrichresult)

} # end of function definition 




