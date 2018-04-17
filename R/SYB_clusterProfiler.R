
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
#' This data is used for cnetplot of enrichment results as well as for KEGG pathway mapping of enriched 
#' KEGG pathways (if any) using the \code{pathview} package.
#' If \code{FCcolumn == NULL} the data in \code{sortcolumn} is also used for the cnetplots and pathway mapping. 
#' Be aware that the legend of the cnetplots will be "Fold Change" anyway.
#' If no quantitative data is provided, sorting and GSEA will be skipped.
#' Optionally, a background list can be provided for enrichment analysis.   
#'
#' 
#' @param genes vector with gene names to analyse or dataframe or character with path to dataframe or list thereof. 
#' Supplied dataframe needs to include a column with identifiers specified in \code{id.column}.
#' @param newheader optional character vector with new header information for \code{genes} dataframe. Only relevant 
#'            if 'genes' is a dataframe (or character string with filepath to a table) with wrong or missing header. 
#'            NULL otherwise.
#' @param backgroundlist optional background IDs for enrichment analysis. Can either be a vector with gene IDs 
#'                 or dataframe or character with path to dataframe. If given, an \code{id.column} is needed as 
#'                 in \code{genes} with type of IDs given in \code{id.type}. 
#'                 If "genome", all ENTREZ IDs from the annotation package of the respective organism
#'                 (denoted in \code{org}) are used as background.
#'                 If NULL, full ID list from \code{genes} is used as background (i.e. that \code{genes} need to 
#'                 contain all genes under investigation without pre-filtering for significance).
#' @param newheaderBackground optional character vector with new header information for \code{backgroundlist}.
#' @param projectfolder character with directory for output files (will be generated if not exisiting).
#' @param projectname optional character prefix for output file names.
#' @param enrichmentCat character vector with categories to be enriched (GO: gene ontology (MF, BP, CC), 
#'                KEGG: KEGG pathways, Reactome: Reactome pathways, DO: Disease ontology).
#' @param maxInputGenes (numeric) max number of top diff regulated elements used for enrichment analysis.(or NULL).
#' @param id.type character with identifier type from annotation package ("ENTREZID" or "SYMBOL")
#'          Gene symbols Will be converted to EntrezIDs prior to enrichment analysis.
#' @param id.column character with column name for identifier variable in \code{genes}. 
#' @param sortcolumn character with column name of quantitative data in \code{genes} used for ordering.
#'             If Null, ranking of genes is omitted and GSEA not possible.
#' @param sortdecreasing (boolean) order parameter for hierarchy of values in \code{sortcolumn}.
#'                 FALSE for increasing values (e.g. p-values), TRUE for decreasing values (e.g. fold changes).
#'                 If FALSE, values in sortcolumn will be transformed prior to GSEA.
#' @param sortcolumn.threshold numeric threshold for \code{sortcolumn} to be included in overepresentation analysis.
#'                       \code{If sortdecreasing=F, value < sortcolumn.threshold 
#'                       else value > sortcolumn.threshold}
#' @param fun.transf.incr.values function definition for transforming values in \code{sortcolumn} prior to GSEA.
#'                        Is applied if \code{sortdecreasing=F}.
#' @param FCcolumn (character) optional column name of foldchanges in \code{genes} if \code{sortcolumn} is used elsewhere.
#'           Used only for cnetplot of enrichment results. Omitted if NULL
#' @param threshold_FC (numeric) Fold change threshold for filtering (threshold interpreted for log transformed foldchange values!)
#'               Only relevant for overrepresentation analysis if an unfiltered gene list is given in \code{genes}
#'               to allow for parallel GSEA.
#' @param pAdjustMethod method for adjusting for multiple testing. 
#'                One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param org character with name of organism ("human", "mouse", "rat").
#' @param enrich.p.valueCutoff numeric p-value threshold for returned enrichment terms.
#' @param enrich.q.valueCutoff numeric q-value threshold for returned enrichment terms.
#' @param nPerm permutation numbers
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param figure.res numeric resolution for png.
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
                         id.type = "ENTREZID",
                         id.column = "ENTREZID",   
                         sortcolumn ="adj.P.Val",  
                         sortdecreasing = FALSE, 
                         sortcolumn.threshold = 0.05,
                         fun.transf.incr.values = function(x) {-log10(x)},
                         FCcolumn = "logFC",
                         threshold_FC= log2(1.5),
                         org = "human",
                         pAdjustMethod = "BH", 
                         enrich.p.valueCutoff = 0.05, 
                         enrich.q.valueCutoff = 0.05,
                         nPerm        = 1000,
                         minGSSize    = 10,
                         maxGSSize    = 500,
                         figure.res = 300) 
{
  
  
    
  ## create output directory
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T)} 
    
  if(!is.null(projectname)) { # add underline to projectname if given
    if(projectname!="") {
      projectname <- paste0(projectname, "_")
    }}
  
  
  
  # define annotation database for organism
  annotationdb <- switch(org, human = "org.Hs.eg.db",
                                  mouse = "org.Mm.eg.db",
                                  rat= "org.Rn.eg.db")
  
  org_alt_name <- switch(org, human = "hsa", # org parameter differs in GO and KEGG functions in clusterProfiler
                         mouse = "mmu",
                         rat= "rno")
  
  ## install/load required packages from CRAN and Bioconductor
  pkg.bioc <- c("clusterProfiler", "DOSE", "ReactomePA", "topGO", "pathview", annotationdb)
  pkg.cran <- c("plyr")
  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  
  if(!is.null(id.column)) {
    if(id.column!="ENTREZID" && "ENTREZID" %in% names(genes)) {
      warning("Column ENTREZIDs found. Consider these data for id.column.")
      }
  }
  
  
  ### process data input

  if(class(genes) != "list") { # make list if genes is single element
    genes <- list(genes=genes)
  }

  # initialise help variable to avoid multiple annotation of background list
  vector.processed <- FALSE
  
  # initialise result list
  enrichresult <- list()
  
  for(ge in names(genes)) {
    
    cat(paste("\nProcessing", ge, "\n"))
    
  # read file if 'genes' is character string with file path
  if(is.character(genes[[ge]]) & length(genes[[ge]])==1) {
    cat("\n\nReading gene list:", genes[[ge]], "\n")
    genes[[ge]] <- read.table(genes[[ge]], header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
  }
  
  if(is.data.frame(genes[[ge]])) {
     
     if (!is.null(newheader)) {# if 'newheader' is defined, it is used as names(DEgenes.unfilt)
        cat("\nNew header added to input file:", newheader, "\n")
        colnames(genes[[ge]]) <- newheader
      }
  
  columns2use <- c(id.column, sortcolumn, FCcolumn) 
  columnsfound <- columns2use[columns2use %in% names(genes[[ge]])]
    if (!all(columns2use %in% names(genes[[ge]]))) {
      warning("column(s) ", columns2use[!(columns2use %in% names(genes[[ge]]))], " not found in data header.")
      }  
  
  genes[[ge]] <- genes[[ge]][,columnsfound, drop=F]  # reduce dataframe to columns of interest
  }
      
  if(is.vector(genes[[ge]], mode="character")) {
    if(!is.null(names(genes[[ge]]))) { # if named vector with IDs as element names and quantitative information as elements
      genes[[ge]] <- data.frame(names(genes[[ge]]), genes[[ge]])
      names(genes[[ge]]) <- c(id.column, sortcolumn)        
        } else { # if unnamed vector with IDs but no quantitative information (no sortcolumn)
            id.column <- id.type
            genes[[ge]] <- data.frame(genes[[ge]], stringsAsFactors = F) # convert single vector to dataframe
            names(genes[[ge]]) <- c(id.column)        
            }
        }
  
  # now 'genes[[ge]]' is a data.frame with 1, 2 or 3 columns (id.column, sortcolumn, FCcolumn)
      genes[[ge]] <- genes[[ge]][!is.na(genes[[ge]][,id.column]), , drop=F] # remove rows without entry in id.column
      if (!is.null(sortcolumn)) {
        if (sortcolumn %in% names(genes[[ge]])) {    
          cat(paste("\nSorting input dataframe for", sortcolumn, ", decreasing =", sortdecreasing, "\n"))
          genes[[ge]] <- genes[[ge]][order(genes[[ge]][,sortcolumn], decreasing=sortdecreasing),] # sorting genes[[ge]] for quantitative variable
          
          # if multiple probes per gene, only the first gene entry is used 
          genes[[ge]] <- genes[[ge]][!duplicated(genes[[ge]][,id.column]),] 
          }
      }
     
 
      
  ## convert IDs to ENTREZ IDs if necessary
    if(id.type!="ENTREZID") {
    #cat(paste("\nConverting", id.type, "from input list to ENTREZIDs\n"))
    
      # entrezids <- bitr(genes[[ge]][,id.column], fromType=id.type, toType="ENTREZID", OrgDb=annotationdb)
      entrezids <- basicAnno(data=genes[[ge]], Symbol.column = id.column, Entrez.column = NULL, org=org)
    
      cat(paste(nrow(entrezids), "of",  length(unique(genes[[ge]][,id.column])), "unique", id.type, "mapped to", sum(!is.na(entrezids$ENTREZID)), "ENTREZIDs\n"))
      
      cat(paste("For dublicated", id.type, "only the first element is used (after sorting if applicable)\n")) 
      genes[[ge]] <- plyr::join(genes[[ge]], entrezids, by=id.column, type="right", match="first") 
      # joined dataframe contains all entries even without available ENTREZIDs and keeps order of 'genes[[ge]]'

    } else {names(genes[[ge]])[names(genes[[ge]])==id.column] <- "ENTREZID"}
    
     # 'genes[[ge]]' now definitively provides a column 'ENTREZID'   
     # remove entries without 'ENTREZID'
     genes[[ge]] <- genes[[ge]][!is.na(genes[[ge]][,"ENTREZID"]),]
     genes[[ge]] <- genes[[ge]][genes[[ge]][,"ENTREZID"] != "",]
     # if multiple probes per gene, only the first gene entry is used 
     # (after sorting, if sortcolumn available).
     genes[[ge]] <- genes[[ge]][!duplicated(genes[[ge]][,"ENTREZID"]),] 
      
      
  
  ### Filtering: Threshold applied to gene list: 'maxInputGenes' as well as 'sortcolumn.threshold' if applicable

     filtgenes <- filterGeneLists(genes[[ge]],
                        newheader=NULL,
                        filtercat1 = sortcolumn,
                        filtercat1.decreasing = sortdecreasing,
                        filtercat1.function = identity,
                        filtercat1.threshold= sortcolumn.threshold,
                        filtercat2 = FCcolumn,
                        filtercat2.decreasing = TRUE,
                        filtercat2.function = abs,
                        filtercat2.threshold= threshold_FC)
     
  cat(paste("\n", min(nrow(filtgenes), maxInputGenes), "unique genes for overrepresentation analysis.\n"))  
  filtgenes <- filtgenes[1:min(nrow(filtgenes), maxInputGenes), , drop=F]
  
  
  ## Processing list for Gene Set Enrichment Analysis (GSEA). Ordered named Vector needed.
  if (!is.null(sortcolumn)) { # quantitative data available?
    gseagenes <- genes[[ge]][,sortcolumn] 
    names(gseagenes) <- genes[[ge]][,"ENTREZID"]
    # for increasing sortcolumn (e.g. p-values). Low values transformed to high values
    if(sortdecreasing==F) {gseagenes <- fun.transf.incr.values(gseagenes)} 
    gseagenes <- sort(gseagenes, decreasing=T) # Vector in decreasing order.
    cat(paste("\n", length(gseagenes), "unique genes for gene set enrichment analysis.\n"))  
  }
  
  
  ###### read background list if available. Otherwise unfiltered ID list from 'genes[[ge]]' used as background Vector.
  if(!is.null(backgroundlist)) {
  cat("\nProcessing background list")
    
    if (all(backgroundlist=="genome")) {
      # all ENTREZIDs from the annotation db of the respective genome (given in org) used as background.
      backgroundlist <- keys(get(annotationdb), keytype="ENTREZID")
      cat(paste("\nUsing", length(backgroundlist), "ENTREZ IDs from", annotationdb, "as background.\n"))      
      } else {
    
          # read file if 'backgroundlist' is character string with file path
          if(is.character(backgroundlist) && length(backgroundlist)==1) {
            cat("\n\nReading background list:", backgroundlist, "\n")
            backgroundlist <- read.table(backgroundlist, header=is.null(newheaderBackground), sep="\t", na.strings = c("", " ", "NA")) # if no newheaderBackground, header must be in file
           }
 
          if(is.vector(backgroundlist, mode=c("character")) && !vector.processed) {
              backgroundlist <- as.data.frame(backgroundlist, stringsAsFactors =F) # convert vector to dataframe for processing in next chunk
              colnames(backgroundlist) <- id.column
              } 

          if(is.data.frame(backgroundlist) ) {
                  if(!is.null(newheaderBackground)) { # if 'newheaderBackground' is defined, it is used as header for backgroundlist
                    cat("\nNew header added to background list:", newheaderBackground, "\n")
                    colnames(backgroundlist) <- newheaderBackground
                    }
                  # convert IDs to ENTREZ IDs if necessary
                  if(id.type!="ENTREZID" && !("ENTREZID" %in% colnames(backgroundlist))) { 
                    # cat(paste("\nConverting", id.type, "from background list to ENTREZIDs\n"))
                    # backgroundlist <- bitr(backgroundlist, fromType=id.type, toType="ENTREZID", annoDb=annotationdb)
                    backgroundlist <- basicAnno(backgroundlist, Symbol.column = id.column, Entrez.column = NULL, org=org)
                    }
               entrez.column.bg <- ifelse(id.type=="ENTREZID", id.column, "ENTREZID")
               backgroundlist <- backgroundlist[!is.na(backgroundlist[,entrez.column.bg]),entrez.column.bg, drop=F] # remove rows without entry in id.column
               backgroundlist <- backgroundlist[backgroundlist[,entrez.column.bg] !="", entrez.column.bg, drop=F] # remove rows with entry "" in id.column
               backgroundlist <- as.character(unique(backgroundlist[,entrez.column.bg]))
               vector.processed <- TRUE # don't process this vector in 2nd run of loop
          } else {
              if(!vector.processed) {
                cat("\nBackground list used as given without any modification (assuming Entrez IDs):")
                print(head(backgroundlist))  
                cat("\n")
                }
              }

    } # end loading and processing of a supplied backgroundlist 
     
    cat(paste("\n", length(backgroundlist), "ENTREZ IDs in background list:\n"))
    print(head(backgroundlist))  
    cat("\n")
    
    } else { # if backgroundlist is NULL, all EntrezIDs from 'genes[[ge]]' used as background instead.

      backgroundlist <- as.character(unique(genes[[ge]][,"ENTREZID"]))
      cat(paste("\n", length(backgroundlist), "ENTREZ IDs used from input object as background:\n"))
      print(head(backgroundlist))  
      cat("\n")
    }
  


  

  #################### Start Enrichment analysis
  
  # remove categories not supported by clusterProfile
  enrichmentCat.found <- enrichmentCat %in% c("GO", "KEGG", "Reactome", "DO")
  enrichmentCat <- enrichmentCat[enrichmentCat.found]
  cat(paste("\nCategories used for enrichment analysis with clusterProfiler:", paste(enrichmentCat, collapse=" "), "\n"))
  
  # initialise result list
  enrichresult[[ge]] <- list()
  
  ##### Gene Ontology Enrichment
  ## Over representation analysis

  if ("GO" %in% enrichmentCat) {

  for (go in c("MF", "BP", "CC")) {
    if(nrow(filtgenes)>1) {


      enrichresult[[ge]][[paste0("Overrep_",go)]] <- clusterProfiler::enrichGO(gene = filtgenes[,"ENTREZID"],
                                    universe      = backgroundlist,
                                    OrgDb      = annotationdb,
                                    ont           = go,
                                    minGSSize    = minGSSize,
                                    maxGSSize    = maxGSSize,
                                    pAdjustMethod = pAdjustMethod, 
                                    pvalueCutoff  = enrich.p.valueCutoff,
                                    qvalueCutoff  = enrich.q.valueCutoff,
                                    readable      = TRUE) # EntrezIDs displayed as gene Symbols
      cat(paste("Overrepresentation analysis for", go, "GO-terms results in", nrow(as.data.frame(enrichresult[[ge]][[paste0("Overrep_",go)]])), "enriched terms.\n"))
    }

  # GO Gene Set Enrichment Analysis
    if (!is.null(sortcolumn)) { # quantitative data available?
      cat(paste("GeneSet Enrichment analysis for", go, "GO-terms results in "))
      enrichresult[[ge]][[paste0("GSEA_",go)]] <- clusterProfiler::gseGO(geneList = gseagenes,
                                                  OrgDb      = annotationdb,
                                                  ont          = go,
                                                  nPerm        = nPerm,
                                                  minGSSize    = minGSSize,
                                                  maxGSSize    = maxGSSize,
                                                  pvalueCutoff = enrich.p.valueCutoff,
                                                  pAdjustMethod = pAdjustMethod, 
                                                  verbose      = FALSE)
  cat(paste(nrow(as.data.frame(enrichresult[[ge]][[paste0("GSEA_",go)]])), "enriched terms.\n"))
  }
  
  } # end go-loop
} # end if GO


##### KEGG pathway enrichment
# KEGG over-representation test


  if ("KEGG" %in% enrichmentCat) {
    if(nrow(filtgenes)>1) {
      enrichresult[[ge]][[paste0("Overrep_","KEGG")]] <- clusterProfiler::enrichKEGG(gene = filtgenes[,"ENTREZID"],
                     organism      = org_alt_name,
                     universe    = backgroundlist, 
                     pvalueCutoff = enrich.p.valueCutoff, 
                     qvalueCutoff  = enrich.q.valueCutoff,
                     pAdjustMethod = pAdjustMethod, 
                     minGSSize    = minGSSize,
                     maxGSSize    = maxGSSize,
                     use_internal_data = FALSE)
    
    cat(paste("Overrepresentation analysis for KEGG pathways results in", nrow(as.data.frame(enrichresult[[ge]][[paste0("Overrep_","KEGG")]])), "enriched pathways\n"))
    }


# KEGG Gene Set Enrichment Analysis
if (!is.null(sortcolumn)) { # quantitative data available?
  cat(paste("GeneSet Enrichment analysis for KEGG pathways results in "))
  enrichresult[[ge]][[paste0("GSEA_","KEGG")]] <- clusterProfiler::gseKEGG(geneList  = gseagenes,
                 organism      = org_alt_name,
                 nPerm        = nPerm,
                 minGSSize    = minGSSize,
                 maxGSSize    = maxGSSize,
                 pvalueCutoff = enrich.p.valueCutoff,
                 pAdjustMethod = pAdjustMethod, 
                 verbose      = FALSE,
                 use_internal_data = FALSE)
  
  cat(paste(nrow(as.data.frame(enrichresult[[ge]][[paste0("GSEA_","KEGG")]])), "enriched pathways\n"))
  }
}


##### Reactome pathway enrichment
# Reactome over-representation test
  
if ("Reactome" %in% enrichmentCat) {
  if(nrow(filtgenes)>1) {
    enrichresult[[ge]][[paste0("Overrep_","Reactome")]] <- ReactomePA::enrichPathway(gene =  filtgenes[,"ENTREZID"],
                      organism      = org,
                      universe    = backgroundlist, 
                      pvalueCutoff = enrich.p.valueCutoff, 
                      qvalueCutoff  = enrich.q.valueCutoff,
                      pAdjustMethod = pAdjustMethod, 
                      minGSSize    = minGSSize,
                      maxGSSize    = maxGSSize,
                      readable     = TRUE)
  
  cat(paste("Overrepresentation analysis for Reactome pathways results in", nrow(as.data.frame(enrichresult[[ge]][[paste0("Overrep_","Reactome")]])), "enriched pathways\n"))
  }



# Reactome Gene Set Enrichment Analysis
  if (!is.null(sortcolumn)) { # quantitative data available?
    cat(paste("GeneSet Enrichment analysis for Reactome pathways results in "))
    enrichresult[[ge]][[paste0("GSEA_","Reactome")]] <- ReactomePA::gsePathway(geneList  = gseagenes,
                      organism      = org,
                      nPerm        = nPerm,
                      minGSSize    = minGSSize,
                      maxGSSize    = maxGSSize,
                      pvalueCutoff = enrich.p.valueCutoff,
                      pAdjustMethod = pAdjustMethod, 
                      verbose      = FALSE)
  
  cat(paste(nrow(as.data.frame(enrichresult[[ge]][[paste0("GSEA_","Reactome")]])), "enriched pathways\n"))
  }

}


  ##### Disease Ontology enrichment
  # DO over-representation test
  if ("DO" %in% enrichmentCat) {
    if(nrow(filtgenes)>1) {
      enrichresult[[ge]][[paste0("Overrep_","DO")]] <- DOSE::enrichDO(gene =  filtgenes[,"ENTREZID"],
                    ont="DO",
                    universe    = backgroundlist, 
                    pvalueCutoff = enrich.p.valueCutoff, 
                    qvalueCutoff  = enrich.q.valueCutoff,
                    pAdjustMethod = pAdjustMethod, 
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    readable     = TRUE)
    
    cat(paste("Overrepresentation analysis for Disease Ontology results in", nrow(as.data.frame(enrichresult[[ge]][[paste0("Overrep_","DO")]])), "enriched terms\n"))
    } 


  # DO Gene Set Enrichment Analysis
  if (!is.null(sortcolumn)) { # quantitative data available?
    cat(paste("GeneSet Enrichment analysis for Disease Ontology results in "))
    enrichresult[[ge]][[paste0("GSEA_","DO")]] <- DOSE::gseDO(geneList  = gseagenes,
                          nPerm        = nPerm,
                          minGSSize    = minGSSize,
                          maxGSSize    = maxGSSize,
                          pvalueCutoff = enrich.p.valueCutoff,
                          pAdjustMethod = pAdjustMethod, 
                          verbose      = FALSE)
  
  cat(paste(nrow(as.data.frame(enrichresult[[ge]][[paste0("GSEA_","DO")]])), "enriched terms.\n"))
  }


}



# Creating result tables and plot 
cat(paste("\nCreating result tables and plots in", file.path(projectfolder), "\n"))


for(usedcat in names(enrichresult[[ge]])) {

    if(nrow(as.data.frame(enrichresult[[ge]][[usedcat]])) >=1) {
  
    cat(paste("\nProcessing", usedcat))
    
    ## create output directory
    subfolder <- ifelse(ge == "genes", projectfolder, paste0(projectfolder, "/", ge))
    projectnamesuffix <- ifelse(ge == "genes", "", paste0(ge, "_"))
    if (!file.exists(file.path(subfolder))) {dir.create(file.path(subfolder), recursive=T)} 
      
    # Result table
    write.table(as.data.frame(enrichresult[[ge]][[usedcat]]), quote=F, row.names=F, sep="\t",
                file=file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_resulttable.txt")))
  
    # Enrichment map
    png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_enrichmentMap.png")), width = 300, height = 300, units = "mm", res=figure.res)
    try(DOSE::enrichMap(enrichresult[[ge]][[usedcat]], n = 30), silent=T)
    dev.off()
    
    # cnet plot (Remark: legend will be 'Fold change' even if other quantitative data are included)
    png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_cnetPlot.png")), width = 300, height = 300, units = "mm", res=figure.res)
    try(DOSE::cnetplot(enrichresult[[ge]][[usedcat]], showCategory = 5, categorySize="pvalue", 
             foldChange= if(!is.null(FCcolumn)) {filtgenes[,FCcolumn]} else {filtgenes[,sortcolumn]}), silent=T)
    dev.off()

    # enriched GO induced graph
    if(grepl("CC|MF|BP", usedcat)) {
      png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_GOgraph.png")), width = 300, height = 300, units = "mm", res=figure.res)
      try(clusterProfiler::plotGOgraph(enrichresult[[ge]][[usedcat]], firstSigNodes = 10, useFullNames = TRUE), silent=T)
      dev.off()
      }
    
    
    
    # Pathview maps for KEGG
    if(grepl("KEGG", usedcat)) {
     
      pathway.dir <- file.path(subfolder, paste0(usedcat, "_pathway_maps"))
      if (!file.exists(pathway.dir)) {dir.create(pathway.dir, recursive=T)} 
      workingdir <- getwd()
      setwd(pathway.dir)
      
      for(k in 1:nrow(enrichresult[[ge]][[usedcat]])) {
        
        if(grepl("gsea", usedcat, ignore.case =T)) { # column name with gene IDs differs for GSEA and overrep analysis
                temp_geneid_column <- "core_enrichment"} else {
                temp_geneid_column <- "geneID"
                }
        
        kegggenes <- enrichresult[[ge]][[usedcat]][k, temp_geneid_column]
        kegggenes <- unlist(strsplit(kegggenes, split="/", fixed=T))
        # 'genes[[ge]]' is a data.frame with 1, 2 or 3 columns (id.column, sortcolumn, FCcolumn)
        kegggenes <- genes[[ge]][genes[[ge]][,"ENTREZID"] %in% kegggenes, ]
        if(!is.null(FCcolumn)) {
          tmp.entrezids <- kegggenes$ENTREZID    
          kegggenes <- kegggenes[,FCcolumn]
          names(kegggenes) <- tmp.entrezids
          kegglimit <- max(abs(kegggenes))
          plot.col.key <- T
            } else {
              kegggenes <- kegggenes[,"ENTREZID"]
              kegglimit <- 1
              plot.col.key <- F
              }
            
        try(
         pathview::pathview(
            # numeric vector named with ENTREZ IDs
            gene.data  = kegggenes,
            pathway.id = enrichresult[[ge]][[usedcat]][k, "ID"],
            species    = org_alt_name,
            gene.idtype = "entrez",
            limit= list(gene= kegglimit, cpd=1),
            plot.col.key = plot.col.key)
         ,silent=F)
            }
      
      setwd(workingdir)
      }
    
  }
}

cat("\n")

} # end ge loop
  
# Detaching libraries not needed any more
detach_package(unique(pks2detach))
  

return(enrichresult)

} # end of function definition 




