
#' Wrapper for gene enrichment analysis using clusterProfiler
#' 
#' Applying over-representation analysis and gene set enrichment analysis to supplied gene list
#' 
#' This function uses genelist and optionally background gene list as input to perform enrichment analysis 
#' using the \code{clusterProfiler} and \code{DOSE} packages. By default, the function performs two kinds 
#' of analysis for all categories given in \code{enrichmentCat}.
#' \itemize{
#'  \item over-representation analysis: \code{DOSE} implements hypergeometric model to assess whether the 
#'  number of selected genes annotated with a term or pathway is larger than expected by chance.
#'  For this approach, usually a cut off has been applied to select genes of interest e.g. from differential
#'  expression analysis. All selected genes are treated with equal priority.
#'  \item gene set enrichment analysis (GSEA): for GSEA all genes under investigation can be used as input.
#'  They are ranked based on a corresponding quantitative value given in \code{sortcolumn}. Given a priori 
#'  defined set of genes S (e.g., genes sharing the same GO term), the goal of GSEA is to determine whether the 
#'  members of S are randomly distributed throughout the ranked gene list (L) or primarily found at the top or bottom.
#'  This approach can handle a situation where the extent of differential expression is small, but evidenced 
#'  in coordinated way in a set of related genes. GSEA aggregates the per gene statistics across genes within 
#'  a gene set, therefore making it possible to detect situations where all genes in a predefined set change 
#'  in a small but coordinated way. If no quantitative data is provided (\code{sortcolumn = NULL}), GSEA is skipped.
#' }
#' The names/IDs are converted to ENTREZ IDs (if necessary) prior to enrichment using the annotation package
#' for the species denoted in \code{org}. A background list with genes under investigation (e.g. expression array content) 
#' can be provided in \code{backgroundlist} to be used as background for over-representation analysis. 
#' Alternatively, all genes of the designated organism can be obtained from the respective
#' annotation package (For GSEA, all genes of the designated organism are used anyway). 
#' Optionally, quantitative data in \code{sortcolumn} can be used for sorting and filtering (using \code{sortcolumn.threshold}
#' or \code{maxInputGenes}) the input gene list for over-representation analysis, if not filtered prior to that.
#' The function generation visualisations for each enrichment result like cnetplots enrichment-maps and dotplots.
#' Enriched KEGG pathway maps (if any) are annotated for input genes using the \code{pathview} package.
#' If quantitative data in \code{sortcolumn} are not fold changes, fold changes may be given additionally in \code{FCcolumn} 
#' for annotation purposes in cnetplots as well as for KEGG pathway mapping, otherwise the data in \code{sortcolumn} 
#' is used for this. Be aware that the legend of the cnetplots will be "Fold Change" anyway.
#' 
#' 
#' @param genes the input data object can be given in several formats. 
#' \itemize{
#'  \item vector with gene names/IDs: all genes are used for over-representation analysis. The gene names/IDs
#'  must be in the format given in \code{id.type}.
#'  \item dataframe with gene names/IDs and quantitative values: genes can be filtered for the quantitative value 
#'  prior to over-representation analysis via the unfiltered list can be used as background list.
#'  Additionally to over-representation analysis, a GSEA is perfomed using the full gene list ranked for the 
#'  quantitative value. Columns of the dataframe must be named according to parameters \code{id.column} and \code{sortcolumn}.
#'  \item character with path to dataframe: as above. The dataframe is loaded from the given path.
#'  \item list of items: Each item is processed separately as indicated above.
#' }
#' @param newheader optional character vector with new header information for \code{genes} dataframe. Only relevant 
#'            if 'genes' is a dataframe (or character string with filepath to a table) with wrong or missing header. 
#'            NULL otherwise.
#' @param backgroundlist the background list of gene names/IDs for enrichment analysis can be given in several formats.
#' \itemize{
#'   \item vector with gene names/IDs: The gene names/IDs must be in the format given in \code{id.type}.
#'   \item dataframe with gene names/IDs: an \code{id.column} is needed as in \code{genes} with type of IDs given in \code{id.type}.
#'   \item "genome": all ENTREZ IDs from the annotation package of the respective organism (denoted in \code{org}) are used as background.
#'   \item NULL: full name/ID list from \code{genes} is used as background (i.e. that \code{genes} need to 
#'                 contain all genes under investigation without pre-filtering).      
#'   }              
#' @param newheaderBackground optional character vector with new header information for \code{backgroundlist}.
#' @param projectfolder character with directory for output files (will be generated if not existing).
#' @param projectname optional character prefix for output file names.
#' @param enrichmentCat character vector with categories to be enriched (\code{GO}: gene ontology (MF, BP, CC), 
#'                \code{KEGG}: KEGG pathways, \code{Reactome}: Reactome pathways, \code{DO}: Disease ontology).
#'                Disease ontology is for human only.
#' @param maxInputGenes (numeric) max number of top diff regulated elements used for over-representation analysis (or NULL).
#' @param id.type character with identifier type from annotation package (\code{"ENTREZID"} or \code{"SYMBOL"})
#'          Gene symbols will be converted to EntrezIDs prior to enrichment analysis.
#' @param id.column character with column name for identifier variable in \code{genes}. 
#' @param sortcolumn character with column name of quantitative data in \code{genes} used for ordering.
#'             If \code{Null}, ranking of genes is omitted and GSEA is not possible.
#' @param highValueHighPriority (logical) priority order of values in \code{sortcolumn}.
#'                 \code{TRUE}: high values have highest priority (e.g. fold changes).
#'                 \code{FALSE}: low values have highest priority (e.g. p-values); 
#'                 If \code{FALSE}, values in \code{sortcolumn} should be transformed prior to GSEA (see \code{fun.transform}).
#' @param sortcolumn.threshold numeric threshold for \code{sortcolumn} to be included in over-representation analysis.
#'                       \code{If highValueHighPriority=F, value < sortcolumn.threshold 
#'                       else value > sortcolumn.threshold}
#' @param fun.transform GSEA needs an input gene list with priority in decreasing order (high values have highest priority). 
#' Since quatitative values given in \code{sortcolumn} may have priority in increasing order (e.g. p-values), these values must 
#' be transformed by a custom function to generate priority in decreasing order prior to GSEA. A suitable function definition 
#' can be given in \code{fun.transform}, e.g. \code{function(x) {-log10(x)}} for p-values or \code{abs} for absolute
#' values of foldchange. 
#' @param FCcolumn (character) optional column name of foldchanges in \code{genes} if \code{sortcolumn} is used for a different data column.
#'           Used only for annotation cnetplot of enrichment results. Omitted if NULL
#' @param threshold_FC (numeric) Fold change threshold for filtering (threshold interpreted for log2 transformed foldchange values!)
#'               Only relevant for over-representation analysis if an unfiltered gene list is given in \code{genes}
#'               to allow for GSEA in parallel.
#' @param pAdjustMethod method for adjusting for multiple testing. 
#'                One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param org character with name of organism ("human", "mouse", "rat").
#' @param enrich.p.valueCutoff numeric p-value threshold for returned enrichment terms.
#' @param enrich.q.valueCutoff numeric q-value threshold for returned enrichment terms.
#' @param nPerm permutation numbers for gene set enrichment analysis
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param figure.res numeric resolution for output png.
#' 
#' 
#' @return List of enrichment-objects defined in \code{DOSE}-package (\code{enrichResult}-object for overrepresentation 
#'   analysis and \code{gseaResult}-objects for gene set enrichment analysis).                          
#'   Enrichment tables and plots are stored in the project folder as side effects.
#' 
#' @author Frank Ruehle
#' 
#' @export 
#' 
#' @note There are three key elements of the GSEA method (taken from 
#' https://bioconductor.org/packages/release/bioc/vignettes/DOSE/inst/doc/GSEA.html):
#' \itemize{
#'  \item Calculation of an Enrichment Score: The enrichment score (ES) represent the degree to which a set S is 
#'  over-represented at the top or bottom of the ranked list L. The score is calculated by walking down the list L, 
#'  increasing a running-sum statistic when we encounter a gene in S and decreasing when it is not. The magnitude 
#'  increment depends on the gene statistics (e.g., correlation of the gene with phenotype). The ES is the of the 
#'  maximum deviation from zero encountered in the random walk; it corresponds to a weighted Kolmogorov-Smirnov-like statistic
#'  \item Esimation of Significance Level of ES: The p-value of the ES is calculated using permutation test. Specifically, 
#'  we permute the gene labels of the gene list L and recompute the ES of the gene set for the permutated data, which 
#'  generate a null distribution for the ES. The p-value of the observed ES is then calculated relative to this null distribution.
#'  \item Adjustment for Multiple Hypothesis Testing: When the entire gene sets were evaluated, \code{DOSE} adjust the estimated 
#'  significance level to account for multiple hypothesis testing and also q-values were calculated for FDR control. 
#' }


    
## Usage 
wrapClusterProfiler <- function (genes, 
                         newheader = NULL, 
                         backgroundlist=NULL, 
                         newheaderBackground = NULL,
                         projectfolder= "clusterProfiler",
                         projectname="", 
                         enrichmentCat = c("GO", "KEGG", "Reactome", "DO"),
                         maxInputGenes = 100,  
                         id.type = "ENTREZID",
                         id.column = "ENTREZID",   
                         sortcolumn ="adj.P.Val",  
                         highValueHighPriority = FALSE, 
                         sortcolumn.threshold = 0.05,
                         fun.transform = function(x) {identity(x)},
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
  
  
  ## additional plotting parameter
  emapplot_showCategory = 20
  cnetplot_showCategory = 5
  dotplot_showCategory = 20
  plotGOgraph_firstSigNodes = 10
  
  
    
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
  pkg.bioc <- c("clusterProfiler", "DOSE", "enrichplot", "ReactomePA", "topGO", "pathview", annotationdb)
  pkg.cran <- c("plyr")

  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)
  
  # warnings
  if(!is.null(id.column)) {
    if(id.column!="ENTREZID" && "ENTREZID" %in% names(genes)) {
      cat("CAUTION: Column ENTREZIDs found. Consider these data for id.column.")
      }
  }
  
  if(!is.null(sortcolumn) & highValueHighPriority == FALSE) {
    cat(paste("\nCAUTION: You indicated high prority for low values in 'sortcolumn',", 
                "but GSEA will assess high values with  higher priority than", 
                "low values. Make sure that you have given a suitable function", 
                "in 'fun.transform' for transforming your data", 
                "(e.g. -log10(x) for p-values).\n"))
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
          cat(paste("\nSorting input dataframe for", sortcolumn, ", decreasing =", highValueHighPriority, "\n"))
          genes[[ge]] <- genes[[ge]][order(genes[[ge]][,sortcolumn], decreasing=highValueHighPriority),] # sorting genes[[ge]] for quantitative variable
          
          # if multiple probes per gene, only the first gene entry is used 
          genes[[ge]] <- genes[[ge]][!duplicated(genes[[ge]][,id.column]),] 
          }
      }
     
 
      
  ## convert IDs to ENTREZ IDs if necessary 
    if(id.type!="ENTREZID") {
      
      entrezids <- basicAnno(data=genes[[ge]], Symbol.column = id.column, Entrez.column = NULL, org=org)
    
      cat(paste(nrow(entrezids), id.type, "(",  length(unique(genes[[ge]][,id.column])), "unique) mapped to", sum(!is.na(entrezids$ENTREZID)), "ENTREZIDs\n"))
      
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
                        filtercat1.decreasing = highValueHighPriority,
                        filtercat1.function = identity,
                        filtercat1.threshold= sortcolumn.threshold,
                        filtercat2 = FCcolumn,
                        filtercat2.decreasing = TRUE,
                        filtercat2.function = abs,
                        filtercat2.threshold= threshold_FC)
     
  filtgenes[,"ENTREZID"] <- as.character(filtgenes[,"ENTREZID"]) # convert numeric EntrezIDs to characters
  cat(paste("\nTop", min(nrow(filtgenes), maxInputGenes), "unique genes selected for overrepresentation analysis:\n"))  
  print(head(filtgenes$ENTREZID))  
    cat("\n")

  filtgenes <- filtgenes[1:min(nrow(filtgenes), maxInputGenes), , drop=F]
  
  
  ## Processing list for Gene Set Enrichment Analysis (GSEA). Named Vector needed in decreasing order.
  if (!is.null(sortcolumn)) { # quantitative data available?

    # Requirements for input gene list according to https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList
    # - numeric vector: fold change or other type of numerical variable
    # - named vector: every number has a name, the corresponding gene ID
    # - sorted vector: number should be sorted in decreasing order
    
    gseagenes <- genes[[ge]][,sortcolumn] 
    names(gseagenes) <- genes[[ge]][,"ENTREZID"]
    # for increasing sortcolumn (e.g. p-values). Low values transformed to high values
    gseagenes <- fun.transform(gseagenes)
    gseagenes <- sort(gseagenes, decreasing=T) # Vector in decreasing order.
    cat(paste("\n", length(gseagenes), "unique genes for gene set enrichment analysis:\n"))
    print(head(names(gseagenes)))  
    cat("\n")
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
 
          if((is.vector(backgroundlist, mode=c("character")) | is.vector(backgroundlist, mode=c("numeric"))) && !vector.processed) {
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
  
  # remove categories not supported by clusterProfiler
  enrichmentCat.found <- enrichmentCat %in% c("GO", "KEGG", "Reactome", "DO")
  enrichmentCat <- enrichmentCat[enrichmentCat.found]
  if(("DO" %in% enrichmentCat) & org !="human") {
    cat("\nDiesease ontology (DO) available for human only!")
    enrichmentCat <- enrichmentCat[enrichmentCat!="DO"]
    if(length(enrichmentCat)==0) {stop("\nNo valid categories selected!")}  
  }
  
  cat(paste("\nCategories used for enrichment analysis with clusterProfiler:", paste(enrichmentCat, collapse=" "), "\n\n"))
  
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
                                    readable      = T) # if readable =T, EntrezIDs displayed as gene Symbols
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


  ##### Disease Ontology enrichment (for human only!)
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
      
    ## Result table
    write.table(as.data.frame(enrichresult[[ge]][[usedcat]]), quote=F, row.names=F, sep="\t",
                file=file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_resulttable.txt")))
  

    ## Enrichment map
    png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_emapPlot_top_", emapplot_showCategory, "_categories.png")), 
        width = 300, height = 300, units = "mm", res=figure.res)
    try(print(enrichplot::emapplot(enrichresult[[ge]][[usedcat]], showCategory = emapplot_showCategory, 
                             color = "p.adjust")) # "p.adjust" refers to enrichment p value of result object
              , silent=T) 
    dev.off()
    
 
    
    ## cnet plot (Remark: legend will be 'Fold change' even if other quantitative data are included)
    cnetplot_anno_column <- "ENTREZID" # by default, gene names are ENTREZIDs, because the enrichment is done with Entrez IDs
    if(grepl("Overrep", usedcat) & grepl("MF|BP|CC|Reactome", usedcat)) {
      # some enrichment functions (see line above) offer to translate enriched ENTREZIDs to SYMBOLs via parameter readable = TRUE
      # For those, SYMBOLs must be used to refer to quantitaive information in cnetplots (foldchange or other data).
      if(id.type == "SYMBOL") {
        cnetplot_anno_column <- id.column 
      } else {
        if(any(grepl("symbol", names(genes[[ge]]), ignore.case = T))) { # if dataframe containes unused symbol annotation
          cnetplot_anno_column <- grep("symbol", names(genes[[ge]]), value=T, ignore.case = T)[1]
        } else { # if SYMBOLs must be derived from ENTREZIDs
          genes[[ge]] <- basicAnno(data=genes[[ge]], Symbol.column = NULL, Entrez.column = "ENTREZID", org=org)
          cnetplot_anno_column <- "SYMBOL" 
        }
      }
    }
    # annotation of cnet plots (either by FCcolumn or sortcolumn)
    if(!is.null(FCcolumn)) {
      cnetplot_annotation <- genes[[ge]][,FCcolumn]
      names(cnetplot_annotation) <- genes[[ge]][,cnetplot_anno_column] #####
    } else {
      if(!is.null(sortcolumn)) {
        cnetplot_annotation <- genes[[ge]][,sortcolumn]
        names(cnetplot_annotation) <- genes[[ge]][,cnetplot_anno_column] ####
      } else {cnetplot_annotation <- NULL}
    }
    # start plot
    png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_cnetPlot_top_", cnetplot_showCategory, "_categories.png")), 
        width = 300, height = 300, units = "mm", res=figure.res)
    try(print(enrichplot::cnetplot(enrichresult[[ge]][[usedcat]], showCategory = cnetplot_showCategory, # categorySize="pvalue", 
             foldChange= cnetplot_annotation))
        , silent=T)
    dev.off()


    
    ## Dot plot of enriched terms
    png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_dotPlot_top_", dotplot_showCategory, "_categories.png")), 
        width = 300, height = 300, units = "mm", res=figure.res)
    try(print(enrichplot::dotplot(enrichresult[[ge]][[usedcat]], showCategory = dotplot_showCategory, color = "p.adjust", 
                                  title = usedcat))
        , silent=T)
    dev.off()
    
    
     ## enriched GO induced graph
    if(grepl("CC|MF|BP", usedcat)) {
      png(file.path(subfolder, paste0(projectname, projectnamesuffix, usedcat, "_GOgraph_first_", plotGOgraph_firstSigNodes, "_nodes.png")), 
          width = 300, height = 300, units = "mm", res=figure.res)
      try(clusterProfiler::plotGOgraph(enrichresult[[ge]][[usedcat]], firstSigNodes = plotGOgraph_firstSigNodes, useFullNames = TRUE), silent=T)
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
        
        kegggenes <- enrichresult[[ge]][[usedcat]][k, temp_geneid_column]  # Entrez IDs
        kegggenes <- unlist(strsplit(kegggenes, split="/", fixed=T))
        # 'genes[[ge]]' is a data.frame with 1, 2 or 3 columns (id.column, sortcolumn, FCcolumn)
        kegggenes <- genes[[ge]][genes[[ge]][,"ENTREZID"] %in% kegggenes, ]
        if(!is.null(FCcolumn)) {
          tmp.entrezids <- kegggenes$ENTREZID    
          kegggenes <- kegggenes[,FCcolumn]
          names(kegggenes) <- tmp.entrezids
          kegglimit <- c(min(kegggenes), max(kegggenes)) # kegglimit <- max(abs(kegggenes))
          plot.col.key <- T
        } else {
          
            if(!is.null(sortcolumn)) {
              tmp.entrezids <- kegggenes$ENTREZID    
              kegggenes <- kegggenes[,sortcolumn]
              names(kegggenes) <- tmp.entrezids
              kegglimit <- c(min(kegggenes), max(kegggenes))  
              plot.col.key <- T

              } else {
                kegggenes <- kegggenes[,"ENTREZID"]
                kegglimit <- 1
                plot.col.key <- F
                }
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
 ## remark: the warning produced by kegggraph or pathview do no harm. See https://support.bioconductor.org/p/96798/
    
       
  }
}

cat("\n")

} # end ge loop
  
# Detaching libraries not needed any more
detach_package(unique(pks2detach))
  

return(enrichresult)

} # end of function definition 




