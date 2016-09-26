
  ## Description
  # Analysis of differentially regulated elements of group comparisons
  
  ## Usage 
  comparefunc <- function( DEgenes.unfilt, 
                           newheader=NULL, 
                           comparisons=names(DEgenes.unfilt), 
                           
                           projectfolder= "GEX",
                           projectname= "", 
                           
                           pvalue = "P.Value", 
                           adj.pvalue = "adj.P.Val", 
                           p.value.threshold = 0.05, 
                           effectsize = "logFC",     
                           effect.threshold = log2(1.5), 
                           Symbol.column = "SYMBOL", 
                           Entrez.column = "ENTREZID",
                           
                           # Volcano plot parameter:
                             volcanoP = adj.pvalue, 
                           
                           # Enrichment analysis parameter:
                            use.clusterProfiler = TRUE,
                            use.davidAccount = FALSE, 
                            use.TFsearch = TRUE,   
                            maxInputGenes = 100,  
                            enrichmentCat = c("GO", "KEGG", "Reactome", "DO"), 
                            org = "human",
                           
                           # Further Subroutines 
                           use.pwOmics = FALSE   # not implemented yet
                           ) {
  
  
  
  ## Arguments
  #  DEgenes.unfilt: dataframe or character with path to dataframe or list of dataframes with unfiltered data
  #                  from diff expression / methylation analysis (e.g. results from groupsLimma-function). 
  #                  Required columns with p-value, adjusted p-value and effect size are defined below.
  #  newheader: optional character vector with new header information for 'DEgenes.unfilt' dataframe. Only relevant 
  #             if 'DEgenes.unfilt' is a dataframe or character string with filepath to a table with wrong or missing header.
  #             NULL otherwise.
  #  comparisons: character vector with respective names of group comparisons in format "groupA-groupB". 
  #               As default, 'DEgenes.unfilt' is expected to be list of dataframes with the respective
  #               group comparison as name of each list element. Otherwise give here the group comparisons 
  #               matching the data in 'DEgenes.unfilt'.  
  #  
  #  projectfolder: character with directory for output files (will be generated if not exisiting).
  #  projectname: optional character prefix for output file names.
  #  
  #  pvalue: character with column name of p values in 'DEgenes.unfilt'
  #  adj.pvalue: character with column name with p value adjusted for multiple testing in 'DEgenes.unfilt'
  #  p.value.threshold: numeric p-value threshold (is applied to p-values corrected for multiple testing in column 'adj.pvalue')
  #  effectsize: character with column name of effectsize in 'DEgenes.unfilt'  
  #  effect.threshold: numeric foldchange threshold
  #  Symbol.column: character with column name of Gene Symbols in 'DEgenes.unfilt'.
  #                 if NULL, SYMBOLs are derived from 'Entrez.column' using the annotation package for 'org'.
  #  Entrez.column: character with column name of ENTREZ IDs in 'DEgenes.unfilt'.
  #                 if NULL, ENTREZ IDs are derived from 'Symbol.column' using the annotation package for 'org'.
  #                 Either 'Symbol.column' or 'Entrez.column' must be specified.
  #  
  #  # Volcano plot parameter:
  #  volcanoP: specify if adjusted (adj.pvalue) or unadjusted (pvalue) p-values shall be used for volcano plots. 
  #  
  #  # Enrichment parameter:
  #  use.clusterProfiler: (boolean) use clusterProfiler for enrichment analysis. 
  #  use.davidAccount: email account for DAVID's Web Service connectivity or DAVIDWebService object.
  #                If NULL or FALSE, DAVID enrichment analysis is skipped.
  #  use.TFsearch: (boolean) use TFsearch.R for TF motive enrichment (human only) 
  #  maxInputGenes: (numeric) max number of top diff regulated elements used for overrepresentation analysis.
  #  enrichmentCat: character vector with DAVID categories to enrich for. E.g. c("GO", "KEGG", "Reactome", "DO")
  #                 for clusterProfiler or c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", 
  #                 "BIOCARTA", "UCSC_TFBS") for DAVID. Combinations possible if clusterProfiler and DAVID are used.
  #  org: character with species name ("human", "mouse", "rat").
  #  
  # # Further Subroutines 
  # use.pwOmics: (boolean) use pwOmics package for analysis (human only) (not implemented yet)
      
  
  ## Details
  # The function uses data.frames as input with differentially regulated elements 
  # (e.g. genes or methylation sites) resulting from a specified group comparison. Data can be delivered
  # as single dataframe, as list of dataframes named by group comparisons or as character containing
  # a filepath to be loaded. Either column names for Symbols or Entrez IDs must be specified. 
  # For each group comparison, volcano plots are generated and stored in the designated output folder.
  # Optional enrichment analsis (overrepresentation or gene set enrichment) is performed on ENTREZ IDs 
  # for specified categories and/or transcription factor binding sites using the clusterProf(), 
  # DavidEnrich() and TFsearch() subroutines.
  
  
  ## Value
  # No value returned. Output files are stored as side-effects.
  
  
  ## Author(s) 
  # Frank R?hle 
    
    
    
    # Load genomic annotation package (includes OrgDb, TxDb and GODb)
    annopkg <- switch(org, human = "Homo.sapiens", 
                      mouse = "Mus.musculus",
                      rat="Rattus.norvegicus")
    attach_package(pkg.bioc=annopkg)
    
    
  
    
    ### Create result directory if not yet existing 
    if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T) }
    if (!file.exists(file.path(projectfolder, "VolcanoPlots"))) {dir.create(file.path(projectfolder, "VolcanoPlots")) }
    
    
      

    # read file if 'DEgenes.unfilt' is character string with file path
    if(is.character(DEgenes.unfilt)) {
      cat("\n\nReading:", DEgenes.unfilt, "\n")
      DEgenes.unfilt <- read.table(DEgenes.unfilt, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
    }
    
    if(is.data.frame(DEgenes.unfilt) & !is.null(newheader)) {# if 'newheader' is defined, it is used as names(DEgenes.unfilt)
      cat("\nNew header added to input file:", newheader, "\n")
      names(DEgenes.unfilt) <- newheader
    }  
      
    # if DEgenes.unfilt is a single dataframe, it is coerced to a named list with 1 entry named after (the only) comparison
    if(is.data.frame(DEgenes.unfilt)) {DEgenes.unfilt <- list(comparisons=DEgenes.unfilt)} # if dataframe, coerced to named list
    names(DEgenes.unfilt) <- comparisons # only relevant if comparisons are defined separately
    
  ### now DEgenes.unfilt is a list with at least 1 dataframes
  
 
  
  
  
  
  #### Start processing gene list
  DEgenes <- list()
  
  if(length(DEgenes.unfilt)!=length(comparisons)) {stop("number of group comparisons does not fit to number of list elements!")}
  for (i in 1:length(comparisons)) {  # start loop for all applied group comparisons
    
      cat("\n\nProcessing group comparison", comparisons[i], "\n")
      
      ### convert SYMBOLs to ENTREZ IDs or vice versa if necessary (no additional rows introduced)
      DEgenes.unfilt[[comparisons[i]]] <- basicAnno(data=DEgenes.unfilt[[comparisons[i]]], 
                                                    Symbol.column = Symbol.column, Entrez.column = Entrez.column, org=org)
    
    # 'Entrez.column' name and 'Symbol.column' name need to be updated if NULL 
    if(is.null(Entrez.column) && any(grepl("entrez", names(DEgenes.unfilt[[comparisons[i]]]), ignore.case=T))) {
          Entrez.column <- grep("entrez", names(DEgenes.unfilt[[comparisons[i]]]), ignore.case=T, value=T)[1] }
    if(is.null(Symbol.column) && any(grepl("symbol", names(DEgenes.unfilt[[comparisons[i]]]), ignore.case=T))) {
            Symbol.column <- grep("symbol", names(DEgenes.unfilt[[comparisons[i]]]), ignore.case=T, value=T)[1] } 

      
      
      ### Filtering gene list
      cat("\nFiltering input list")  
      DEgenes[[comparisons[i]]] <- filterGeneLists(DEgenes.unfilt[[comparisons[i]]],
                                   newheader=NULL,
                                   filtercat1 = adj.pvalue,
                                   filtercat1.decreasing = FALSE,
                                   filtercat1.function = identity,
                                   filtercat1.threshold=p.value.threshold,
                                   filtercat2 = effectsize,
                                   filtercat2.decreasing = TRUE,
                                   filtercat2.function = abs,
                                   filtercat2.threshold= effect.threshold)
 
   
            
 if (nrow(DEgenes[[i]])>0) { # check if there are diff regulated elements
      
      #### Volcano Plots
      cat("\nWrite Volcano plot to", file.path(projectfolder, "VolcanoPlots", paste0("VolcanoPlot_", projectname, "_", comparisons[i], ".pdf")))
      pdf(file.path(projectfolder, "VolcanoPlots", paste0("VolcanoPlot_", projectname, "_", comparisons[i], ".pdf")), width = 14, height = 7) 
        plot(DEgenes.unfilt[[i]][,effectsize], -log10(DEgenes.unfilt[[i]][,volcanoP]), pch=16, cex=0.3, main=comparisons[i], xlab="logFC", ylab="-log10 p-value")
        abline(h=-log10(p.value.threshold),lty=c(1))
        abline(v=c(effect.threshold, -effect.threshold),lty=2) 

        # unadjusted pvalue < threshold in orange
          if(volcanoP==pvalue) {
            DEgenes.Punadj <- filterGeneLists(DEgenes.unfilt[[comparisons[i]]], # temporary dataframe filtered for unadjusted p-values
                                                         newheader=NULL,
                                                         filtercat1 = pvalue,
                                                         filtercat1.decreasing = FALSE,
                                                         filtercat1.function = identity,
                                                         filtercat1.threshold=p.value.threshold,
                                                         filtercat2 = effectsize,
                                                         filtercat2.decreasing = TRUE,
                                                         filtercat2.function = abs,
                                                         filtercat2.threshold= effect.threshold)
            
            
            points(DEgenes.Punadj[,effectsize], -log10(DEgenes.Punadj[,volcanoP]),pch=16, col="orange", cex=0.3) 
          }
  
        # adjusted pvalue < threshold in red
          points(DEgenes[[i]][,effectsize], -log10(DEgenes[[i]][,volcanoP]),pch=16, col="red", cex=0.3)  
          dev.off()
        
     } else {cat(paste("\n\nno differentially regulated elements for comparison:",comparisons[i],"\n"))} 
   } # end loop comparisons    
  

  
  
  ########### clusterProfiler Enrichment Analysis
  
  if (use.clusterProfiler==TRUE) {
    cat("\n\nEnrichment analysis with clusterProfiler\n")
    clusterProfResults <- list()
    for (i in 1:length(comparisons)) {  # start loop for all applied group comparisons

      clusterProfResults[[comparisons[i]]] <-  clusterprof(DEgenes.unfilt[[comparisons[i]]], 
                           newheader = NULL, 
                           backgroundlist=DEgenes.unfilt[[comparisons[i]]], 
                           newheaderBackground = NULL,
                           projectfolder=file.path(projectfolder, "clusterProfiler"),
                           projectname=comparisons[i], 
                           enrichmentCat = enrichmentCat,
                           maxInputGenes = maxInputGenes,  
                           id.type = "ENTREZID",
                           id.column = Entrez.column,   
                           sortcolumn = adj.pvalue,  
                           sortdecreasing = FALSE, 
                           sortcolumn.threshold = p.value.threshold ,
                           FCcolumn = effectsize,
                           threshold_FC= effect.threshold,
                           fun.transf.incr.vales = function(x) {-log10(x)},
                           org = org,
                           pAdjustMethod = "BH", 
                           enrich.p.valueCutoff = 0.05, 
                           enrich.q.valueCutoff = 0.05,
                           nPerm        = 1000,
                           minGSSize    = 10) 
    }
   }
  
  
           
      
   ## ######## DAVID enrichment
      if (!is.null(use.davidAccount)) {
        if (use.davidAccount != FALSE) {
        # Contact DAVIDWebService via try()-Function to avoid interruption if service is temporarily unavailable
        cat("\nTry to connect to DAVIDWebService\n")
        david <- NULL # initialise david object. Necessary if connection to DAVIDWebService failes.
        try2connect <- try(david <- if(class(use.davidAccount)=="DAVIDWebService") {use.davidAccount} else {DAVIDWebService$new(email=use.davidAccount)}, silent=F)
        
        if (class(david)=="DAVIDWebService") {
        
            DavidResults <- list()
            for (i in 1:length(comparisons)) {  # start loop for all applied group comparisons
              
              if (nrow(DEgenes[[comparisons[i]]])>0) { # check if there are diff regulated elements
                cat("\nProcessing group comparison", comparisons[i], "\n")    
             
                
                # id.type.david <- switch(id.type, ENTREZID="ENTREZ_GENE_ID", SYMBOL="OFFICIAL_GENE_SYMBOL", UNIGENE="UNIGENE")
                org.david <- switch(org, human="Homo sapiens", mouse="Mus musculus", rat="Rattus norvegicus")
                
                DavidResults[[comparisons[i]]] <- DavidEnrich(DEgenes[[comparisons[i]]], # vector with gene names or dataframe or character with path to dataframe
                                         newheader = NULL, # only relevant if DEgenes.unfilt is a character string indicating to a table with wrong or no header.
                                         backgroundlist=DEgenes.unfilt[[comparisons[i]]], 
                                         newheaderBackground = NULL,
                                         davidAccount = david, # email account for DAVID's Web Service connectivity or DAVIDWebService object
                                         projectfolder= file.path(projectfolder, "DAVID_Enrichment"),
                                         projectname=comparisons[i], # prefix for output file names
                                         maxInputGenes = maxInputGenes,  # max number of diff expressed probes used for enrichment analysis
                                         davidCat = enrichmentCat, # categories to enrich for
                                         id.type = "ENTREZ_GENE_ID",
                                         id.column = Entrez.column,   # column name for id.type
                                         sortcolumn = adj.pvalue,  # column to be sorted to select top entries
                                         sortdecreasing = FALSE, # drecreasing parameter for ordering
                                         sortcolumn.threshold = NULL, # provided data is already filtered
                                         org = org.david) 
                 
                } # end of if diff genes
              } # end loop comparisons 
      
        } else { # if no connection to DAVIDWebService
              cat("\nConnection to DAVID's Web Service failed or no feature IDs found. DAVID enrichment analysis is omitted.\n")
        }
  
        
  }} #### end of if use DAVID 
      
      
  
  
  
 
  
  ####### Transcription Factor Binding Sites
  
  #### PWMenrich  
 
  if(use.TFsearch) {
    for (i in 1:length(comparisons)) {  # start loop for all applied group comparisons
      if (nrow(DEgenes[[i]])>0) { # check if there are diff regulated elements
         
     cat("\n\nEnrichment Analysis for binding motivs of Transcription Factors\n")    
     cat("\nProcessing group comparison", comparisons[i], "\n")    
        
    TFreport <-  TFsearch(DEgenes[[i]], 
                          name.organism="hsapiens", 
                          projectfolder=file.path(projectfolder, "TFBS"),
                          projectname=paste(projectname, comparisons[i], sep="_"),
                          applyFilter = FALSE,
                          filtercat1 = pvalue,
                          filtercat1.decreasing = FALSE,
                          filtercat1.threshold=p.value.threshold, 
                          filtercat2 = effectsize,
                          filtercat2.decreasing = TRUE,
                          filtercat2.threshold= effect.threshold,
                          PromLookup = TRUE,
                          Entrez.col = Entrez.column, 
                          PromSeqUpstreamTSS = 2000,
                          PromSeqDownstreamTSS = 200)
                            
        } # end of if diff genes
    } # end loop comparisons 
  } # end TFsearch


  
  
  ######## pwOmics (needs diff genes and diff proteins!)
  ## Read in of protein data not implemented yet!
  if(use.pwOmics) {
    
    for (i in 1:length(comparisons)) {  # start loop for all applied group comparisons
      if (nrow(DEgenes[[i]])>0) { # check if there are diff regulated elements
        
        cat("\n\npwOmics Analysis\n")    
        cat("\nProcessing group comparison", comparisons[i], "\n")    
        
        inputpwOmics <- list(P=list(),
                             G=list(DEgenes.unfilt[[i]][!is.na(DEgenes.unfilt[[i]][,Symbol.column]),Symbol.column], 
                                    list(DEgenes[[i]][!is.na(DEgenes[[i]][,Symbol.column]),c(Symbol.column,effectsize)])))
        
        applypwOmics(inputpwOmics, 
                     projectdir=file.path(projectfolder, "pwOmics"),
                     projectname=paste(projectname, comparisons[i], sep="_"))
        
      } # end of if diff genes
    } # end loop comparisons 
  } # end pwOmics
  
  
  
 
  
  
} # end function definition
