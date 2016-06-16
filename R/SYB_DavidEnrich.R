

  ## Description
  # Wrapper for DAVID enrichment analysis using the RDAVIDWebService
  
  ## Usage 
  DavidEnrich <- function (genes, 
                           newheader = NULL, 
                           backgroundlist=NULL, 
                           newheaderBackground = NULL,
                           davidAccount = NULL, 
                           projectfolder=pipepar[["outdir"]],
                           projectname="", 
                           maxDavid = 100,  
                           davidCat = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", "UCSC_TFBS"), 
                           id.type = "ENTREZ_GENE_ID",
                           id.column = "ENTREZID",   
                           sortcolumn ="adj.P.Val",  
                           sortdecreasing = FALSE, 
                           sortcolumn.threshold = 0.05,
                           org = "Homo sapiens",
                           threshold4enrich = 0.05, 
                           mincount=2L 
                           ) 
                           {

  
  ## Arguments
  # genes: vector with gene names to analyse or dataframe or character with path to dataframe. Supplied
  #        dataframe needs to include a column with identifiers specified in 'id.column'.
  # newheader: optional character vector with new header information for 'genes' dataframe. Only relevant 
  #            if 'genes' is a dataframe or character string with filepath to a table with wrong or missing header. 
  #            NULL otherwise.
  # backgroundlist: optional dataframe containing background list for enrichment analysis. 
  #                 If NULL, species-specific background from DAVID is used instead.
  # newheaderBackground: optional character vector with new header information for 'backgroundlist' dataframe.
  # davidAccount: email account for DAVID's Web Service connectivity or DAVIDWebService object.
  #               Registration for RDAVIDWebService is required at "https://david.ncifcrf.gov/webservice/register.htm".
  # projectfolder: character with directory for output files (will be generated if not exisiting).
  # projectname: optional character prefix for output file names.
  # maxDavid: (numeric) max number of top diff regulated elements used for enrichment analysis.
  # davidCat: character vector with DAVID categories to enrich for. 
  # id.type: character with identifier DAVID uses to enrich for (e.g. "ENTREZ_GENE_ID" or "OFFICIAL_GENE_SYMBOL")
  # id.column: character with respective column name for id.type in DEgenes.unfilt'.
  # sortcolumn: character with column name to be sorted for selecting top entries.
  # sortdecreasing: (boolean) drecreasing parameter for order of 'sortcolumn'.
  # sortcolumn.threshold: numeric threshold for 'sortcolumn' to be included in overepresentation analysis.
  #                       If sortdecreasing=F, value < sortcolumn.threshold 
  #                       else value > sortcolumn.threshold
  # org: character with species name DAVID uses for enrichment (e.g."Homo sapiens" or "Mus musculus").
  # threshold4enrich: numeric p-value threshold for enrichment terms.
  # mincount: numeric minimum gene count in enrichment terms.  
    
    
    
  ## Details
  # Function uses genelist and optionally background gene list as input to perform enrichment analysis 
  # using the RDAVIDWebService. Gene list can be supplied as vector with ids given in id.type, dataframe
  # with or character string with path of data file. If a vector is given, sorting parameter will be ignored.
  # Enrichment categories are defined in 'davidCat' and species is defined in 'org'.
  # If no connection to RDAVIDWebService can be established, NULL is returned.
  
    
  ## value
  # Functional annotation chart as R object. 
  # Additionally, a functional annotation chart is stored in the project folder as side effect.
   
 
  ## Author(s) 
  # Frank Rühle 
    
    
  
    
  # create output directory
  if (!file.exists(file.path(projectfolder, "DAVID_Enrichment"))) {dir.create(file.path(projectfolder, "DAVID_Enrichment"), recursive=T)} 
  
  
  
  
  
  ## install/load required packages from CRAN and Bioconductor
  source("http://bioconductor.org/biocLite.R")
  pcksBioc <- c("RDAVIDWebService", "GSEABase", "KEGG.db", "KEGGgraph")
  
  for (lib in pcksCRAN) { # CRAN packages
    if (!is.element(pcksCRAN[p], installed.packages()[,1])) {
      install.packages(pcksCRAN[p], dependencies=TRUE)}
    require(pcksCRAN[p], character.only = TRUE) 
  }
 
  
  
  
  
  
  # read file if 'genes' is character string with file path
  if(is.character(genes) & length(genes)==1) {
    cat("\n\nReading gene list:", genes, "\n")
    genes <- read.table(genes, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
  }
  if(is.data.frame(genes) & !is.null(newheader)) {# if 'newheader' is defined, it is used as names(DEgenes.unfilt)
    cat("\nNew header added to input file:", newheader, "\n")
    colnames(genes) <- newheader
  }
  if(is.vector(genes)) {genes <- unique(genes)[1:min(maxDavid, length(unique(genes)))]
    } else { # if genes is a data.frame
      genes <- genes[!is.na(genes[,id.column]),]
      genes <- genes[order(genes[,sortcolumn], decreasing=sortdecreasing),]
      
      if(!is.null(sortcolumn.threshold)) {
        if(sortdecreasing==F) { # for increasing sortcolumn (e.g. p-values) value < threshold
          genes <- genes[genes[,sortcolumn] < sortcolumn.threshold, ]
          cat(paste("\nGenelist filtererd for", sortcolumn, "<", sortcolumn.threshold, "(", nrow(filtgenes), "genes remaining)\n"))
          
        } else { # for decreasing sortcolumn (e.g. fold changes) value > threshold
          genes <- genes[genes[,sortcolumn] > sortcolumn.threshold, ]
          cat(paste("\nGenelist filtererd for", sortcolumn, ">", sortcolumn.threshold, "(", nrow(filtgenes), "genes remaining)\n"))
        }}
        
      genes <- unique(genes[,id.column])[1:min(maxDavid, length(unique(genes[,id.column])))]
      }
      
  
  
  # read background list if available. Otherwise species-specific background from DAVID
  if(!is.null(backgroundlist)) {
    # read file if 'Dgenes' is character string with file path
    if(is.character(backgroundlist) & length(backgroundlist)==1) {
      cat("\n\nReading background list:", backgroundlist, "\n")
      backgroundlist <- read.table(backgroundlist, header=is.null(newheaderBackground), sep="\t", na.strings = c("", " ", "NA")) # if no newheaderBackground, header must be in file
    }
    if(is.data.frame(backgroundlist) & !is.null(newheaderBackground)) {# if 'newheaderBackground' is defined, it is used as names(DEgenes.unfilt)
      cat("\nNew header added to input file:", newheaderBackground, "\n")
      colnames(backgroundlist) <- newheaderBackground
      }
    if(is.vector(backgroundlist)) {backgroundlist <- unique(backgroundlist)
    } else { # if genes is a data.frame
      backgroundlist <- unique(backgroundlist[!is.na(backgroundlist[,id.type]),id.type])
    }
  }
  
 
  ### remove categories not supported by DAVID
  davidCat.found <- davidCat %in% c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", "UCSC_TFBS")
  davidCat <- davidCat[davidCat.found]
  cat(paste("\nCategories used for enrichment analysis with DAVID:", paste(davidCat, collapse=" "), "\n"))
  
  
  
  
  david <- NULL # initialise david object. Necessary if connection to DAVIDWebService failes.
  if(class(davidAccount)=="DAVIDWebService") {
    david <- davidAccount} else {
    # Contact DAVIDWebService via try()-Function to avoid interruption if service is temporarily unavailable
    cat("\nTry to connect to DAVIDWebService\n")
    try2connect <- try(david <- DAVIDWebService$new(email=davidAccount), silent=F)
    }
  
  if (class(david)=="DAVIDWebService") {
        
        cat("\n\nDAVID enrichment for categories:", davidCat, "\n")    
        setAnnotationCategories(david, davidCat)
        setCurrentSpecies(david, species=org)  # or species=1 (if human is at first position)?
        ForegroundGenes <- addList(david, genes, idType=id.type,  listName="genelist", listType="Gene")
        if(!is.null(backgroundlist)) {
          BackgroundGenes <- addList(david, backgroundlist, idType=id.type, listName="background", listType="Background")
          }
        getFunctionalAnnotationChartFile(david, fileName=file.path(projectfolder, "DAVID_Enrichment", paste0("Enrich_", projectname, ".txt")), threshold=threshold4enrich, count=mincount)
        davidResult <- getFunctionalAnnotationChart(david, threshold=threshold4enrich, count=mincount)
    
    } else { # if no connection to DAVIDWebService
    cat("\nConnection to DAVID's Web Service failed or no IDs found. DAVID enrichment analysis is omitted.\n")
    return()
    }

  
  # Detaching libraries not needed any more
  for (lib in pcksCRAN) {detach_package(lib)}
  
  return(davidResult) # Get functional annotation chart as R object.
} # end of function definition 






  ## Description
  # Wrapper for DAVID enrichment analysis using the RDAVIDWebService
  
  ## Usage 
  DavidEnrich <- function (genes, 
                           newheader = NULL, 
                           backgroundlist=NULL, 
                           newheaderBackground = NULL,
                           davidAccount = NULL, 
                           projectfolder=pipepar[["outdir"]],
                           projectname="", 
                           maxDavid = 100,  
                           davidCat = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", "UCSC_TFBS"), 
                           id.type = "ENTREZ_GENE_ID",
                           id.column = "ENTREZID",   
                           sortcolumn ="adj.P.Val",  
                           sortdecreasing = FALSE, 
                           sortcolumn.threshold = 0.05,
                           org = "Homo sapiens",
                           threshold4enrich = 0.05, 
                           mincount=2L 
                           ) 
                           {

  
  ## Arguments
  # genes: vector with gene names to analyse or dataframe or character with path to dataframe. Supplied
  #        dataframe needs to include a column with identifiers specified in 'id.column'.
  # newheader: optional character vector with new header information for 'genes' dataframe. Only relevant 
  #            if 'genes' is a dataframe or character string with filepath to a table with wrong or missing header. 
  #            NULL otherwise.
  # backgroundlist: optional dataframe containing background list for enrichment analysis. 
  #                 If NULL, species-specific background from DAVID is used instead.
  # newheaderBackground: optional character vector with new header information for 'backgroundlist' dataframe.
  # davidAccount: email account for DAVID's Web Service connectivity or DAVIDWebService object.
  #               Registration for RDAVIDWebService is required at "https://david.ncifcrf.gov/webservice/register.htm".
  # projectfolder: character with directory for output files (will be generated if not exisiting).
  # projectname: optional character prefix for output file names.
  # maxDavid: (numeric) max number of top diff regulated elements used for enrichment analysis.
  # davidCat: character vector with DAVID categories to enrich for. 
  # id.type: character with identifier DAVID uses to enrich for (e.g. "ENTREZ_GENE_ID" or "OFFICIAL_GENE_SYMBOL")
  # id.column: character with respective column name for id.type in DEgenes.unfilt'.
  # sortcolumn: character with column name to be sorted for selecting top entries.
  # sortdecreasing: (boolean) drecreasing parameter for order of 'sortcolumn'.
  # sortcolumn.threshold: numeric threshold for 'sortcolumn' to be included in overepresentation analysis.
  #                       If sortdecreasing=F, value < sortcolumn.threshold 
  #                       else value > sortcolumn.threshold
  # org: character with species name DAVID uses for enrichment (e.g."Homo sapiens" or "Mus musculus").
  # threshold4enrich: numeric p-value threshold for enrichment terms.
  # mincount: numeric minimum gene count in enrichment terms.  
    
    
    
  ## Details
  # Function uses genelist and optionally background gene list as input to perform enrichment analysis 
  # using the RDAVIDWebService. Gene list can be supplied as vector with ids given in id.type, dataframe
  # with or character string with path of data file. If a vector is given sorting parameter will be ignored.
  # Enrichment categories are defined in 'davidCat' and species is defined in 'org'.
  # If no connection to RDAVIDWebService can be established, NULL is returned
  
    
  ## value
  # Functional annotation chart as R object. 
  # Additionally, a functional annotation chart is stored in the project folder as side effect.
   
 
  ## Author(s) 
  # Frank Rühle 
    
    
  
    
  # create output directory
  if (!file.exists(file.path(projectfolder, "DAVID_Enrichment"))) {dir.create(file.path(projectfolder, "DAVID_Enrichment"), recursive=T)} 
  
  
  
  
  
  ## install/load required packages from CRAN and Bioconductor
  source("http://bioconductor.org/biocLite.R")
  pcksBioc <- c("RDAVIDWebService", "GSEABase", "KEGG.db", "KEGGgraph")
  
  for (lib in pcksCRAN) { # CRAN packages
    if (!is.element(pcksCRAN[p], installed.packages()[,1])) {
      install.packages(pcksCRAN[p], dependencies=TRUE)}
    require(pcksCRAN[p], character.only = TRUE) 
  }
 
  
  
  
  
  
  # read file if 'genes' is character string with file path
  if(is.character(genes) & length(genes)==1) {
    cat("\n\nReading gene list:", genes, "\n")
    genes <- read.table(genes, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
  }
  if(is.data.frame(genes) & !is.null(newheader)) {# if 'newheader' is defined, it is used as names(DEgenes.unfilt)
    cat("\nNew header added to input file:", newheader, "\n")
    colnames(genes) <- newheader
  }
  if(is.vector(genes)) {genes <- unique(genes)[1:min(maxDavid, length(unique(genes)))]
    } else { # if genes is a data.frame
      genes <- genes[!is.na(genes[,id.column]),]
      genes <- genes[order(genes[,sortcolumn], decreasing=sortdecreasing),]
      
      if(!is.null(sortcolumn.threshold)) {
        if(sortdecreasing==F) { # for increasing sortcolumn (e.g. p-values) value < threshold
          genes <- genes[genes[,sortcolumn] < sortcolumn.threshold, ]
          cat(paste("\nGenelist filtererd for", sortcolumn, "<", sortcolumn.threshold, "(", nrow(filtgenes), "genes remaining)\n"))
          
        } else { # for decreasing sortcolumn (e.g. fold changes) value > threshold
          genes <- genes[genes[,sortcolumn] > sortcolumn.threshold, ]
          cat(paste("\nGenelist filtererd for", sortcolumn, ">", sortcolumn.threshold, "(", nrow(filtgenes), "genes remaining)\n"))
        }}
        
      genes <- unique(genes[,id.column])[1:min(maxDavid, length(unique(genes[,id.column])))]
      }
      
  
  
  # read background list if available. Otherwise species-specific background from DAVID
  if(!is.null(backgroundlist)) {
    # read file if 'Dgenes' is character string with file path
    if(is.character(backgroundlist) & length(backgroundlist)==1) {
      cat("\n\nReading background list:", backgroundlist, "\n")
      backgroundlist <- read.table(backgroundlist, header=is.null(newheaderBackground), sep="\t", na.strings = c("", " ", "NA")) # if no newheaderBackground, header must be in file
    }
    if(is.data.frame(backgroundlist) & !is.null(newheaderBackground)) {# if 'newheaderBackground' is defined, it is used as names(DEgenes.unfilt)
      cat("\nNew header added to input file:", newheaderBackground, "\n")
      colnames(backgroundlist) <- newheaderBackground
      }
    if(is.vector(backgroundlist)) {backgroundlist <- unique(backgroundlist)
    } else { # if genes is a data.frame
      backgroundlist <- unique(backgroundlist[!is.na(backgroundlist[,id.type]),id.type])
    }
  }
  
 
  ### remove categories not supported by DAVID
  davidCat.found <- davidCat %in% c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", "UCSC_TFBS")
  davidCat <- davidCat[davidCat.found]
  cat(paste("\nCategories used for enrichment analysis with DAVID:", paste(davidCat, collapse=" "), "\n"))
  
  
  
  
  david <- NULL # initialise david object. Necessary if connection to DAVIDWebService failes.
  if(class(davidAccount)=="DAVIDWebService") {
    david <- davidAccount} else {
    # Contact DAVIDWebService via try()-Function to avoid interruption if service is temporarily unavailable
    cat("\nTry to connect to DAVIDWebService\n")
    try2connect <- try(david <- DAVIDWebService$new(email=davidAccount), silent=F)
    }
  
  if (class(david)=="DAVIDWebService") {
        
        cat("\n\nDAVID enrichment for categories:", davidCat, "\n")    
        setAnnotationCategories(david, davidCat)
        setCurrentSpecies(david, species=org)  # or species=1 (if human is at first position)?
        ForegroundGenes <- addList(david, genes, idType=id.type,  listName="genelist", listType="Gene")
        if(!is.null(backgroundlist)) {
          BackgroundGenes <- addList(david, backgroundlist, idType=id.type, listName="background", listType="Background")
          }
        getFunctionalAnnotationChartFile(david, fileName=file.path(projectfolder, "DAVID_Enrichment", paste0("Enrich_", projectname, ".txt")), threshold=threshold4enrich, count=mincount)
        davidResult <- getFunctionalAnnotationChart(david, threshold=threshold4enrich, count=mincount)
    
    } else { # if no connection to DAVIDWebService
    cat("\nConnection to DAVID's Web Service failed or no IDs found. DAVID enrichment analysis is omitted.\n")
    return()
    }

  
  # Detaching libraries not needed any more
  for (lib in pcksCRAN) {detach_package(lib)}
  
  return(davidResult) # Get functional annotation chart as R object.
} # end of function definition 




