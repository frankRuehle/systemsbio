
########### Filter gene lists

filterGeneLists <- function(genes,
                            newheader=NULL,
                            filtercat1 = "adj.P.Val",
                            filtercat1.decreasing = FALSE,
                            filtercat1.function = abs,
                            filtercat1.threshold = 0.05,
                            filtercat2 = "logFC",
                            filtercat2.decreasing = TRUE,
                            filtercat2.function = abs,
                            filtercat2.threshold = log2(1.5)) {
  
  # Filters dataframe or character directing to dataframe for up to two categories (e.g. p-value and foldchange).
  # Columns for filtercriteria must be included.
  # Values in filtercriteria may be transformed e.g. to ABSOLUTE values by filtercat.function = abs, which is needed for log foldchanges. 
  # No tranformation if filtercat.function == identity".
  
  
  #     genes: dataframe or character with path directing to table with genelist and columns to filter for.
  #     newheader: NULL if 'genes' already supplied with header. Character vector with new header otherwise.
  #     filtercat1: column name of first category to filter 'genes' (e.g. p-values). Skipped if NULL.
  #     filtercat1.decreasing (boolean): direction to order and filter filtercat1
  #     filtercat1.function: select transforming function for filter category1 (no quotes). e.g. abs for absolute values, identity for no transformation
  #     filtercat1.threshold: Threshold for filtercat1 or 'top123' for top Hits
  #     filtercat2: column name of second category to filter 'genes' (e.g. effect size). Skipped if NULL.
  #     filtercat2.decreasing (boolean): direction to order and filter filtercat2
  #     filtercat2.function: select transforming function for filter category2 (no quotes). E.g. abs for foldchanges
  #     filtercat2.threshold: Threshold for filtercat2 or 'top123' for top Hits
  
  # value: dataframe filtered for desired criteria and sorted for filtercat1 (if not Null).
  
  
  # read file if 'genes' is character string with file path
  if(is.character(genes)) {
    cat("\n\nReading:", genes, "\n")
    genes <- read.table(genes, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
  }
  if(!is.null(newheader)){# if 'newheader' is defined, it is used as names(genes)
    cat("\nNew header added to input file:", newheader, "\n")
    names(genes) <- newheader
  }  
  
  
  
  # apply filtercat2  
  if (!is.null(filtercat2)) {
    if(filtercat2 %in% names(genes)) {
      genes <- genes[order(filtercat2.function(genes[,filtercat2]), decreasing=filtercat2.decreasing),]
      
      if(is.character(filtercat2.threshold) & grepl("top", filtercat2.threshold, ignore.case=T)) { # select top hits
        toptreshold <- min(nrow(genes), as.numeric(sub("top", "", filtercat2.threshold, ignore.case=T)))
        cat("\nTop", toptreshold, "entries selected according to", filtercat2, "(decreasing =",filtercat2.decreasing, ")\n")
        genes <- genes[1:toptreshold,]
      } else {
        cat("\nFiltering Data for", filtercat2, ", decreasing =",filtercat2.decreasing, ", treshold =", filtercat2.threshold, "\n")
        if(filtercat2.decreasing==FALSE) { # less than threshold
          genes <- genes[!is.na(genes[,filtercat2]) & filtercat2.function(genes[,filtercat2]) < filtercat2.threshold, ]
        } else { # greater than threshold
          genes <- genes[!is.na(genes[,filtercat2]) & filtercat2.function(genes[,filtercat2]) > filtercat2.threshold, ]
        }
      } 
    }}
  
  
  # apply filtercat1
  if (!is.null(filtercat1)) {
    if(filtercat1 %in% names(genes)) {
      genes <- genes[order(filtercat1.function(genes[,filtercat1]), decreasing=filtercat1.decreasing),]
      
      if(is.character(filtercat1.threshold) & grepl("top", filtercat1.threshold, ignore.case=T)) { # select top hits
        toptreshold <- min(nrow(genes), as.numeric(sub("top", "", filtercat1.threshold, ignore.case=T)))
        cat("\nTop", toptreshold, "entries selected according to", filtercat1, "(decreasing =",filtercat1.decreasing, ")\n")
        genes <- genes[1:toptreshold,]
      } else {
        cat("Filtering Data for", filtercat1, ", decreasing =",filtercat1.decreasing, ", treshold =", filtercat1.threshold, "\n")
        if(filtercat1.decreasing==FALSE) { # less than threshold
          genes <- genes[!is.na(genes[,filtercat1]) & filtercat1.function(genes[,filtercat1]) < filtercat1.threshold, ] 
          
        } else { # greater than threshold
          genes <- genes[!is.na(genes[,filtercat1]) & filtercat1.function(genes[,filtercat1]) > filtercat1.threshold, ]
        }
      } 
    }}
  
  
  return(genes)
  
} # end function definition 'filterGeneLists'


