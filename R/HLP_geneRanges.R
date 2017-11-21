#' Help functions for handling gene ranges
#' 
#' Keep all meta data when overlapping gene ranges or converting to dataframe.
#' 
#' By default, meta data is lost when GRanges are merged or overlapped.
#' \code{subsetByOverlaps.keepAllMeta} returns overlap of two GRanges objects with 
#' merged meta data from both. Meta data stored in \code{CompressedCharacterList}
#' is collapsed to a single column. 
#' \code{granges2df} generates a data.frame from an GRanges containing the meta data.
#' Meta data stored in \code{CompressedCharacterList} is collapsed to a single column.
#' 
#' @param gr1,gr2 GRanges object
#' @param write.ranges.tofile character with file path. If given, a data.frame is generated
#' from the returned GRanges object and written to this destination. Omitted if NULL.
#' @param addStart numeric with mumber of bases to be added to the start coordinate
#' when converting a GRanges object to a dataframe. E.g. necessary for switching between 
#' 0-based and 1-based genome coordinate systems.
#' 
#' @return \code{subsetByOverlaps.keepAllMeta} returnes GRanges object containing 
#' overlap and meta data from input ranges. 
#' \code{granges2df} returns a data.frame with genomic coordinates and meta data 
#' from input object.
#' 
#' @author Frank Ruehle
#' 
#' @export subsetByOverlaps.keepAllMeta 
#' @export granges2df  




#### subsetByOverlaps function which keeps meta data from both objects
subsetByOverlaps.keepAllMeta <- function(gr1, gr2, write.ranges.tofile = NULL, addStart=0) {
  
  ranges <- subsetByOverlaps(gr1, gr2) # query, subject
  
  hits <- findOverlaps(ranges, gr2) 
  idx <- unique(subjectHits(hits)) 
  
  names <- CharacterList(split(names(gr2)[subjectHits(hits)], queryHits(hits))) # row names of gr2 object added to meta data
  names <- sapply(names, function(x) {paste(unique(x), collapse="; ") })
  mcols(ranges) <- DataFrame(mcols(ranges), names)
  
  for(m in names(mcols(gr2))) { # meta columns of gr2 summarized and added to meta data of gr1
    meta <- mcols(gr2)[subjectHits(hits),m]
    if (is.factor(meta)) {meta <- as.character(meta)}
    meta <- CharacterList(split(meta, queryHits(hits)))
    mcols(ranges) <- DataFrame(mcols(ranges), metaname=meta)
    names(mcols(ranges))[names(mcols(ranges))=="metaname"] <- m
    if(class(mcols(ranges)[,m])=="CompressedCharacterList") {
      mcols(ranges)[,m] <-  sapply(mcols(ranges)[,m], function(x) {paste(unique(x), collapse="; ") })
    } else {
      mcols(ranges)[,m] <-  mcols(ranges)[,m]
    }
  }
  
  if(!is.null(write.ranges.tofile)) {
    df <- granges2df(ranges, addStart=addStart) 
    write.table(df, write.ranges.tofile, sep="\t", quote = F, row.names = F)
  }
  
  return(ranges)
} 




#### Convert GRanges object to dataframe
granges2df <- function(gr1, addStart=0) {
  
  if(is.null(names(gr1))) {names(gr1) <- 1:length(gr1)}
  df <- data.frame(names=names(gr1),
                   seqnames=seqnames(gr1),
                   start=start(gr1) + addStart,  # -1: BED uses 0-based coordinates
                   end=end(gr1),
                   strand=strand(gr1))
  
  if(addStart!=0) {cat(addStart, "bp added to start coordinate.\n")}
  
  for(m in names(mcols(gr1))) {
    if(class(mcols(gr1)[,m])=="CompressedCharacterList") {
      df[,m] <-  sapply(mcols(gr1)[,m], function(x) {paste(unique(x), collapse="; ") })
    } else {
      df[,m] <-  mcols(gr1)[,m]
    }
  }
  
  return(df)   
}
