#' Process LNCipedia bed-file
#' 
#' Load and/or modify lncRNA data from LNCipedia.
#' 
#' LNCipedia date is loaded from bed file or provided as dataframe or GRanges object.
#' For reducing complexity in plotting purposes, all transcripts belonging to a lncRNA gene are given 
#' with respective gene coordinates. For each transcript, exon ccordinates are calculated from
#' "blockStarts" and "blockSizes" data from LNCipedia.
#' 
#' @param LNCipedia Either character with filepath to LNCipedia bed-file (no header), dataframe 
#' or GRanges object made from LNcipedia bed-file.
#' @param collapseTranscripts2Genes logical. If TRUE, all transcripts of a gene are assigned with the coordinates
#' of the corresponding gene. Recommended for plotting to reduce complexity.
#' @param makeExonRanges logical. If True, an GRanges object for exons is generated using the 
#' "blockStarts" and "blockSizes" data from LNCipedia. Transcript or gene coordinates are stored in metadata.
#' @param addBases numeric. Add bases to both sides of the ranges.
#' 
#' @return GRanges object containing LNCipedia transcripts or exons
#' 
#' @author Frank Ruehle
#' 
#' @export processLNCipedia 



processLNCipedia <- function(LNCipedia, collapseTranscripts2Genes=T, makeExonRanges=T, addBases=0) {
  
  
  ## read LNCipedia
  if (is.character(LNCipedia)) { # read Lncipedia bed-file
    LNCipedia <- read.table(LNCipedia, header=F, stringsAsFactors = F)
    colnames(LNCipedia) <- c("chrom", "start", "end", "name", "score", "strand", "thickStart", 
                             "thickEnd", "itemRgb", "blockCount",  "blockSizes", "blockStarts")
      }
  
  if (is.data.frame(LNCipedia)) {
      # check if all required columns present in data.frame
      requiredColumns <- c("chrom", "start", "end", "name", "strand", "blockCount",  "blockSizes", "blockStarts")
      if(!all(requiredColumns %in% colnames(LNCipedia))) {
        stop("Columns missing:", paste(requiredColumns[!requiredColumns %in% colnames(LNCipedia)], collapse=" "), "\n")}
      
      LNCipedia$chrom <- gsub("chr", "", LNCipedia$chrom, ignore.case=T)
      # mage GRanges object from dataframe
      LNCipedia <- makeGRangesFromDataFrame(LNCipedia, keep.extra.columns = TRUE)
  }
  
  # Now, LNCipedia should be GRanges object
  if(class(LNCipedia)!="GRanges") {stop("Input object is neither character, data.frame nor GRanges object!\n")}
  
  if(is.null(names(LNCipedia))) {names(LNCipedia) <- factor(mcols(LNCipedia)$name)} # transcript name
  if(!("geneName" %in% names(mcols(LNCipedia)))) {mcols(LNCipedia)$geneName <- gsub(":.*$", "", mcols(LNCipedia)$name)} # gene names without transcript suffix
  
  
  ## prune Lncipedia transcripts to genes (by ":Number"-suffix)
    if(collapseTranscripts2Genes) {
      cat("\nCollapsing lncRNA transcripts to genes\n")
      names(LNCipedia) <- factor(mcols(LNCipedia)$geneName)
      for(i in unique(names(LNCipedia))) {
        tempgr <- range(LNCipedia[names(LNCipedia) == i]) # collapse all transcripts per gene
        start(LNCipedia[names(LNCipedia) == i]) <- start(tempgr) # all transcripts assigned with gene coordinates
        end(LNCipedia[names(LNCipedia) == i]) <- end(tempgr)
      }
  } # end prune transcripts
  
  
 
  
  ##makeLNCipediaExons
  if(makeExonRanges) {
    cat("\nPreparing lncRNA exons.\n")
    gr.exons <- GRanges()
    for(i in 1:length(LNCipedia)) {
      gr.exons <- c(gr.exons, 
                    GRanges(seqnames = seqnames(LNCipedia)[i],
                            ranges = IRanges(start= start(LNCipedia)[i] + as.numeric(unlist(strsplit(as.character(LNCipedia$blockStarts[i]), split=","))), 
                                             end  = start(LNCipedia)[i] + as.numeric(unlist(strsplit(as.character(LNCipedia$blockStarts[i]), split=","))) +
                                                    as.numeric(unlist(strsplit(as.character(LNCipedia$blockSizes[i]), split=","))), 
                                             names = paste0(names(LNCipedia)[i], "_exon",1:LNCipedia$blockCount[i])),
                            strand = strand(LNCipedia)[i]) )
      }

    mcols(gr.exons) <- mcols(LNCipedia)[Rle(1:length(LNCipedia),LNCipedia$blockCount),]
    mcols(gr.exons)$startGene <- start(LNCipedia)[as.vector(Rle(1:length(LNCipedia),LNCipedia$blockCount))]
    mcols(gr.exons)$endGene <- end(LNCipedia)[as.vector(Rle(1:length(LNCipedia),LNCipedia$blockCount))]
    mcols(gr.exons)$exonName <- names(gr.exons)
  
    LNCipedia <- gr.exons
  } # end makeExonRanges 
  
 
  ## append 10 Bases to each side of the exons
    if(addBases>0) {
      cat("\n", addBases, "Bases added to each side of the ranges")
      start(LNCipedia) <- start(LNCipedia) - addBases
      end(LNCipedia) <- end(LNCipedia) + addBases
    } # end append bases
    
   
 ## reduce ranges, e.g. final regions for sequencing
 # reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.  
  reduceRanges = FALSE
  if(reduceRanges) {
    strand(LNCipedia) <- "*"
    LNCipedia <- reduce(LNCipedia)
    cat("\nTotal regions length with collapsed strands for NGS:", sum(width(LNCipedia)), "bp")
  } # end reduce ranges  
  
 
  return(LNCipedia)
}

