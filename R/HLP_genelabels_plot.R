#' Plot genes as arrows
#' 
#' Plotting of genes and exones incl. gene labels
#' 
#' Shows an arrow (orientation = strand) and the according gene symbol in the plot.
#' Can also mark exons with a thickening of the gene line.
#'
#' @param genes a dataframe of genes to plot, with columns "start", "end", "mean", "name", "strand".
#'               \code{strand = -1} means reverse strand, everything else forward.
#' @param levels integer giving desired number of layers to plot multiple genes.
#' @param exons a data frame containing columns "exon_chrom_start" and "exon_chrom_end" with start and end positions of exons.
#' @param xstart an integer number defining the smallest base position of the plot packet (~xlim).
#' @param xstop an integer number defining the largest base position of the plot packet (~xlim).
#' @param ystart numeric y position of the line marking the gene.
#' @param ygenesize numeric y extentsion.
#' @param scale.factor character expansion value as in e.g 'xyplot'. Is also used for other scaling calculations (e.g. exon size, lwd) here.
#' @param arrow.scale.factor character expansion value for arrows.
#' @param truncate When TRUE, gene start / end positions will be truncated to xstart and xstop. Exons with exon_chrom_start > xstop
#'         or exon_chrom_end < xstart are dicarded (that are those that lie completely beyond xstart and xstop).
#' @param genecol character giving gene color.
#' @param genenamecol character giving color of gene names. 
#' 
#' @author Frank Ruehle
#' 
#' @export 


regionalplot.genelabels <- function(
  genes,
  levels,
  exons,
  xstart,
  xstop,
  ystart,
  ygenesize,
  scale.factor,
  arrow.scale.factor = scale.factor,
  truncate = TRUE,
  genecol = "forestgreen", 
  genenamecol = "black"
) 

{
   if(nrow(genes) > 0) {
    neworder <- order(genes$mean)
    genes <- genes[neworder, ]
   # genes$ypos <- rep(-1 * (ygenesize/3.5 + ygenesize * 0:(levels-1)) + ystart, nrow(genes))[1:nrow(genes)] 
    genes$ypos <- rep(-1 * ( ygenesize * 0:(levels-1)) + ystart, nrow(genes))[1:nrow(genes)] 
    
    
    # transfer y position to exon
    if(!is.null(exons) && nrow(exons) > 0)
      exons <- merge(exons, genes, by.x = c("chromo", "genestart", "geneend"), by.y = c("chromo", "start", "end"))
    if(truncate) {
      if(!is.null(exons) && nrow(exons) > 0)
        exons <- exons[exons$end > xstart & exons$start < xstop, ]
      genes[genes$start < xstart, "start"] <- xstart
      genes[genes$end > xstop, "end"] <- xstop
      # do not show gene name when mean is outside window
      # genes[genes$mean < xstart | genes$mean > xstop, "name"] <- ""
    }
    

  # gene arrows - forward and reverse
  if(nrow(genes) > 0)
      arrows(
        x0 = ifelse(genes$strand == -1, genes$end, genes$start), # if reverse, else forward
        y0 = genes$ypos,
        x1 = ifelse(genes$strand == -1, genes$start, genes$end),
        y1 = genes$ypos,
        col = if(length(genecol)==1) {genecol} else {genecol[neworder]},
        length = 0.05 * arrow.scale.factor^2,
        lwd = 2*arrow.scale.factor
      )
    

    # y position of exons == gene arrow positions
    if(!is.null(exons) && nrow(exons) > 0) {
      rect(
        xleft = exons$start,
        ybottom = exons$ypos - 0.2 / scale.factor^2,
        xright = exons$end,
        ytop = exons$ypos + 0.2 / scale.factor^2,
        border = "black",
        col = "black"
      )
    }
  
    # gene names
    text(
      x = genes$mean,
      y = genes$ypos,
      labels = genes$name,
      col = if(length(genenamecol)==1) {genenamecol} else {genenamecol[neworder]},
      cex = 0.6 * scale.factor,
      pos = 1, offset= 0.7
      
    )
  }
}
