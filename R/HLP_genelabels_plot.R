#' Plot genes as arrows
#' 
#' Plotting of genes and exones incl. gene labels
#' 
#' Shows an arrow (orientation = strand) and the according gene symbol in the plot.
#' Can also mark exons with a thickening of the gene line.
#'
#' @param genes a dataframe of genes to plot, with columns \code{start}, \code{end}, \code{name}, \code{strand}.
#'               \code{strand = -1 or -} means reverse strand, everything else forward.
#' @param levels integer giving desired number of layers to plot multiple genes. Genes are assigned to layers via gene start position.
#' If NULL, an additional column \code{ypos} is expected in \code{genes} giving the y-coordinates for the genes to plot.
#' @param exons a data frame containing columns \code{exon_chrom_start} and \code{exon_chrom_end} with start and end positions of exons.
#' Additionally, columns \code{start} and \code{end} with respective gene coordinates needed as in \code{genes}. 
#' Exon plotting omitted if NULL.
#' @param xstart an integer number defining the smallest base position of the plot packet.
#' @param xstop an integer number defining the largest base position of the plot packet.
#' @param ystart numeric y position of the line marking the gene.
#' @param ygenesize numeric y extentsion.
#' @param scale.factor character expansion value as in e.g \code{xyplot}. Is also used for other scaling calculations (e.g. exon size, lwd) here.
#' @param arrow.scale.factor character expansion value for arrows.
#' @param truncate When TRUE, gene start / end positions will be truncated to xstart and xstop. Exons with exon_chrom_start > xstop
#'         or exon_chrom_end < xstart are dicarded (that are those that lie completely beyond xstart and xstop).
#' @param genecol character giving gene color.
#' @param genenamecol character giving color of gene names.
#' 
#' @return no value returned. Genes are plotted in active graphics device.
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

     if(!all(grepl("1", genes$strand))) {
     genes$strand <-  sapply(as.character(genes$strand), function(x) {switch(x, "+" = 1, "-" = -1)}) # recode +/- to 1/-1 if necessary
     }
     
     # calculate mean value of gene poistion if not supplied
     if(!"mean" %in% names(genes)) {
       genes$mean <- genes$start+(genes$end-genes$start)/2 
     }
     
     # re-order gene according to coordinates 
     neworder <- order(genes$mean)
     genes <- genes[neworder, ]
     
    if(!is.null(levels)) { # otherwise, column "ypos" expected in genes
      ## genes$ypos <- rep(-1 * (ygenesize/3.5 + ygenesize * 0:(levels-1)) + ystart, nrow(genes))[1:nrow(genes)] # ADAMTS13
      genes$ypos <- rep(-1 * ( ygenesize * 0:(levels-1)) + ystart, nrow(genes))[1:nrow(genes)] 
    }

    if(!"ypos" %in% names(genes)) {stop("no y-coordinate for genes specified!\n")}
      


    # transfer y position to exon
    if(!is.null(exons) && nrow(exons) > 0)
      exons <- merge(exons, genes, by.x = c("start", "end"), by.y = c("start", "end")) # merge gene and exon data
    if(truncate) {
      if(!is.null(exons) && nrow(exons) > 0)
        exons <- exons[exons$end > xstart & exons$start < xstop, ]
      genes[genes$start < xstart, "start"] <- xstart
      genes[genes$end > xstop, "end"] <- xstop
      # do not show gene name when mean is outside window
      # genes[genes$mean < xstart | genes$mean > xstop, "name"] <- ""
    }
    

  # plot gene arrows - forward and reverse
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
    

    # plot exons as rectangles. y position of exons == gene arrow positions
    if(!is.null(exons) && nrow(exons) > 0) {
      rect(
        xleft = exons$exon_chrom_start,
        ybottom = exons$ypos - 0.2 / scale.factor^2,
        xright = exons$exon_chrom_end,
        ytop = exons$ypos + 0.2 / scale.factor^2,
        border = "black",
        col = "black"
      )
    }
  
    # plot gene names
    text(
      x = genes$mean,
      y = genes$ypos,
      labels = genes$name,
      col = if(length(genenamecol)==1) {genenamecol} else {genenamecol[neworder]},
      cex = 0.5 * scale.factor,
      pos = 1, offset= 0.7
      
    )
  }
} # end function def




###########################################
#' Non-overlapping gene layers for plotting
#' 
#' Determine gene layers for plotting of genes in regional plot without overlap
#' 
#' This is a help function for regional plots. For a given set of genes, the functions
#' determines separate layers necessary to plot the genes without overlapping. 
#' 
#' @param gr GRanges object containing gene coordinates.
#' @param ystart numeric start y-coordinate for first layer. 
#' @param ysize numeric added to each new layer.
#' @param minDistance numeric minimum distance for two genes to be considered as non-overlapping.
#' If \code{minDistance} < 1, it is interpreted as ratio of the width of the full \code{gr} object.
#' @param gene_name_column character giving column name of gene names in \code{gr}. 
#' If not NULL, overlapping gene names are considered additionally to gene coordinates.
#' @param prot_domain_column character giving column name in \code{gr} with logical 
#' values indicating genes to be plotted with protein domains. For those, needed space in layer
#' below the respective gene layer is considered when determining non-overlapping gene layers.
#' Omitted if NULL.
#' @param scale.factor numeric scale factor as used in \code{regionalplot.genelabels}.
#' Required to determine character extension of gene names (if \code{gene_name_column} is not NULL).
#' @param units character indicating units for \code{strwidth}-function, either \code{"user"}, \code{"inches"} 
#' or \code{"figure"}. Measurements in \code{"user"} units are only available after \code{plot.new} has been called.
#' Required to determine character extension of gene names (if \code{gene_name_column} is not NULL).
#' @param sortGeneStart logical. If TRUE, gene levels are sorted according to minimum start coordinate.  
#' Useful to keep free space for legend in a plot.
#' 
#' @return GRanges object with column \code{ypos} added to meta data.
#' 
#' @author Frank Ruehle
#' 
#' @export 


determineNonOverlapGenelayers <- function(gr, ystart = 0, ysize = 1, minDistance = 0.1, 
                                          gene_name_column = NULL, prot_domain_column = NULL, 
                                          scale.factor=1, units="user",
                                          sortGeneStart= FALSE) {

  gr$ypos <- ystart # initialise first layer
  
  if(is.null(prot_domain_column)) {
    prot_domain_column <- "prot_dom_hlp"
    mcols(gr)[,prot_domain_column] <- FALSE
    }
  
 
  # make additional GRanges objext with string width of gene names
  sw <- gr
  
  if(!is.null(gene_name_column)) { # if gene_name_column = NULL, gr (=sw) just used twice in while condition
  mcols(sw) <- mcols(sw)[, c(gene_name_column, prot_domain_column, "ypos")]
  sw$mean <- start(sw) + (end(sw) - start(sw))/2
  start(sw) <- sw$mean - 0.5* strwidth(sw$name, units = units, cex = 0.5 * scale.factor)
  end(sw) <- sw$mean + 0.5* strwidth(sw$name, units = units, cex = 0.5 * scale.factor)
  }


  # interpret minDistanve either as ratio or as absolute number
  if(minDistance < 1) {minDistance <- minDistance * width(range(gr, ignore.strand=T))}
  
  if(length(gr)==1) {return(gr)} # if just one gene present

  for (i in 2:length(gr)){

    level.hlp <- ystart
    
    temp.gr <- gr[i]
    start(temp.gr) <- start(temp.gr) - minDistance
    end(temp.gr) <- end(temp.gr) + minDistance
   
    temp.sw <- sw[i]
     
    while( # condition for gene width and string width
      !is.na(findOverlaps(range(temp.gr), reduce(gr[1:i-1][gr[1:i-1]$ypos==level.hlp], ignore.strand=T), select="first")) | # overlap gene arrows?
      !is.na(findOverlaps(range(temp.gr), reduce(gr[1:i-1][gr[1:i-1]$ypos== (level.hlp - ysize * mcols(temp.gr)[,prot_domain_column])], ignore.strand=T), select="first")) | # overlap protein domains to plot?
      !is.na(findOverlaps(range(temp.sw), reduce(sw[1:i-1][sw[1:i-1]$ypos==level.hlp], ignore.strand=T), select="first"))   # overlap name characters?
      
       ) {
        gr[i]$ypos <- gr[i]$ypos + ysize
        sw[i]$ypos <- sw[i]$ypos + ysize
        level.hlp <- level.hlp + ysize
      }
  }  

  
  if(isTRUE(sortGeneStart)) {
    # reorder gene layers to keep space free for legend

    minStart <- tapply(ranges(gr), gr$ypos, function(x) {min(start(x))}) # gene layer with minimum start coordinate
    newOrder <- sort(minStart, decreasing = T)

    gr$ypos <- factor(gr$ypos) 
    levels(gr$ypos) <- names(newOrder) # assign factor levels in requested order
    gr$ypos <- as.numeric(as.character(gr$ypos))
   }  

  mcols(gr)[, names(mcols(gr))=="prot_dom_hlp"] <- NULL # remove temporary column from meta data
  
  
  return(gr)
} # end function def





