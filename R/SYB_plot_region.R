
#' Plot p-values in genomic context
#' 
#' \code{plot.region} reads p-value data e.g. from association analysis and prepares a regional plot of a given chromosomal region of interest.
#' 
#' Up to 5 dataframes can be commited in \code{data} and are plotted in one diagram. 
#' If functional information for variants is available, respective variants which fulfill the 
#' regular expressionn in \code{EFFECT2highlight} are highligted in red. Additionally, variants 
#' given in \code{variant2highlight} are highligthed by filled symbols (e.g. the leading SNP of interest).
#' If given, recombination rates for that region are added to the plot with a separate y-axis.
#' Gene information for the specified region is downloaded from biomaRt and/or LNCipedia and is plotted 
#' beneath the diagram. 
#' Modified graphical parameters are resetted at the end of the function. Nevertheless, 
#' this function is not compatible with using \code{par(mfrow())} for multiple plots.
#' 
#' @param region Character with region of interest of form "chr1:20000-30000" or "chr1:20000".
#' @param window numeric with window size in bp to plot. Only applied when a single basepair position 
#' is given in \code{region}. The window is centered around this position.
#' @param title character with title to be used in plot.
#' @param data named list of dataframes containing p-vales to be plotted. Required columns of each data frame are "CHR", "POS" and "P".
#' An additionally column "EFFECT" with functional characterisation of the locus may be given optionally. 
#' @param variant2highlight character vector with variant to be highligthed as filled symbols. For this, an additionally
#' column \code{ID} is required within \code{data}. Omitted if NULL, If "centered", the centered SNP is highlighted if available
#' (applicable e.g. if \code{region} is of form "chr1:20000"). Vertical lines are added for highligted SNPs or to the center
#' of the plot, respectively.
#' @param EFFECT2highlight character vector with regular expressions (case insensitive). If an \code{EFFECT} column with 
#' functional annotation is given within a dataset, variants containing any of these expressions are highlighted in red.
#' @param recombination.rate character with path to file or to folder containing recombination rates to be plotted. 
#' Alternatively, a dataframe object can be supplied. Omitted if NULL.
#' @param biomaRt biomaRt object to be used for gene annotation. If NULL, biomaRt annotation is skipped. 
#' @param hgnc.symbols.only logical. If TRUE, only genes plotted with annotated HGNC Symbol. If FALSE, non-annotated
#' genes in the plot are labeled with Ensembl gene id if available.
#' @param LNCipedia character with path to LNCipedia bed file. If given lncRNA genes are added to the plot. Omitted if NULL
#' @param numberOfRowsForGenePlotting numeric number of rows used for plotting genes
#' 
#' @return no value returned. Figure is plotted in current graphics device.
#' 
#' @author Frank Ruehle
#' 
#' @export plot.region





plot.region <- function (region, 
                         window=50000, 
                         title = NULL,
                         data = list(),
                         variant2highlight = "centered",
                         EFFECT2highlight = c("miss", "frame", "splice", "start", "stop"),
                         recombination.rate = "V:/proj/public_data/genetic_map_human/NCBI_genetic_map_GRCh37",
                         biomaRt,
                         hgnc.symbols.only = TRUE,
                         LNCipedia = NULL,
                         numberOfRowsForGenePlotting = 3
                          ) {

 
  
  # load required packages. 
  pkg.bioc <- c("GenomicRanges")
  pks2detach <- attach_package(pkg.bioc = pkg.bioc)
  
  

  cex.plot <- 1.5
  cex.legend <- 1.3
  
  # define plot characters for different data sets
  if(length(data)>5) stop("Maximum of 5 datasets allowed.")
  dottype <- c(21,4,24,22,23)
  names(dottype)[1:length(data)] <- names(data)
  dottypeHighlighted <- c(21,8,24,22,23)
  names(dottypeHighlighted)[1:length(data)] <- names(data)
  
  # Named character vector with regular expressions for gene biotype color coding 
  # (name = color, value = regexp) in order of priority high to low,
  # i.e. if multiple biotypes available per gene, first biotype in vector is used.
  # gray = "other" is appended to the vector for all remaining biotypes not found by the reg exp.
  gene.color.coding <- c(black ="protein_coding", orange = "miRNA", forestgreen ="lincRNA|antisense", 
                         brown ="snRNA", blue ="pseudogene")
  
  
  ## extract genomic coordinates from region
  chr <- gsub(":.*$", "", region) # extract chromosome 
  chr <- gsub("chr", "", chr)
  
  range <- gsub("^.*:", "", region) # extract chromosomal range from region
  rangeStartPos <- as.numeric(gsub("-.*$", "", range))
  
  
  if(grepl("-", region)) { # define start and stop positions
    start <- rangeStartPos
    end <- as.numeric(gsub("^.*-", "", range))
  } else {
    start <- rangeStartPos - 0.5* window
    end <- rangeStartPos + 0.5* window
  }
  
  
  range.gr <- GRanges(seqnames = chr, # selected region as genomic range
                      ranges = IRanges(start = start, end = end, names = "selected_range"),
                      strand = "*")
  
  
  #### Get recombination rate
  if (!is.null(recombination.rate)) {
   cat("Get recombination rate\n") 
      if(!is.data.frame(recombination.rate)) {
        
        if(file.info(recombination.rate)$isdir) { # find file for respective chromosome
          recombination.rate <- file.path(recombination.rate, grep(paste0("chr", chr, "\\D"), list.files(recombination.rate), ignore.case = T, value = T)[1])
        } 
        
        rate <- read.table(recombination.rate, header=T, sep="\t",stringsAsFactors = F)
      } else {rate <- recombination.rate}  
      
      column.chr <- grep("chr", names(rate), value=T, ignore.case = T)
      column.pos <- grep("pos", names(rate), value=T, ignore.case = T)
      column.rate <- grep("rate", names(rate), value=T, ignore.case = T)
      
      if(length(column.chr) == 1) { rate <- rate[gsub("chr", "", rate[,column.chr], ignore.case=T) == chr,] }
      rate <- rate[rate[,column.pos] > start & rate[,column.pos] < end, ]
    
    rate.legend <- "Rec."; rate.lty <- 1; rate.pch <- NA ; rate.col <- "blue"; rate.lwd=2
  
  } else {
    rate.legend <- NULL; rate.lty <- NULL; rate.pch <- NULL; rate.col <- NULL; rate.lwd=NULL 
  }
  
  #### get genes to plot
    if (!is.null(biomaRt)) {
    cat("Download gene coordinates from biomart\n")
    
      attrGenes = list(id = "entrezgene", ensemblID = "ensembl_gene_id", name = "hgnc_symbol", chr = "chromosome_name", 
                     startbp = "start_position", endbp = "end_position", strand = "strand", type='transcript_biotype')
      
      genes <- getBM(attributes = attrGenes, filters = c("chromosome_name", "start", "end"),
                           values = list(chr, start, end), mart = biomaRt)
      
      
      if(hgnc.symbols.only) { # plot only genes with annotated hgnc symbol
        genes <- genes[!is.na(genes$hgnc_symbol) & genes$hgnc_symbol!="",]
      } else { # use ensembl gene id as gene label if no symbol available
        genes$hgnc_symbol <- ifelse(is.na(genes$hgnc_symbol), genes$ensembl_gene_id, genes$hgnc_symbol)
      }
      
      # apply color code to gene biotypes. Gray is default for all other types
      genes$col <- rep("gray", nrow(genes))   
      if(nrow(genes) >0) {
        for (g in rev(gene.color.coding)) {
          genes[grepl(g, genes$transcript_biotype, ignore.case = T), "col"] <- names(gene.color.coding)[gene.color.coding == g]
          }
      }
      # define biotype priority 
      # e.g.c("protein_coding", "processed_pseudogene", "lincRNA", "miRNA", "snRNA", "retained_intron", "nonsense_mediated_decay")
      gene.color.coding <- c(gene.color.coding, gray= "other") # append color entry for all remaining biotypes 
      
      genes$col <- ordered(genes$col, levels = names(gene.color.coding))
      # order gene biotypes by priority and remove gene duplicates with lower biotype priority.
      genes <- genes[order(genes$col),]
      genes <- genes[is.na(genes$hgnc_symbol) | !duplicated(genes$hgnc_symbol),]
      
      # Prepare genes dataframe for gene plotting with regionalplot.genelabels(). 
      genes <- data.frame(start=genes$start_position, end=genes$end_position, 
                          mean=genes$start_position+(genes$end_position-genes$start_position)/2,
                          name=genes$hgnc_symbol, strand=genes$strand, chrom=genes$chromosome_name,
                          col=genes$col)
  
  }  ### end genes biomaRt
  
  
  ### get LNCipedia coordinates
  if (!is.null(LNCipedia)) {
    cat("Get lncRNA coordinates from LNCipedia\n") 
    if (is.character(LNCipedia)) {
      LNCipedia <- read.table(LNCipedia, header=F, stringsAsFactors = F)
    }
    
    colnames(LNCipedia) <- c("chrom", "start", "end", "name", "score", "strand", "thickStart", 
                           "thickEnd", "itemRgb", "blockCount",  "blockSizes", "blockStarts")
    LNCipedia$chrom <- gsub("chr", "", LNCipedia$chrom, ignore.case=T)
    #LNCipedia <- LNCipedia[LNCipedia$chrom == chr, ]
    
    LNCipedia.ranges <- makeGRangesFromDataFrame(LNCipedia, keep.extra.columns = TRUE)
    names(LNCipedia.ranges) <- factor(gsub(":.*$", "", elementMetadata(LNCipedia.ranges)$name)) # gene names without transcript suffix
  
    LNCipedia.ranges <- subsetByOverlaps(LNCipedia.ranges, range.gr)
    
    if(length(LNCipedia.ranges) >0) { # check if LNCipedia entries existing in desired region
      
        ## prune Lncipedia transcripts to genes (by ":Number"-suffix)
        LncipediaCollapseTranscripts <- TRUE
        if(LncipediaCollapseTranscripts == T) {
          
          LNCipedia.ranges.genes <- LNCipedia.ranges[0] # initialise new GRange object
          for(i in unique(names(LNCipedia.ranges))) {
             # tempgr <- reduce(LNCipedia.ranges[names(LNCipedia.ranges) == i], drop.empty.ranges=FALSE, min.gapwidth=1L) # collapse overlapping transcripts
             tempgr <- range(LNCipedia.ranges[names(LNCipedia.ranges) == i]) # collapse all transcripts per gene
             mcols(tempgr)$name <- i
             LNCipedia.ranges.genes <- c(LNCipedia.ranges.genes, tempgr)
           }
  
          LNCipedia.ranges <- LNCipedia.ranges.genes
          names(LNCipedia.ranges) <- 1:length(LNCipedia.ranges)
        } # end prune transcripts
  
      start(LNCipedia.ranges) <- start(LNCipedia.ranges) +1 # conformity with ensembl annotation
          
      LNCipedia <- as.data.frame(LNCipedia.ranges)  
      names(LNCipedia)[names(LNCipedia)=="seqnames"] <- "chrom"
      LNCipedia$strand <-  sapply(as.character(LNCipedia$strand), function(x) {switch(x, "+" = 1, "-" = -1)}) # recode +/- to 1/-1
      LNCipedia$col <- "forestgreen"
      LNCipedia$mean <- LNCipedia$start+(LNCipedia$end-LNCipedia$start)/2
      
      if (!is.null(biomaRt)) {
        genes <- merge(genes, LNCipedia, by=c("chrom", "start", "end",  "name", "strand", "col", "mean"), all=T)
        genes <- genes[order(genes$width), ] # by this, LNCipedia entries are prefered against ensembl annotation in later duplicate removal
        } else {
          genes <- LNCipedia}
      
    } else {LNCipedia <- NULL} # end if length(LNCipedia.ranges) >0. Skip LNCipedia otherwise
  } ## end LNCipedia
  
  if (!is.null(biomaRt) || !is.null(LNCipedia)) { # purify gene table if applicable
    genes <- genes[!duplicated(genes$name), ]
    genes <- genes[!duplicated(genes[ ,c("chrom", "start", "end", "strand")]), ]
    genes <- genes[order(genes$start), ]
    print(genes[,c("chrom", "start", "end", "name", "strand", "col")]) # print genes included in plot
  }

  # purify supplied datasets and define plot y-Dimensions as max of p-vales from all plotted datasets 
  stopy <- 1 # initialise stop coordinate of y-axis (used if no datasets available but recombination rate has to be plotted).
  for (d in names(data)) { 
    data[[d]]$CHR <- as.character(data[[d]]$CHR)
    data[[d]] <- data[[d]][!is.na(data[[d]]$CHR), ]
    data[[d]] <- data[[d]][!is.na(data[[d]]$POS), ]
    data[[d]] <- data[[d]][!is.na(data[[d]]$P), ]
    data[[d]] <- data[[d]][data[[d]]$CHR == chr, ]
    data[[d]] <- data[[d]][data[[d]]$POS > start & data[[d]]$POS < end, ]
    if (length(na.omit(data[[d]]$P)) > 0) {stopy.temp <- ceiling(max(-log10(data[[d]]$P), na.rm = T))} else {stopy.temp <- 0}
   stopy <- max(stopy, stopy.temp, na.rm = T)
  }
  
  
  
  
  ########################################
  #par(mfrow=c(2,1)) # plot in 2 rows and 1 column 
  ## start plotting
  
  # backup current graphical parameter
  par.old.fig <- par("fig")
  par.old.mar <- par("mar")
  par.old.oma <- par("oma")
  
  par(fig=c(0, 1, 0.4, 1))
  par(mar=c(1,5,3,5)) # number of margin lines: bottom, left, top, right   
  

  plot(c(start, end), c(0, ceiling(stopy)), type = "n",  xlab = "", ylab = "-log10(p)", axes=T, 
       main = title, cex= cex.plot,  cex.axis= cex.plot, cex.lab = cex.plot, cex.main= cex.plot, lwd=2)   
  
  

  
  # initialise for legend if needed
  if(!is.null(EFFECT2highlight)) { # initialise for legend if needed
    func.legend <- "Func. var."; func.lty <- F; func.pch <- dottype[1]; func.col= "red"; func.lwd=2
      } else {func.legend <- NULL; func.lty <- NULL; func.pch <- NULL; func.col= NULL; func.lwd= NULL}
  
  
  for (d in names(data)) { # plot each dataset
    
    ### plot p values from data as circles. Functional variants in red.
    data.plot <- data[[d]]
    
    data.plot$dot.col <- rep("black", nrow(data.plot))
    #dot.col <- "black"
    if(!is.null(EFFECT2highlight)) {
      if("EFFECT" %in% names(data.plot)) {
        data.plot$dot.col <- ifelse(grepl(paste(EFFECT2highlight,collapse="|"), data.plot$EFFECT, ignore.case = T),
                          func.col, "black")
        } else {warning(paste("no EFFECT column found in", d))}
    }
    
    points(data.plot$POS, -log10(data.plot$P), cex=1, lwd= 2, pch=dottype[d] , col= data.plot$dot.col)
    
    # plot again SNPs to highlight
    if(!is.null(variant2highlight)) {
      
      if (variant2highlight =="centered") {
        data.highlight <- data.plot[data.plot$POS == rangeStartPos, , drop=F]
        abline(v=rangeStartPos, lty=2, col="darkgrey") # prints line even if there is no centered SNP in the dataset
        }
        else {data.highlight <- data.plot[data.plot$ID %in% variant2highlight, , drop=F]}
      
      if(nrow(data.highlight)>0) {
        points(data.highlight$POS, -log10(data.highlight$P), cex=1, lwd= 2, pch=dottypeHighlighted[d], 
               col= data.highlight$dot.col, bg= data.highlight$dot.col) # plot highligted variant
        abline(v=data.highlight$POS, lty=2, col="darkgrey") # print lines through SNP(s)
        }
    }
    
  }
  
  # plot legend for upper panel
  # legend for functional annotation and recombination rate is added if needed only
  legend("topleft", bty= "n", cex = cex.legend,    # ncol = 2, 
         title= "",
         legend= c(names(data), func.legend, rate.legend),
         lty=c(rep(F, length(data)), func.lty, rate.lty), 
         pch=c(dottype[names(data)], func.pch, rate.pch),
         lwd= c(rep(2, length(data)), func.lwd, rate.lwd),
         col=c(rep("black", length(data)), func.col, rate.col)
        )
  # add region coordinates as legend title. This done separetly because otherwise the legend itself
  # is centered below the (long) title and not left-adjusted.
  legend("topleft", legend="", title= paste0("chr", chr, ":", start, "-", end), bty= "n", cex = cex.legend)


  ##### Include racombination rate as Lineplot
  if(!is.null(recombination.rate)) {
    
    rate.max <- max(rate[,column.rate])
 
    # cat("\n max rate: ", max(rate[,column.rate]))
     par(new = T)
     plot(c(start, end), c(0, rate.max),  type = "n", xlab = "", ylab = "", axes=F)   
     axis(side=4, cex.axis= cex.plot) # , padj=-1.5
     mtext(text= "Recombination Rate (cM/Mb)", side=4,  line=3, cex= cex.plot, col="black", lwd=2) #, padj=0.3


     points(rate[,column.pos], rate[,column.rate], type = "l", lwd=1.5, col= rate.col)
    
  } 


  # plot genes 
  if (!is.null(biomaRt) || !is.null(LNCipedia)) {
    
    par(fig=c(0,1,0,0.4), new=TRUE)
    par(mar=c(1,5,1,5)) # number of margin lines: bottom, left, top, right   
    ymax <- (numberOfRowsForGenePlotting / 2) 
    plot(c(start, end), c(0, ymax), type = "n",  xlab = "", ylab = "", axes=F)   
    
    
    # call plot function
    regionalplot.genelabels(
      genes=genes,
      levels=numberOfRowsForGenePlotting,
      exons = NULL,
      xstart=start,
      xstop=end,
      ystart=ymax -0.1,
      ygenesize= 0.5,
      scale.factor=2,
      truncate = F,
      genecol = as.character(genes$col),
      genenamecol = as.character(genes$col)
    )
  

    # prepare legend for gene biotypes (print relevant biotypes only)
    gene.color.coding.applied <- gene.color.coding[names(gene.color.coding) %in% genes$col] 
    gene.color.coding.applied <- gene.color.coding.applied[order(nchar(gene.color.coding.applied))] # order by character length
    gene.color.coding.applied <- gsub("[\\^\\$]", "", gene.color.coding.applied) # remove special characters from legend
    # An additional transparent plot is called for the legend to utilize the margin area
    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    
    legend("bottomleft", bty= "n", cex = cex.legend, 
           legend= gene.color.coding.applied,
           pch= "-", lwd= 3,
           col= names(gene.color.coding.applied)
         )
    
  }
  
  # re-set graphical parameter
  par(fig = par.old.fig)
  par(mar = par.old.mar)
  par(oma = par.old.oma)
  

# Detaching libraries not needed any more
detach_package(unique(pks2detach))

} # end of function def



