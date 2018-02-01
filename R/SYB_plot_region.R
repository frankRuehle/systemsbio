
#' Plot p-values in regional genomic context
#' 
#' \code{plot.region} reads p-value data e.g. from association analysis and prepares a regional plot of a given chromosomal region of interest.
#' 
#' Up to 5 dataframes can be committed in \code{data} and are plotted in one diagram. 
#' If functional information for variants is available, respective variants which fulfill the 
#' regular expression in \code{EFFECT2highlight} are highlighted by color. Additionally, variants 
#' given in \code{variant2highlight} are highlighted by filled symbols, text annotation and vertical lines 
#' (e.g. for the leading SNP of interest).
#' If given, recombination rates for that region are added to the plot using a separate y-axis.
#' Gene information for the specified region is downloaded from biomaRt and/or LNCipedia and is plotted 
#' beneath the diagram. Genes can be selected to include corresponding protein domain data for plotting. 
#' Modified graphical parameters are resetted at the end of the function. 
#' Nevertheless, this function can not be used \code{par(mfrow())} for multiple plots.
#' 
#' @param region Character with region of interest of form \code{chr1:20000-30000} or \code{chr1:20000}.
#' @param region_ext numeric with region size extension in bp to plot. Half of the extension is added 
#' to both sides of the given \code{region}.
#' @param title character with title to be used in plot.
#' @param data named list of dataframes containing p-vales to be plotted. Required columns of each 
#' data frame are \code{"CHR"}, \code{"POS"} and \code{"P"}.
#' An additionally column \code{"EFFECT"} with functional characterisation of the locus may be given optionally. 
#' @param lines_pvalue_threshold named character with p-values to be plotted as threshold lines. 
#' Line color is given by vector names (e.g. \code{lines_pvalue_threshold = c(blue=0.05, red=0.01)}).
#' @param variant2highlight character vector with variants to be highlighted as filled symbols. For this, an additionally
#' column \code{ID} is required within \code{data}. If the vector contains color names or effect names, all according variants 
#' are also highlighted. Numbers \code{>=1} are interpreted as BP position to highlight while numbers \code{<1} are 
#' interpreted as p-value threshold with all SNPs highlighted with \code{p < threshold}.
#' If \code{variant2highlight = "centered"}, the centered SNP is highlighted if available 
#' (applicable if \code{region} is of form \code{"chr1:20000"}). 
#' Vertical lines are added for the highlighted SNPs. Omitted if \code{NULL}.
#' @param EFFECT2highlight named character vector with regular expressions (name = color, value = regexp) 
#' in order of priority low to high. Regular expressions is case insensitive.
#' If an \code{EFFECT} column with functional annotation is given within a dataset, variants with functional annotation 
#' corresponding to these expressions are highlighted by colors given as vector names.
#' If no \code{EFFECT} column is given, exonic SNPs can be highlighted according to overlapping gene exons.
#' For this, an \code{EFFECT} column is created if not yet existing and \code{exonic} as well as the respective 
#' gene biotype is appended to the entries of the \code{EFFECT} column for exonic SNPs. 
#' @param recombination.rate character with path to file or to folder containing recombination rates to be plotted. 
#' Alternatively, a dataframe object can be supplied. Omitted if \code{NULL}.
#' @param biomaRt biomaRt object to be used for gene annotation. If \code{NULL}, biomaRt annotation is skipped. 
#' @param hgnc.symbols.only logical. If \code{TRUE}, only Ensemble genes plotted with annotated HGNC Symbol. 
#' If \code{FALSE}, non-annotated genes in the plot are labeled with Ensembl gene id if available.
#' @param LNCipedia character with path to LNCipedia bed file to plot lncRNA genes. Omitted if \code{NULL}.
#' @param gene.color.coding named character vector with regular expressions for gene biotype color coding
#' (name = color, value = regexp) in order of priority low to high, i.e. if multiple biotypes available per gene, 
#' the last biotype in the vector is used. Regular expressions are case insensitive. 
#' \code{gray = "other"} is appended to the vector for all remaining biotypes not found by the reg exp.
#' @param numberOfRowsForGenePlotting numeric number of rows used for plotting genes. If \code{"auto"}, function determines
#' appropriate number of rows itself.
#' @param plot.protein.domains named character vector with file path to protein domain annotation data for a selected gene (Omitted if \code{NULL}). 
#' Vector names are used as gene name of the selected gene (e.g. \code{GeneXY = "filepath_to_protein_data_of_GeneXY"}).
#' Domains are plotted as symbols below the respective gene. The protein length is scaled to length of the plotted gene.
#' Domain positions and width are scaled accordingly. Arrows indicate the respective genomic start and stop positions for each domain. 
#' The respective txt-file may be generated by the function \code{makeDomainsFromExons} and contains the following columns:
#' \itemize{
#'   \item BP_start: genomic start position for corresponding protein domain AA position.
#'   \item BP_end: genomic end position for corresponding protein domain AA position.
#'   \item feature_length: domain/feature length in AA. Used for domain scaling in the plot.
#'   \item protein_length: total protein length in AA. Used for domain scaling in the plot.
#'   \item domain_name_plot: domain/feature name to be plotted.
#'   \item symbol_plot (optional): Shape to be used for plotting (either \code{"ellipse"}, \code{"rectangle"} or \code{"circle"}). 
#'   Default is \code{"ellipse"}.
#'   \item domain_height_extension (optional): height factor for symbol height. These factors are scaled respective to each other. Default is \code{1}.
#'   \item domain_color (optional): color for domain symbol and name to be plotted. Default is \code{"black"}.
#'   \item label_pos (optional): label position in the domain plot may be adjusted in case of overlapping labels 
#'   (Values of 1, 2, 3 and 4 indicate positions below, left, above and right of the domain center coordinates). Default is \code{3}.
#'   \item assignArrows2Gene (optional): Indicate if arrows shall be plotted from genomic coordinates to protein domain (default is \code{TRUE}).
#'   May be set to \code{FALSE} for very short features within other domains, e.g. "active site". 
#'   }
#' @param cex.plot numeric character extension plot axes.
#' @param cex.legend numeric character extension plot legends.
#' @param gene.scale.factor numeric extension factor used for gene and exon plotting.
#' 
#' @return no value returned. Figure is plotted in the current graphics device.
#' 
#' @author Frank Ruehle
#' 
#' @export plot.region





plot.region <- function (region, 
                         region_ext=50000, 
                         title = NULL,
                         data = list(),
                         lines_pvalue_threshold = NULL,
                         variant2highlight = "centered",
                         EFFECT2highlight = c(green="splice", red="miss", red="frame", red="start|stop"),
                         recombination.rate = NULL,
                         biomaRt,
                         hgnc.symbols.only = TRUE,
                         LNCipedia = NULL,
                         gene.color.coding = c(lightgreen ="pseudogene", brown ="snRNA", forestgreen ="ncRNA|antisense",
                                                orange = "miRNA", darkblue ="protein_coding"),
                         numberOfRowsForGenePlotting = "auto",
                         plot.protein.domains = NULL,
                         cex.plot = 1.1,
                         cex.legend = 1,
                         gene.scale.factor = 2
                          ) {

 
  
  # load required packages. 
  pkg.cran <- c("TeachingDemos", "plotrix")
  pkg.bioc <- c("GenomicRanges")
  pks2detach <- attach_package(pkg.cran = pkg.cran, pkg.bioc = pkg.bioc)
  
  
  # define plot characters for different data sets
  if(length(data)>5) stop("Maximum of 5 datasets allowed.")
  dottype <- c(21,4,24,22,23)
  names(dottype)[1:length(data)] <- names(data)
  dottypeHighlighted <- c(21,8,24,22,23)
  names(dottypeHighlighted)[1:length(data)] <- names(data)
  

  
  ## extract genomic coordinates from region
  chr <- gsub(":.*$", "", region) # extract chromosome 
  chr <- gsub("chr", "", chr, ignore.case = T)
  
  region <- gsub("^.*:", "", region) # extract chromosomal range from region
  regionStartPos <- as.numeric(gsub("-.*$", "", region)) # either start pos or center of region
  
  
  if(grepl("-", region)) { # define start and stop positions
    start <- regionStartPos - 0.5* region_ext
    end <- as.numeric(gsub("^.*-", "", region)) + 0.5* region_ext
  } else {
    start <- regionStartPos - 0.5* region_ext
    end <- regionStartPos + 0.5* region_ext
  }
  
  
  region.gr <- GRanges(seqnames = chr, # selected region as genomic range
                      ranges = IRanges(start = start, end = end, names = "selected_range"),
                      strand = "*")
  
  
  #### Get recombination rate
  if (!is.null(recombination.rate)) {
      if(!is.data.frame(recombination.rate)) {
        
        if(file.info(recombination.rate)$isdir) { # if directory, find file for respective chromosome
          recombination.rate <- file.path(recombination.rate, grep(paste0("chr", chr, "\\D"), list.files(recombination.rate), ignore.case = T, value = T)[1])
        } 
        
        rate <- read.table(recombination.rate, header=T, sep="\t",stringsAsFactors = F)
        cat(paste("Load recombination rate from", recombination.rate), "\n")
        } else {rate <- recombination.rate}  
      
      column.chr <- grep("chr", names(rate), value=T, ignore.case = T)
      column.pos <- grep("pos", names(rate), value=T, ignore.case = T)
      column.rate <- grep("rate", names(rate), value=T, ignore.case = T)
      
      if(length(column.chr) == 1) { rate <- rate[gsub("chr", "", rate[,column.chr], ignore.case=T) == chr,] }
      rate <- rate[rate[,column.pos] > start & rate[,column.pos] < end, ]
    
    rate.legend <- "Rec."; rate.lty <- 1; rate.pch <- NA ; rate.col <- "blue"; rate.lwd=2; rate.bg=par("bg")
  
  } else {
    rate.legend <- NULL; rate.lty <- NULL; rate.pch <- NULL; rate.col <- NULL; rate.lwd = NULL; rate.bg=NULL 
  }
  
  
  
  #### get Ensembl genes to plot
    if (!is.null(biomaRt)) {
    cat("Download gene coordinates from biomart\n")
    
      attrGenes = list(ensemblID = "ensembl_gene_id", name = "hgnc_symbol", chr = "chromosome_name", 
                     startbp = "start_position", endbp = "end_position", strand = "strand", type='transcript_biotype')
      
      genes <- getBM(attributes = attrGenes, 
                     filters = c("chromosome_name", "start", "end"),
                     values = list(chr, start, end), mart = biomaRt)

      # get all exons from all transcripts
      attrExons <- c("ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", # exon coordinates
                     "ensembl_gene_id", "start_position", "end_position", 
                    type='transcript_biotype', strand = "strand")  # corresponding gene coordinates for merging
      
      exons <- getBM(attributes = attrExons, 
                     filters = c("chromosome_name", "start", "end"),
                     values = list(chr, start, end), mart = biomaRt)
      
            names(exons)[names(exons)=="start_position"] <- "startGene" # gene start coordinate
            names(exons)[names(exons)=="end_position"] <- "endGene" # gene end coordinate
            names(exons)[names(exons)=="chromosome_name"] <- "chrom" # gene start coordinate
            names(exons)[names(exons)=="transcript_biotype"] <- "biotype" # transcript_biotype
            names(exons)[names(exons)=="ensembl_exon_id"] <- "exon_id" 
            names(exons)[names(exons)=="ensembl_gene_id"] <- "gene_id" 
                      
            exons$strand <-  sapply(as.character(exons$strand), function(x) {switch(x, "1" = "+", "-1" = "-")}) # recode 1/-1 to +/- 
            exon.ranges <- makeGRangesFromDataFrame(exons, keep.extra.columns = T,
                                                 ignore.strand = FALSE,
                                                 start.field = "exon_chrom_start",
                                                 end.field = "exon_chrom_end",
                                                 strand.field = "strand",
                                                 starts.in.df.are.0based = FALSE)
            names(exon.ranges) <- exon.ranges$exon_id
            
            

      if(hgnc.symbols.only) { # plot only genes with annotated hgnc symbol
        cat("Ensemble genes removed without hgnc symbol:\n")    
        print(genes[is.na(genes$hgnc_symbol) | genes$hgnc_symbol!="",])
        
        genes <- genes[!is.na(genes$hgnc_symbol) & genes$hgnc_symbol!="",]
        
            } else { # use ensembl gene id as gene label if no symbol available
            genes$hgnc_symbol <- ifelse(is.na(genes$hgnc_symbol) | genes$hgnc_symbol=="", genes$ensembl_gene_id, genes$hgnc_symbol)
        }
                 
 
      # make GRanges object from genes (required for determining non-overlapping gene layers)
      genes$strand <-  ifelse(grepl("-", genes$strand), "-", "+") # recode -1/+1 to -/+
      gene.ranges <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE,
                                              ignore.strand = FALSE,
                                              start.field = "start_position",
                                              end.field = "end_position",
                                              strand.field = "strand",
                                              starts.in.df.are.0based = FALSE)
      names(mcols(gene.ranges))[names(mcols(gene.ranges)) == "hgnc_symbol"] <- "name"
      names(mcols(gene.ranges))[names(mcols(gene.ranges)) == "ensembl_gene_id"] <- "gene_id"
      names(mcols(gene.ranges))[names(mcols(gene.ranges)) == "transcript_biotype"] <- "biotype"
      names(gene.ranges) <- mcols(gene.ranges)$name
      
    }  ### end genes biomaRt
  
  
  ### get LNCipedia coordinates
  if (!is.null(LNCipedia)) {
    cat("Get lncRNA data from LNCipedia\n") 

    # load LNCipedia bed file
    LNCipedia.ranges <- processLNCipedia(LNCipedia, collapseTranscripts2Genes=F, makeExonRanges=F, addBases= c(1,0)) # conformity with ensembl annotation
 
    LNCipedia.ranges <- subsetByOverlaps(LNCipedia.ranges, region.gr)

    if(length(LNCipedia.ranges) >0) { # check if LNCipedia entries exist in desired region

      mcols(LNCipedia.ranges)$biotype <- "LNCipedia_lncRNA"

      # prune LNCipedia transcripts to genes
      # Must be done BEFORE preparing LNCipedia.exon.ranges. Otherwise "startGene" and "endGene" coordinates 
      # in exon table refere to lncRNA transcripts and not to genes. Will cause problems in assigning exons to genes.
      LNCipedia.ranges <- processLNCipedia(LNCipedia.ranges, collapseTranscripts2Genes=T, makeExonRanges=F) 
      #names(mcols(LNCipedia.ranges))[names(mcols(LNCipedia.ranges)) == "geneName"] <- "gene_id" # for compatibility with exon.ranges

      # prepare LNCipedia exon ranges before pruning transcripts to genes
      LNCipedia.exon.ranges <- processLNCipedia(LNCipedia.ranges, collapseTranscripts2Genes=F, makeExonRanges=T) 
      LNCipedia.exon.ranges <- LNCipedia.exon.ranges[, c("exonName", "geneName", "startGene",  "endGene", "biotype")]
      names(mcols(LNCipedia.exon.ranges))[names(mcols(LNCipedia.exon.ranges)) == "exonName"] <- "exon_id" # for compatibility with exon.ranges
      names(mcols(LNCipedia.exon.ranges))[names(mcols(LNCipedia.exon.ranges)) == "geneName"] <- "gene_id" # for compatibility with exon.ranges

      # Remove duplicate entries from LNCipedia gene table coming from collapsing transcripts to genes.
      # Must be done AFTER preparing LNCipedia.exon.ranges. Otherwise exon coordinate information is lost for the removed transcripts.
      LNCipedia.ranges <- LNCipedia.ranges[!duplicated(LNCipedia.ranges)] 
      
    if (!is.null(biomaRt)) { # merge Ensembl and LNCipedia entries
      gene.ranges <- c(gene.ranges[, c("name", "biotype")], LNCipedia.ranges[, c("name", "biotype")]) # no need for "gene_id"
 
      exon.ranges <- c(exon.ranges, LNCipedia.exon.ranges)

             } else {
               gene.ranges <- LNCipedia.genes
                exon.ranges <- LNCipedia.exon.ranges
                }
    
        
    } else {LNCipedia <- NULL} # end if length(LNCipedia.ranges) >0. Skip LNCipedia otherwise
  } ## end LNCipedia
  

  if (!is.null(biomaRt) || !is.null(LNCipedia)) { # annotate variants and purify gene table if applicable
 
    # apply color code to gene biotypes (Ensembl and Lncipedia). Gray is default for all other types
    mcols(gene.ranges)$col <- rep("gray", length(gene.ranges))   
    if(length(gene.ranges) >0) {
      for (g in gene.color.coding) {
        mcols(gene.ranges)[grepl(g, mcols(gene.ranges)$biotype, ignore.case = T), "col"] <- names(gene.color.coding)[gene.color.coding == g]
      }
    }
    
    # append "gray" to gene.color.coding for all other biotypes
    if(!("other" %in% gene.color.coding)) {gene.color.coding <- c(gray= "other", gene.color.coding)} 

    # order gene biotypes by priority and remove gene duplicates with lower biotype priority.
    gene.ranges <- gene.ranges[order(ordered( mcols(gene.ranges)$col, levels = names(gene.color.coding)), decreasing = T),]
    
    options("showHeadLines"=Inf) # show all ranges
    
    gene.ranges.removed <- gene.ranges[duplicated(gene.ranges)]
    cat("\nDuplicate gene entries removed:\n")
    print(gene.ranges.removed) # print genes excluded from plot
    
    gene.ranges <- gene.ranges[!duplicated(gene.ranges)]
    cat("\nGenes included in plot:\n")
    print(gene.ranges) # print genes included in plot
    
    options("showHeadLines"=NULL) # Revert to default value
    }

  
  
  # purify supplied datasets and define plot y-Dimensions as max of p-vales from all plotted datasets 
  stopy <- 1 # initialise stop coordinate of y-axis (used if no datasets available but recombination rate has to be plotted).
  for (d in names(data)) { 
    data[[d]]$CHR <- as.character(data[[d]]$CHR)
    data[[d]]$CHR <- sub("chr", "", data[[d]]$CHR, ignore.case = T)
    data[[d]] <- data[[d]][!is.na(data[[d]]$CHR), , drop=F]
    data[[d]] <- data[[d]][!is.na(data[[d]]$POS), , drop=F]
    data[[d]] <- data[[d]][!is.na(data[[d]]$P), , drop=F]  # or set missing p-values to 1.  data[[d]][is.na(data[[d]]$P), "P"] <- 1
    data[[d]] <- data[[d]][data[[d]]$CHR == chr, , drop=F]
    data[[d]] <- data[[d]][data[[d]]$POS > start & data[[d]]$POS < end, , drop=F]
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
  
  
  ## plot horizontal threshold lines if selected (color given as names in vector)
  if(!is.null(lines_pvalue_threshold)) {
     abline(h = -log10(lines_pvalue_threshold), lty=2, col= names(lines_pvalue_threshold))
   }
    
  ### plot each dataset
  colors.applied <- character()
  for (d in names(data)) { 
    
    ### plot p values from data. Functional variants in color.
    data.plot <- data[[d]]
    
    data.plot$dot.col <- rep("black", nrow(data.plot)) # initialise black
    if(!is.null(EFFECT2highlight)) {
      
      ## prune EFFECT2highlight for identical colors given in function call
      uniqueEffectcol <- unique(names(EFFECT2highlight))

      EFFECT2highlight.h <- character()
      
      for (i in uniqueEffectcol) {
        EFFECT2highlight.h[i] <- paste(EFFECT2highlight[names(EFFECT2highlight)==i], collapse="|")
      }
      EFFECT2highlight <- EFFECT2highlight.h
      
 
      # assign colors for function to datasets
      if(!("EFFECT" %in% names(data.plot))) {
        data.plot$EFFECT <- "" 
        warning(paste("no EFFECT column found in", d, "-> Highlight exonic SNPs only\n")) 
        }
      
         # find overlaps of SNP positions and exon data (Ensembl and LNCipedia)
         # If overlapping with any exon, merge "exonic" to SNP EFFECT column
         # If overlapping with any exon, merge respective biotype from exon.ranges to SNP EFFECT column
          data.plot$EFFECT <- as.character(data.plot$EFFECT)
          data.plot.ranges <-  GRanges(seqnames = data.plot$CHR,
                                      ranges = IRanges(data.plot$POS, end = data.plot$POS, names = data.plot$ID), 
                                      strand = "*",
                                      P = data.plot$P,
                                      EFFECT = data.plot$EFFECT)
          
          hits <- findOverlaps(data.plot.ranges, exon.ranges, select = "first") # Ensembl first, LNCipedia second
           
          data.plot$EFFECT[!is.na(hits)] <- paste(data.plot$EFFECT[!is.na(hits)], "exonic", sep="_")
          data.plot$EFFECT[!is.na(hits)] <- paste(data.plot$EFFECT[!is.na(hits)], mcols(exon.ranges)$biotype[na.omit(hits)], sep="_")
          
            # find regular expression for functional color coding in EFFECT column
            for(i in uniqueEffectcol) {
            data.plot$dot.col <- ifelse(grepl(EFFECT2highlight[i], data.plot$EFFECT, ignore.case = T),
                              i, data.plot$dot.col)
            }
         
    }
    
    colors.applied <- c(colors.applied, unique(data.plot$dot.col)) # list colors needed for legend
    
    data.plot <- data.plot[order(ordered(data.plot$dot.col, levels = c("black", names(EFFECT2highlight))), decreasing = F),] # plot colored dots at last
    
    points(data.plot$POS, -log10(data.plot$P), cex=1, lwd= 2, pch=dottype[d] , col= data.plot$dot.col)
    
    # plot again SNPs to highlight
    if(!is.null(variant2highlight)) {
      
      if (variant2highlight[1] == "centered") {
        data.highlight <- data.plot[data.plot$POS == round(start+(end-start)/2, 0), , drop=F]
        # abline(v=round(start+(end-start)/2, 0), lty=2, col="darkgrey") # prints line even if there is no centered SNP in the dataset
        }
        else {data.highlight <- rbind(data.plot[data.plot$ID %in% variant2highlight, , drop=F], # highlight SNPs by rsids
                                      data.plot[data.plot$dot.col %in% variant2highlight, , drop=F], # highlight SNPs by color code
                                      data.plot[data.plot$POS %in% variant2highlight, , drop=F]) # highlight SNPs by position
        
               if("EFFECT" %in% colnames(data.plot)) { # highlight SNPs by Effect column if present
                 data.highlight <- rbind(data.highlight, data.plot[grepl(paste(variant2highlight, collapse = "|"), data.plot$EFFECT, ignore.case=T), , drop=F])
               }
          
              variant2highlight.threshold <- as.numeric(na.omit(variant2highlight[suppressWarnings(as.numeric(variant2highlight)) < 1]))
              if (length(variant2highlight.threshold)>0) { # highlight SNPs by p-value threshold
                data.highlight <- rbind(data.highlight, data.plot[data.plot$P < min(variant2highlight.threshold), , drop=F])
              }
        
          data.highlight <- unique(data.highlight)              
       }
      
      if(nrow(data.highlight)>0) {
        # text(data.highlight$POS, -log10(data.highlight$P), labels= data.highlight$ID, cex=cex.legend,  
        #        col= data.highlight$dot.col, pos = 4, offset= 0.4) # plot ID of highligted variant without background shadow      
        abline(v = data.highlight$POS, lty=2, col="darkgrey") # print lines through SNP(s)
        points(data.highlight$POS, -log10(data.highlight$P), cex=1, lwd= 2, pch=dottypeHighlighted[d], 
               col= data.highlight$dot.col, bg= data.highlight$dot.col) # plot highligted variant
         shadowtext(data.highlight$POS, -log10(data.highlight$P), labels= data.highlight$ID, cex=cex.legend, pos = 4, offset= 0.4,
                    col = data.highlight$dot.col, bg = "white", r = 0.2) # plot ID of highligted variant using background shadow      
      }
    }
    
  }

  
  
  #### plot legend for upper panel  
    if(!is.null(EFFECT2highlight)) { # initialise for legend if needed

      # legend for needed functional annotation only        
      EFFECT2highlight <- EFFECT2highlight[names(EFFECT2highlight) %in% colors.applied]
           if(length(EFFECT2highlight) > 0) { 
              func.legend <- EFFECT2highlight 
              func.lty <- rep(F, times=length(names(EFFECT2highlight))) 
              func.pch <- rep(dottype[1], times=length(names(EFFECT2highlight))) 
              func.col <- names(EFFECT2highlight)
              func.bg  <- names(EFFECT2highlight)
              func.lwd <- rep(2, times=length(names(EFFECT2highlight)))
           } else {func.legend <- NULL 
           func.lty <- NULL 
           func.pch <- NULL 
           func.col <- NULL 
           func.bg <- NULL 
           func.lwd <- NULL}
    
      } else {func.legend <- NULL 
    func.lty <- NULL 
    func.pch <- NULL 
    func.col <- NULL 
    func.bg <- NULL 
    func.lwd <- NULL}
  

  # legend for functional annotation and recombination rate is added if needed only
  legend("topleft", bty= "n", cex = cex.legend,    # ncol = 2, 
         title = "",
         legend = c(names(data), func.legend, rate.legend),
         lty = c(rep(F, length(data)), func.lty, rate.lty), 
         pch = c(dottype[names(data)], func.pch, rate.pch),
         lwd = c(rep(2, length(data)), func.lwd, rate.lwd),
         col = c(rep("black", length(data)), func.col, rate.col),
         bg = c(rep(par("bg"), length(data)), func.bg, rate.bg)
        )
  # add region coordinates as legend title. This done separetly because otherwise the legend itself
  # is centered below the (long) title and not left-adjusted.
  legend("topleft", legend="", title= paste0("chr", chr, ":", start, "-", end), bty= "n", cex = cex.legend)

 
  ##### Include recombination rate as lineplot
  if(!is.null(recombination.rate)) {
    
    # rate.max <- max(rate[,column.rate]) # scale recombination rate to maximum
    rate.max <- 100 # as in LocusZoom
    
    
    # cat("\n max rate: ", max(rate[,column.rate]))
     par(new = T)
     plot(c(start, end), c(0, rate.max),  type = "n", xlab = "", ylab = "", axes=F)   
     axis(side=4, cex.axis= cex.plot) # , padj=-1.5
     mtext(text= "Recombination Rate (cM/Mb)", side=4,  line=3, cex= cex.plot, col="black", lwd=2) #, padj=0.3

     points(rate[,column.pos], rate[,column.rate], type = "l", lwd=1.5, col= rate.col)
  } 


  ### plot genes 
  if (!is.null(biomaRt) || !is.null(LNCipedia)) {
    
    par(fig=c(0,1,0,0.4), new=TRUE)
    par(mar=c(1,5,1,5)) # number of margin lines: bottom, left, top, right   
    ystart =0 # start position for first gene layer
    ygenesize= 0.5 # height of each gene layer
 

     # process domain of selected gene if required 
    if(!is.null(plot.protein.domains)) {
      if(!any(names(plot.protein.domains) %in% gene.ranges$name)) {stop("Genes in plot.protein.domains not fond in region!\n")}
      
      gene.ranges$prot_dom <- gene.ranges$name %in% names(plot.protein.domains) 
      

      domains_list <- list() # initial list object for domain objects
    
    for (pd in names(plot.protein.domains)) {
        domains <- read.table(plot.protein.domains[pd], sep="\t", header=T, stringsAsFactors = F) # read domain data

      # check column names
      columns.requested <- c("domain_name_plot", "BP_start", "BP_end", "protein_length", "feature_length")
      if(!all(columns.requested %in% names(domains))) {
        stop(paste("columns missing in", plot.protein.domains[i], ": "), columns.requested[!columns.requested %in% names(domains)])}
      
      # Initialise additional columns if not given yet
      if (!("plotFeature" %in% names(domains))) {domains$plotFeature <- T}
      if (!("domain_height_extension" %in% names(domains))) {domains$domain_height_extension <- 1}
      if (!("domain_color" %in% names(domains))) {domains$domain_color <- "black"}
      if (!("symbol_plot" %in% names(domains))) {domains$symbol_plot <- "ellipse"}
      if (!("label_pos" %in% names(domains))) {domains$label_pos <- 3}
      if (!("assignArrows2Gene" %in% names(domains))) {domains$assignArrows2Gene <- T}
      
      
      # domain.min <- min(c(domains$BP_start, domains$BP_end)) 
      # domain.max <- max(c(domains$BP_start, domains$BP_end)) 
      # domain.totalLength <- domain.max - domain.min
      # domain.count <- nrow(domains)
      gene_start <- start(gene.ranges[pd])
      gene_end <- end(gene.ranges[pd])
      gene.totalLength <- width(gene.ranges[pd])
      AA.totalLength <- max(domains$protein_length)
      
      # remove features not to be plotted and order domains for start position
      domains <- domains[domains$plotFeature == T, ]
      domains <- domains[order(domains$BP_start, decreasing = F),]

            # ######### Previous domain plotting design with adjacent, non-overlapping symbols
            # # determine domain borders: #domains+1
            # domainBorders <- c(domain.min, domain.min + cumsum((domains$domain_width_extension * domain.totalLength) / sum(domains$domain_width_extension))) 
            # 
            # # determine start and end positions for each protein domain depending on strand of the gene
            # if (as.logical(strand(gene.ranges[gene.ranges$name== pd]) == "-")) { 
            #   domains$domstart <- domainBorders[-1]
            #   domains$domend <- domainBorders[-length(domainBorders)]
            #   
            # } else {
            #   domains$domstart <- domainBorders[-length(domainBorders)]
            #   domains$domend <- domainBorders[-1]
            # }
            #  ########### End previous plotting style with adjacent, non-overlapping symbols
      
      ## New domain plot design allows overlapping symbols
      # determine start and end positions for each protein domain depending on strand of the gene
      if (as.logical(strand(gene.ranges[gene.ranges$name== pd]) == "-")) { 
      
         domains$domstart <- gene_end - ((domains$start-1) * gene.totalLength) / AA.totalLength
         domains$domend   <- gene_end - ((domains$end-1) * gene.totalLength) / AA.totalLength
         } else {
           domains$domstart <- gene_start + ((domains$start-1) * gene.totalLength) / AA.totalLength
           domains$domend   <- gene_start + ((domains$end-1) * gene.totalLength) / AA.totalLength
         }
      ################### end of new plot design
      
      
      domains$xsymbol <- domains$domstart + 0.5*(domains$domend - domains$domstart) # x-mean for function symbols()
      domains$width <- abs(domains$domend - domains$domstart)   # for function symbols()
      domains$radius <- 0.5 * abs(domains$domend - domains$domstart) # for function symbols()
      domains$height <- ygenesize * domains$domain_height_extension / max(domains$domain_height_extension)  # consider different rectangle heigths 
      
 
     domains_list[[pd]] <- domains
        } # end pd loop
     } # end plot.protein.domains
    
       
    if(numberOfRowsForGenePlotting=="auto") {
      gene.ranges <- determineNonOverlapGenelayers(gene.ranges, ystart=ystart, ysize= ygenesize, minDistance=0.05,
                                             gene_name_column = "name", 
                                             prot_domain_column = if(!is.null(plot.protein.domains)) {"prot_dom"} else {NULL},
                                             scale.factor = gene.scale.factor, units="user", 
                                             sortGeneStart = is.null(plot.protein.domains)) # sort needs fixing
          levels <- NULL # will use newly created column 'ypos' as levels
          ymax <- max(gene.ranges$ypos)
          
        if(!is.null(plot.protein.domains)) {
          if(any(mcols(gene.ranges)[mcols(gene.ranges)$ypos == min(mcols(gene.ranges)$ypos), "prot_dom"] == TRUE)) {ystart <- ystart - ygenesize}
        }
          
      } else {ymax <- (numberOfRowsForGenePlotting / ygenesize) 
            levels=  numberOfRowsForGenePlotting  
            }
   
    

    

   ######## start plot gene panel
      plot(c(start, end), c(ystart-0.5*ygenesize, ymax), type = "n",  xlab = "", ylab = "", axes=F)   # axes
 
    ##### plot protein domains and arrows if requested 
    if (!is.null(plot.protein.domains)) {     
       
       for (pd in names(plot.protein.domains)) {  
              
            ystart_domain <-  mcols(gene.ranges[pd])$ypos  # gene.ranges[gene.ranges$name == pd, "ypos"]
            domains <- domains_list[[pd]]
            gene_start <- start(gene.ranges[pd])
            gene_end <- end(gene.ranges[pd])


            domainellipse <- domains[grepl("elli", domains$symbol_plot, ignore.case = T), ,drop=F]     
            if(nrow(domainellipse) >=1) {
              draw.ellipse(x = domainellipse$xsymbol, 
                         y = rep(ystart_domain - ygenesize, times=nrow(domainellipse)),
                         a = domainellipse$width * 0.5, # x-radius
                         b = domainellipse$height * 0.5, # y-radius
                         border = domainellipse$domain_color,
                         lwd = 1.3)
            }

            domainrect <- domains[grepl("rect", domains$symbol_plot, ignore.case = T), ,drop=F]     
            if(nrow(domainrect) >=1) {
            symbols(x= domainrect$xsymbol,
                    y= rep(ystart_domain - ygenesize, times=nrow(domainrect)),
                    rectangles= as.matrix(cbind(domainrect$width, domainrect$height)),
                    fg=domainrect$domain_color,
                    inches=FALSE,
                    add=TRUE,
                    lwd = 1.3)
            }
              
            domaincirc <- domains[grepl("circ", domains$symbol_plot, ignore.case = T), ,drop=F]     
            if(nrow(domaincirc) >=1) {
              symbols(x= domaincirc$xsymbol,
                    y= rep(ystart_domain - ygenesize, times=nrow(domaincirc)),
                    circles=domaincirc$radius,
                    fg=domaincirc$domain_color,
                    inches=FALSE,
                    add=TRUE,
                    lwd = 1.3)
            }
            
            # Continuous line for full protein length
            arrows(gene_start, ystart_domain - ygenesize, # arrow start: gene
                   gene_end, ystart_domain - ygenesize, # arrow tip: protein domain
                   angle = 90, code=3, lty=1, lwd=1.5, length=0.1, col="darkgrey")

            # domain labels
            text(domains$xsymbol, ystart_domain - ygenesize, labels= as.character(domains$domain_name_plot), 
                 cex = 0.7, pos= domains$label_pos, col= domains$domain_color) # adj=0.5
            
            ## plot arrows between gene and protein domain
            # Every domain is represented by 2 arrows indicating the corresponding genomic start and stop position.
            # Ideally, if genomic start and end positions of two neighboring domains are located closely
            # within a single exon, the two arrows belonging to this domain boundary are overlapping
            # and appear as single arrow for this boundary in the plot.
 
                  ## arrow positions for previous design with adjacent protein domains
                  # arrowstartx <- rbind(domains$BP_start, domains$BP_end) # n domains need 2*n arrows
                  # arrowstartx <- apply(arrowstartx, 2, sort) # sort start and stop coordinate per domain (necessary if minus strand)
                  # arrowstartx <- as.vector(arrowstartx)
                  # 
                  # domainBorders.within <- domainBorders[-c(1,length(domainBorders))] # start, middle domain borders twice (start,stop), end
                  # arrowtipsx <- c(domainBorders[1], as.vector(rbind(domainBorders.within, domainBorders.within)), domainBorders[length(domainBorders)])
                  #
                  # arrowtipsy <- rbind(domains$height[-length(domains$height)], domains$height[-1])
                  # arrowtipsy <- apply(arrowtipsy, 2, max) # height for arrow tip is maximum of two adjacent domain heights
                  # arrowtipsy <- c(domains$height[1], as.vector(rbind(arrowtipsy, arrowtipsy)), domains$height[length(domains$height)])
                  # arrowtipsy <- ystart_domain - ygenesize + 0.5*as.vector(arrowtipsy)
                  ## end arrow positions for previous design with adjacent protein domains
              
            domains.plotArrows <- domains[domains$assignArrows2Gene, ]
            arrowstartx <- c(domains.plotArrows$BP_start, domains.plotArrows$BP_end)
            arrowtipsx <- c(domains.plotArrows$domstart, domains.plotArrows$domend)
            
            #arrowtipsy <- c(domains.plotArrows$height, domains.plotArrows$height) # consider domain height for y-pos of arrow tips
            arrowtipsy <- 0.5 # arrow tips all at the same y-position
            arrowtipsy <- ystart_domain - ygenesize + 0.5*as.vector(arrowtipsy)
            
            arrows(arrowstartx, ystart_domain, # arrow start: gene
                   arrowtipsx, arrowtipsy, # arrow tip: protein domain
                   lty=5, lwd=1.2, length=0.08, col="darkgrey")
            } # end pd loop
      
   }   # end if plot.protein.domains


 
      # convert GRanges object to dataframe for plotting 
      genes <- granges2df(gene.ranges)

      # prepare dataframe with exon coordinates (must include gene coordinates for merging with gene dataframe!)
      #exons.df <- granges2df(exon.ranges)  
      exons <- granges2df(exon.ranges) # as.data.frame throws error in case of duplicated row names
      names(exons)[names(exons)=="seqnames"] <- "chrom"
      names(exons)[names(exons)=="start"] <- "exon_chrom_start"
      names(exons)[names(exons)=="end"] <- "exon_chrom_end"
      names(exons)[names(exons)=="startGene"] <- "start"
      names(exons)[names(exons)=="endGene"] <- "end"
      
    # call plot function
    regionalplot.genelabels(
      genes = genes,
      levels = levels,
      exons = exons,
      xstart = start,
      xstop = end,
      ystart = 0,
      ygenesize= ygenesize,
      scale.factor = gene.scale.factor,
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
  
 
      
  } # end if(is.null(biomaRt|LNCipedia))
  
  # re-set graphical parameter
  par(fig = par.old.fig)
  par(mar = par.old.mar)
  par(oma = par.old.oma)
  

# Detaching libraries not needed any more
detach_package(unique(pks2detach))

} # end of function definition



