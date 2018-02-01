

#' Transcription factor binding site enrichment
#' 
#' This function uses the \code{PWMEnrich}-package for TFBS enrichment.
#' 
#' The function takes input genes in \code{sequences} and looks up all promotor sequences by (unique) entrezIDs refering to human genome build hg19.
#' If \code{PromLookup == FALSE}, look up of promotor sequences is omitted and sequences of interest must be given in \code{sequences}
#' as coordinates instead. The dataset may be filtered for designated filter criteria if desired. 
#' These sequences are transferred to \code{motifEnrichment} from \code{PWMEnrich}-package to identify enriched Transcription factor 
#' binding motivs for the input sequences. Optionally, preselected motivs given in \code{SearchSelMotifs} are looked up in \code{sequences}. 
#' All result tables and plots are stored in the project folder.
#' 
#' @param sequences dataframe or character with file path to dataframe or named list containing dataframes. The dataframe must either contain
#' a column with Entrez IDs or sequence coordinates if \code{PromLookup == FALSE}. Required columns for latter case are defined 
#' in \code{PromLookup} below.
#' @param newheader NULL if \code{sequences} already supplied with header. Character vector with new header otherwise.
#' @param annoColumn character or vector of characters. Column name(s) of \code{sequences}-object with sequence annotation to maintain.
#' @param name.organism currently human data only (hg19).
#' @param projectfolder output directory.
#' @param projectname character prefix for output name.
#' @param figure.res numeric resolution for png.
#' @param applyFilter (boolean) If TRUE, sequences are filtered for applied categories and thresholds. Filter Values converted to ABSOLUTE values.
#'     Optional Filtering criteria (Ignored if \code{applyFilter=FALSE}):
#' @param filtercat1 column name of first category to filter \code{sequences} (e.g. p-values). 
#' @param filtercat1.decreasing (boolean) direction to order \code{filtercat1}.
#' @param filtercat1.function select transforming function for filter category1 (no quotes). e.g. \code{abs} for absolute values, \code{identity} for no transformation
#' @param filtercat1.threshold Threshold for \code{filtercat1} or 'top123' for top Hits
#' @param filtercat2 column name of second category to filter \code{sequences} (e.g. effect size).
#' @param filtercat2.decreasing (boolean) direction to order \code{filtercat2}.
#' @param filtercat2.function select transforming function for filter category2 (no quotes). E.g. \code{abs} for foldchanges
#' @param filtercat2.threshold Threshold for \code{filtercat2} or 'top123' for top Hits
#' @param PromLookup (boolean) if TRUE, all promotor sequences corresponding to genes in \code{sequences} are downloaded. 
#' Therefore a column with EntrezIDs is requried (column name given in \code{Entrez.col}). if FALSE, Sequences are 
#' downloaded according to given coordinates in \code{sequences}. Therefore columns for \code{chromosome}, \code{start}, 
#' \code{stop} and \code{strand} information are required! Additional meta columns allowed.
#' @param id.type character with identifier type from annotation package ("ENTREZID" or "SYMBOL")
#'          Gene symbols Will be converted to EntrezIDs prior to enrichment analysis.
#' @param id.column character with column name for identifier variable in \code{sequences}. 
#' @param PromSeqUpstreamTSS definition of promotor regions to download upstream to TSS.
#' @param PromSeqDownstreamTSS: definition of promotor regions to download downstream to TSS.
#' @param SearchSelMotivs Character Vector of selected motives to search in \code{sequences}. Omitted if NULL.
#' @param motif.min.score minimum score to match motif pwm to target sequence (ignored if SearchSelMotifs = NULL).
#' 
#' @return groupReport of motifEnrichment results.
  
#' @seealso MotifDb, PWMEnrich, PWMEnrich.Hsapiens.background
#' 
#' @author Frank Ruehle
#' 
#' @export 
#' 


TFsearch <- function(sequences, 
                     newheader = NULL,  
                     annoColumn = NULL,  
                     name.organism="hsapiens", 
                     projectfolder= "GEX/TFBS",
                     projectname="",
                     figure.res = 300,
                     applyFilter = FALSE,
                     filtercat1 = "adj.P.Val",
                     filtercat1.decreasing = FALSE,
                     filtercat1.function = abs,
                     filtercat1.threshold= 0.05,
                     filtercat2 = "logFC",
                     filtercat2.decreasing = TRUE,
                     filtercat2.function = abs,
                     filtercat2.threshold= log2(1.5),
                     PromLookup = TRUE,
                     id.type = "ENTREZID",
                     id.column = "ENTREZID",   
                     PromSeqUpstreamTSS = 2000,
                     PromSeqDownstreamTSS = 200,
                     SearchSelMotifs = NULL,
                     motif.min.score = 0.9) 
{
  
  
  
  # create output directory
  if (!file.exists(file.path(projectfolder))) {dir.create(file.path(projectfolder), recursive=T)}
  
  if(!is.null(projectname)) { # add underline to projectname if given
    if(projectname!="") {
      projectname <- paste0(projectname, "_")
    }}
  
  
  # check organism
  if (!grepl("human|sapiens", name.organism, ignore.case = TRUE)) {stop("TFsearch for human data only")} 
  
  # load required libraries
  pkg.cran <- NULL
  pkg.bioc <- c("MotifDb", "motifStack", "GenomicFeatures", "PWMEnrich", "PWMEnrich.Hsapiens.background", "seqLogo", "AnnotationHub", 
                      "Homo.sapiens", "BSgenome.Hsapiens.UCSC.hg19") # "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene",
  pks2detach <- attach_package(pkg.cran=pkg.cran, pkg.bioc=pkg.bioc)    
  
  data(PWMLogn.hg19.MotifDb.Hsap)
  genome <- Hsapiens
  
   
  if(class(sequences) != "list") { # make list if sequences is single element
    sequences <- list(sequences=sequences)
  }
  
  # initialise result list
  report <- list()
  
  for(sq in names(sequences)) {
    
    cat(paste("\nProcessing", sq, "\n"))
  
  
  # read file if 'sequences[[sq]]' is character string with file path
  if(is.character(sequences[[sq]]) && length(sequences[[sq]])==1) {
    cat("\n\nReading:", sequences[[sq]], "\n")
    sequences[[sq]] <- read.table(sequences[[sq]], header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
    }
  if(!is.null(newheader)){# if 'newheader' is defined, it is used as names(sequences[[sq]]). e.g. c("chr", "start", "stop", "customname", "custom2", "strand"),
    cat("\nNew header added to input file:", newheader, "\n")
    names(sequences[[sq]]) <- newheader
    }  
    
  if(is.vector(sequences[[sq]], mode="character")) {
    sequences[[sq]] <- dataframe(sequences[[sq]])
    names(sequences[[sq]]) <- id.column
    }
   
    
      
     
  # Optional Filtering of 'sequences[[sq]]' for supplied categories
  if(applyFilter) {
    sequences[[sq]] <- filterGeneLists(sequences[[sq]],
                                newheader=newheader,
                                filtercat1 = filtercat1,
                                filtercat1.decreasing = filtercat1.decreasing,
                                filtercat1.function = filtercat1.function,
                                filtercat1.threshold=filtercat1.threshold,
                                filtercat2 = filtercat2,
                                filtercat2.decreasing = filtercat2.decreasing,
                                filtercat2.function = filtercat2.function,
                                filtercat2.threshold= filtercat2.threshold)
    } # end of if(applyFilter)
  
  
  
  ## Lookup of promotor sequences 
  if(PromLookup) { # if TRUE, lookup promotor sequences
     if(!(id.type %in% c("ENTREZID", "SYMBOL"))) {stop("Either ENTREZID or SYMBOL required for lookup of promotor sequences!")}     
 
    ## convert IDs to ENTREZ IDs if necessary
    if(id.type!="ENTREZID") {

      entrezids <- basicAnno(data=sequences[[sq]], Symbol.column = id.column, Entrez.column = NULL, org=name.organism)
      
      cat(paste(nrow(entrezids), " of ", length(unique(sequences[[sq]][,id.column])), "unique", id.type, "mapped to ENTREZIDs\n"))
      cat(paste("For dublicated", id.type, "only the first element is used)\n")) 
      sequences[[sq]] <- plyr::join(sequences[[sq]], entrezids, by=id.column, type="right", match="first") 
      # joined dataframe contains only entries with available ENTREZIDs, but keeps order of 'genes[[ge]]'
    } else {names(sequences[[sq]])[names(sequences[[sq]])==id.column] <- "ENTREZID"}
    
    
    
          
#      ### accessor functions for TranscriptDb objects:
#      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#      columns(txdb)
#      keytypes(txdb) # "GENEID"   "TXID"     "TXNAME"   "EXONID"   "EXONNAME" "CDSID"    "CDSNAME"
#      select(txdb, keys = c("100033416", "100033417"), columns=c("TXID", "TXNAME", "EXONID", "EXONNAME", "CDSID", "CDSNAME"), keytype="GENEID")
#      head(transcripts(txdb)) # get all transcripts
#      head(transcripts(txdb, vals <- list(tx_chrom = "chr15", tx_strand = "+"))) # subsetting transcripts
#      head(exons(txdb)) # get all exons
#      head(cds(txdb)) # get all coding sequences
     
     # Vector with Entrez IDs (only unique ENTREZ IDs)

        seqEG <- as.character(unique(sequences[[sq]][!is.na(sequences[[sq]][,"ENTREZID"]),"ENTREZID"]))

        # Obtain the coordinates of the transcripts in a gene ranges object.
        # (extract a GRangesList object that groups the transcripts by the genes that they are associated with)
        grl <- GenomicFeatures::transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")
        grl <- GenomeInfoDb::keepStandardChromosomes(grl, pruning.mode= "coarse")

        cat("\n\n", length(setdiff(seqEG,names(grl))), "Entrez IDs from analysis not found in Transcription database:\n")
        print(setdiff(seqEG,names(grl)))
        
        # filter for ENTREZ IDs available in grl, other wise error when subsetting grl with entrez IDs
        seqEG <- seqEG[seqEG %in% names(grl)] 
        grl <- grl[seqEG]
        
        
        # download promotor sequences
        cat("\nGet Promotor Sequences from", PromSeqUpstreamTSS, "bp upstream and", PromSeqDownstreamTSS, "downstream from TSS.\n")
        seqs2search <- GenomicFeatures::getPromoterSeq(grl, genome, upstream=PromSeqUpstreamTSS, downstream=PromSeqDownstreamTSS)
        seqs2search <- unlist(seqs2search) # need DNAStringSet instead of DNAStringSetList
        names(seqs2search) <- sub("\\..*$", "", names(seqs2search))  # solve names of form: GeneID.GeneID      
        cat("\nFound", length(seqs2search), "Promotor sequences for", length(seqEG), "Entrez IDs.\n")
        
        # edit meta colums still present in unlisted DNAStringSet        
        SymbolAnnotation <- select(org.Hs.eg.db, keys=names(seqs2search), columns=c("SYMBOL"), keytype="ENTREZID")
        mcols(seqs2search)$SYMBOL.orgdb <- SymbolAnnotation$SYMBOL # annotate Entrez IDs with Symbols by org.Hs.eg.db
        if(!is.null(annoColumn)) { # annotate meta data by columns given in annoColumn (if any)
          mcols(seqs2search) <- DataFrame(plyr::join(as.data.frame(mcols(seqs2search)), 
                             data.frame(geneID= sequences[[sq]][,"ENTREZID"], sequences[[sq]][,annoColumn, drop=F]), by="geneID", type="left", match="first"))
        }



  } else { # if PromLookup==FALSE, sequences are downloaded according to given start and stop coordiates in 'sequences[[sq]]'
  
    # create gene ranges from input 'sequences' 
    if(is.data.frame(sequences[[sq]])) { 
      gr <-  GenomicRanges::makeGRangesFromDataFrame(sequences[[sq]],
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo = NULL,
                                      seqnames.field = grep("chr", names(sequences[[sq]]), value=T, ignore.case=T),
                                      start.field = grep("(begin)|(start)", names(sequences[[sq]]), value=T, ignore.case=T),
                                      end.field = grep("(end)|(stop)", names(sequences[[sq]]), value=T, ignore.case=T),
                                      strand.field = grep("strand", names(sequences[[sq]]), value=T, ignore.case=T),
                                      starts.in.df.are.0based=FALSE)
    }
    
    # remove duplicates if necessary
    if(any(duplicated(gr))) {
      cat("\nRemove", sum(duplicated(gr)),"duplicated entries")
      gr <- gr[!duplicated(gr)]
      }
    
    # get sequence information from gene ranges
    seqs2search <- getSeq(genome, gr)
    
    # annotate meta data by columns given in annoColumn (if any)
    if(!is.null(annoColumn)) { 
    mcols(seqs2search) <- DataFrame(matrix(nrow=length(seqs2search), ncol=0))
    mcols(seqs2search)[annoColumn] <- mcols(gr)[annoColumn]
    }
  }

  
######## motif Enrichment using the lognormal affnity background distribution
# The main function is motifEnrichment() which took our sequences and calculated motif enrichment
# using the lognormal affnity background distribution. We then applied the groupReport function 
# to calculate the enrichment over the whole group of sequences. This produced a ranked list of 
# motifs according to the estimated P-values.
# motifEnrichment() needs either a single sequence (an object of class DNAString), or a list of DNAString objects, 
# or a DNAStringSet object.

cat("\n\nStarting Motif Enrichment with lognormal distribution background of PWMs pre-computed in 'PWMLogn.hg19.MotifDb.Hsap'\n")
res = PWMEnrich::motifEnrichment(seqs2search, pwms=PWMLogn.hg19.MotifDb.Hsap, 
                      score = "autodetect",
                      bg = "autodetect") # bg = "PWMEnrich.Hsapiens.background" and bg ="autodetect" gives same results when pwms=PWMLogn.hg19.MotifDb.Hsap
# ### a MotifEnrichmentResults object containing a subset following elements:
# "score" - scoring scheme used
# "bg" - background correction used
# "params" - any additional parameters
# "sequences" - the set of sequences used
# "pwms" - the set of pwms used
# "sequence.nobg" - per-sequence scores without any background correction. 
# "sequence.bg" - per-sequence scores after background correction. For "logn" and "pval" the P-value (smaller is better).
# "group.nobg" - aggregate scores for the whole group of sequences without background correction. 
#               --> 'raw.score' in groupReport
# "group.bg" - aggregate scores for the whole group of sequences with background correction. 
#              For "logn" and "pval", the P-value for the whole group (smaller is better).
#              --> 'p.value' in groupReport
# "sequence.norm" - (only for "logn") the length-normalized scores for each of the sequences. 
# "group.norm" - (only for "logn") similar to sequence.norm, but for the whole group of sequences.

### warning generated if 'N's are present in sequence:
# In .Call2("PWM_score_starting_at", pwm, subject, starting.at,  ... : 'subject' contains letters not in [ACGT] ==> assigned weight 0 to them

report[[sq]] <- PWMEnrich::groupReport(res, top=0.05)  # top: what proportion of top motifs should be examined in each individual sequence (by default 5%)
# ### a MotifEnrichmentReport object containing a table with the following columns:
# 'rank' - The rank of the PWM's enrichment in the whole group of sequences together
# 'target' - The name of the PWM's target gene, transcript or protein complex.
# 'id' - The unique identifier of the PWM (if set during PWM creation).
# 'raw.score' - The raw score before P-value calculation
# 'p.value' - The P-value of motif enrichment (if available)
# 'top.motif.prop' - The proportion (between 0 and 1) of sequences where the motif is within top proportion of enrichment motifs.

# The last column indicates the breadth of enrichment using a 5% ranking threshold. This column helps to differentiate cases 
# where the motif enrichment is strongly focused to a small subset of sequences (in which case breadth is small), 
# versus being more widespread but weaker (in which case breadth is bigger). We can also sort by this column:
# report[[sq]] = groupReport(res, by.top.motifs=TRUE)
  
### Special conditions for Human Sequences:
# the only major difference is that the P-value of groups of sequences (i.e. groupReport()) is replaced with an 
# average of log(P-values) of individual sequences. To get most complete results we recommend examining motif enrichment 
# based on this score and breadth (last column of groupReport)

### Comments from Robert Stojnic (package author) in mails from 23.05.2016 and 25.05.2016
# Raw score is an average of affinity scores, i.e. exp(PWM score) over the whole sequence. The numerical values 
# is itself not very informative unless it's compared to the distribution over the genome. 
# The P-value is the P-value of over-representation. Some motifs will have P-values very close to 1, 
# so they are effectively under-represented. However, the statistics is designed to model the over-representation tail, 
# so these under-representation P-values are not very accurate. 
# If you have the latest version of the page, what you will be seeing for human is the P-value, although the 
# calculation is done on log(P-values) because it's more numerically stable.






cat("\nWriting PWMEnrich report to:\n", file.path(projectfolder, paste0(projectname, sq, "_PWMEnrich_report.txt")), "\n")
reportdf <- as.data.frame(report[[sq]])
write.table(reportdf, file=file.path(projectfolder, paste0(projectname, sq, "_PWMEnrich_report.txt")), sep="\t", quote=F, row.names=F)

#pdf(file.path(projectfolder, paste0(projectname, sq, "_PWMEnrich_top.pdf")), paper="a4") 
png(file.path(projectfolder, paste0(projectname, sq, "_PWMEnrich_top.png")), width = 210, height = 297, units = "mm", res=figure.res) 
PWMEnrich::plot(report[[sq]][1:10], fontsize=10, id.fontsize=10)
dev.off()


  

######## Optional Search for preselected motivs given in 'SearchSelMotivs'
if(!is.null(SearchSelMotifs)) {
  if (!file.exists(file.path(projectfolder, "Search_Selected_Motifs"))) {dir.create(file.path(projectfolder, "Search_Selected_Motifs"))}
  cat("\nSearch for selected motivs in input sequences", SearchSelMotifs)
  
  # check if motifs of SearchSelMotifs exist in databse
  motifexists <- sapply(SearchSelMotifs, function(mot) length(MotifDb::query(MotifDb::query(MotifDb, name.organism), mot)) > 0)
  if (!any(motifexists)) {stop("no motifs found in database")}  
  cat("\n", sum(motifexists), "of", length(motifexists), "motifs found in database.")
  if (!all(motifexists)) {cat("\nNot found:", SearchSelMotifs[!motifexists], "\n")}
    
  # non-existing motifs removed  
  SearchSelMotifs <- SearchSelMotifs[motifexists]
  
  cat("\nPlotting sequence logos and writing count tables to:\n",file.path(projectfolder, "Search_Selected_Motifs"), "\n")
  
  # initialise parameter for motif loop
    listpfm <- list() # position frequency matrix (PFM)
    listpcm <- list()
    listhits <- list()
    countmotif <- data.frame(matrix(nrow=length(seqs2search), ncol=0))
    
  for(motif in SearchSelMotifs) { # find hits in sequences for every motiv
    # create instances of the pfm class for every motiv (object classes from motifStack package)
    listpfm[[motif]] <- new("pfm", mat=MotifDb::query (MotifDb::query (MotifDb, name.organism), motif)[[1]], name=motif) # selects first matrix entry
    listpcm[[motif]] <-  round(100 * listpfm[[motif]]@mat)  
  
      # plot sequence logo
      #tiff(file.path(projectfolder, "Search_Selected_Motifs", paste0(motif, "_SequenceLogo.tiff")), width = 7016 , height = 4960, res=600, compression = "lzw")
      png(file.path(projectfolder, "Search_Selected_Motifs", paste0(projectname, sq, "_", motif, "_SequenceLogo.png")), width = 140, height = 100, units = "mm", res=figure.res) 
      seqLogo::seqLogo(listpfm[[motif]]@mat)
      dev.off()  
      
      # All of the matches in 'seqs2search' to each binding motif are found at once with this command:
      listhits[[motif]] <- sapply(seqs2search, function(pseq) matchPWM(listpcm[[motif]], pseq, min.score=motif.min.score))
      # Now count their lengths:
      countmotif[motif] <- sapply(listhits[[motif]], length)
  }
  
  
  # Create a data.frame from this information:
  entrezIDs <- names(seqs2search)
  # check if meta information exists
  metacolumns <- c("SYMBOL.orgdb", "tx_name", annoColumn)
  metacolumnsexist <- metacolumns %in% names(mcols(seqs2search))
  metacolumns <- metacolumns[metacolumnsexist]
  
  tbl.gata  <- data.frame(ENTREZ=entrezIDs, mcols(seqs2search)[metacolumns], countmotif[,,drop=F])
  write.table(tbl.gata, row.names=F, quote=F, sep="\t", 
              file=file.path(projectfolder, "Search_Selected_Motifs", paste0("Motifcount_", projectname, sq, ".txt")))
  
  } # end of if(!is.null(SearchSelMotifs))
  
} # end sq loop

# Detaching libraries not needed any more
detach_package(unique(pks2detach))


return(report[[sq]])

} # end of TFsearch




