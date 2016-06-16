

## Description
# Get TF Binding sites from input genes

## Usage
TFsearch <- function(sequences, 
                     newheader = NULL, # c("chr", "start", "stop", "customname", "custom2", "strand"), 
                     annoColumn = NULL, # e.g. "customname" 
                     name.organism="hsapiens", 
                     projectfolder=pipepar[["outdir"]],
                     projectname="",
                     applyFilter = FALSE,
                        filtercat1 = "adj.P.Val",
                        filtercat1.decreasing = FALSE,
                        filtercat1.function = abs,
                        filtercat1.threshold=pipepar[["threshold_p"]],
                        filtercat2 = "logFC",
                        filtercat2.decreasing = TRUE,
                        filtercat2.function = abs,
                        filtercat2.threshold= log2(pipepar[["ThresholdFC"]]),
                     PromLookup = TRUE,
                       Entrez.col = "ENTREZID", 
                       PromSeqUpstreamTSS = 2000,
                       PromSeqDownstreamTSS = 200,
                     SearchSelMotifs = NULL,
                     motif.min.score = 0.9 
                    ) {

  
  ## Arguments
  # sequences: dataframe (or character with file path) with column containing Entrez IDs or sequence coordinates.
  #            Required columns as defined in 'PromLookup' (see below)
  # newheader: NULL if 'sequences' already supplied with header. Character vector with new header otherwise.
  # annoColumn (character or vector of characters): column name(s) of 'sequences'-object with sequence annotation to maintain
  # name.organism: currently human data only (hg19)
  # projectfolder: output directory
  # projectname: character merged to output name
  # applyFilter (boolean): If TRUE, sequences are filtered for applied categories and thresholds. Filter Values converted to ABSOLUTE values.
  #     Optional Filtering criteria (Ignored if applyFilter==FALSE):
  #     filtercat1: column name of first category to filter 'sequences' (e.g. p-values). 
  #     filtercat1.decreasing (boolean): direction to order filtercat1
  #     filtercat1.function: select transforming function for filter category1 (no quotes). e.g. abs for absolute values, identity for no transformation
  #     filtercat1.threshold: Threshold for filtercat1 or 'top123' for top Hits
  #     filtercat2: column name of second category to filter 'sequences' (e.g. effect size).
  #     filtercat2.decreasing (boolean): direction to order filtercat2
  #     filtercat2.function: select transforming function for filter category2 (no quotes). E.g. abs for foldchanges
  #     filtercat2.threshold: Threshold for filtercat2 or 'top123' for top Hits
  # PromLookup (boolean): if TRUE, all promotor sequences corresponding to genes in 'sequences' are downloaded. Therefore a
  #                         column with EntrezIDs is requried (column name given in 'Entrez.col')!
  #                       if FALSE, Sequences are downloaded according to given coordinates in 'sequences'. Therefore columns 
  #                         for chromosome, start, stop and strand information required! Additional meta columns allowed.
  #   Entrez.col: name of column containing Entrez IDs
  #   PromSeqUpstreamTSS, PromSeqDownstreamTSS: definition of promotor regions to download respective to TSS
  # SearchSelMotivs: Character Vector of selected motives to search in 'sequences' seperately to enrichment or NULL
  # motif.min.score: minimum score to match motif pwm to target sequence (ignored if SearchSelMotifs = NULL)
  
  
  ## Details
  # Function takes input sequences, optionally looks up all promotor sequences by (unique) entrezID, 
  # or looks up sequences defined by coordinates, optionally filters 'sequences' for given 
  # categories and calls motifEnrichment() to identify enriched Transcription factor binding motivs. 
  # Optionally tables generated for in 'SearchSelMotifs' preselected Motivs looked up in 'sequences'.
  
  
  ## Value
  # groupReport of motifEnrichment results.
  
  
  ## Author(s) 
  # Frank Rühle 
  
  
  
  # create output directory
  if (!file.exists(file.path(projectfolder, "TFBS"))) {dir.create(file.path(projectfolder, "TFBS"), recursive=T)}
  
  # check organism
  if (!grepl("human|sapiens", name.organism, ignore.case = TRUE)) {stop("TFsearch for human data only")} 
  
  # load required libraries
  scriptpackages <- c("MotifDb", "GenomicFeatures", "PWMEnrich", "PWMEnrich.Hsapiens.background", "seqLogo", "AnnotationHub", "org.Hs.eg.db", 
                      "Homo.sapiens", "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19")
  for (lib in scriptpackages) {library(lib, character.only = TRUE)}
  data(PWMLogn.hg19.MotifDb.Hsap)
  genome <- Hsapiens
  
        
  
  # read file if 'sequences' is character string with file path
  if(is.character(sequences)) {
    cat("\n\nReading:", sequences, "\n")
    sequences <- read.table(sequences, header=is.null(newheader), sep="\t", na.strings = c("", " ", "NA")) # if no newheader, header must be in file
    }
  if(!is.null(newheader)){# if 'newheader' is defined, it is used as names(sequences)
    cat("\nNew header added to input file:", newheader, "\n")
    names(sequences) <- newheader
    }  
    
  
  # Optional Filtering of 'sequences' for supplied categories
  if(applyFilter) {
    sequences <- filterGeneLists(sequences,
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
  
  
  
  ## Lookup of promotor sequences (optional)
  if(PromLookup) { # if TRUE, lookup promotor sequences
     if(is.null(Entrez.col)) {stop("EntrezIDs required for lookup of promotor sequences!")}     
       
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
        seqEG <- as.character(unique(sequences[!is.na(sequences[,Entrez.col]),Entrez.col]))
                          
        # Obtain the coordinates of the transcripts in a gene ranges object.
        # (extract a GRangesList object that groups the transcripts by the genes that they are associated with)
        grl <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")
        grl <- keepStandardChromosomes(grl)
        
        cat("\n\n", length(setdiff(seqEG,names(grl))), "Entrez IDs from analysis not found in Transcription database:\n")
        print(setdiff(seqEG,names(grl)))
        
        # filter for ENTREZ IDs available in grl, other wise error when subsetting grl with entrez IDs
        seqEG <- seqEG[seqEG %in% names(grl)] 
        grl <- grl[seqEG]
        
        
        # download promotor sequences
        cat("\nGet Promotor Sequences from", PromSeqUpstreamTSS, "bp upstream and", PromSeqDownstreamTSS, "downstream from TSS.\n")
        seqs2search <- getPromoterSeq(grl, genome, upstream=PromSeqUpstreamTSS, downstream=PromSeqDownstreamTSS)
        seqs2search <- unlist(seqs2search) # need DNAStringSet instead of DNAStringSetList
        names(seqs2search) <- sub("\\..*$", "", names(seqs2search))  # solve names of form: GeneID.GeneID      
        cat("\nFound", length(seqs2search), "Promotor sequences for", length(seqEG), "Entrez IDs.\n")
        
        # edit meta colums still present in unlisted DNAStringSet        
        SymbolAnnotation <- select(org.Hs.eg.db, keys=names(seqs2search), columns=c("SYMBOL"), keytype="ENTREZID")
        mcols(seqs2search)$SYMBOL.orgdb <- SymbolAnnotation$SYMBOL # annotate Entrez IDs with Symbols by org.Hs.eg.db
        if(!is.null(annoColumn)) { # annotate meta data by columns given in annoColumn (if any)
          mcols(seqs2search) <- DataFrame(plyr::join(as.data.frame(mcols(seqs2search)), 
                             data.frame(geneID= sequences[,Entrez.col], sequences[,annoColumn, drop=F]), by="geneID", type="left", match="first"))
         
                  
        }



  } else { # if PromLookup==FALSE, sequences are downloaded according to given start and stop coordiates in 'sequences'
  
    # create gene ranges from input 'sequences' 
      gr <-  makeGRangesFromDataFrame(sequences,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo = NULL,
                                    seqnames.field = grep("chr", names(sequences), value=T, ignore.case=T),
                                    start.field = grep("(begin)|(start)", names(sequences), value=T, ignore.case=T),
                                    end.field = grep("(end)|(stop)", names(sequences), value=T, ignore.case=T),
                                    strand.field = grep("strand", names(sequences), value=T, ignore.case=T),
                                    starts.in.df.are.0based=FALSE)
    
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
res = motifEnrichment(seqs2search, pwms=PWMLogn.hg19.MotifDb.Hsap, 
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

report = groupReport(res, top=0.05)  # top: what proportion of top motifs should be examined in each individual sequence (by default 5%)
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
# report = groupReport(res, by.top.motifs=TRUE)
  
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






cat("\nWriting PWMEnrich report to:\n", file.path(projectfolder, "TFBS", paste0(projectname, "_PWMEnrich_report.txt")), "\n")
reportdf <- as.data.frame(report)
# names(reportdf)[names(reportdf)=="p.value"] <- "avg.log2.p" 
# reportdf$p.value <- 2^reportdf$avg.log2.p  # das macht keinen Sinn! 
write.table(reportdf, file=file.path(projectfolder, "TFBS", paste0(projectname, "_PWMEnrich_report.txt")), sep="\t", quote=F, row.names=F)

pdf(file.path(projectfolder, "TFBS", paste0(projectname, "_PWMEnrich_top.pdf")), paper="a4") 
plot(report[1:15], fontsize=7, id.fontsize=5)
dev.off()






######## Optional Search for preselected motivs given in 'SearchSelMotivs'
if(!is.null(SearchSelMotifs)) {
  if (!file.exists(file.path(projectfolder, "TFBS", "Search_Selected_Motifs"))) {dir.create(file.path(projectfolder, "TFBS", "Search_Selected_Motifs"))}
  cat("\nSearch for selected motivs in input sequences", SearchSelMotifs)
  
  # check if motifs of SearchSelMotifs exist in databse
  motifexists <- sapply(SearchSelMotifs, function(mot) length(MotifDb::query(MotifDb::query(MotifDb, name.organism), mot)) > 0)
  if (!any(motifexists)) {stop("no motifs found in database")}  
  cat("\n", sum(motifexists), "of", length(motifexists), "motifs found in database.")
  if (!all(motifexists)) {cat("\nNot found:", SearchSelMotifs[!motifexists], "\n")}
    
  # non-existing motifs removed  
  SearchSelMotifs <- SearchSelMotifs[motifexists]
  
  cat("\nPlotting sequence logos and writing count tables to:\n",file.path(projectfolder, "TFBS", "Search_Selected_Motifs"), "\n")
  
  # initialise parameter for motif loop
    listpfm <- list()
    listpcm <- list()
    listhits <- list()
    countmotif <- data.frame(matrix(nrow=length(seqs2search), ncol=0))
    
  for(motif in SearchSelMotifs) { # find hits in sequences for every motiv
    # create instances of the pfm class for every motiv
    listpfm[[motif]] <- new("pfm", mat=MotifDb::query (MotifDb::query (MotifDb, name.organism), motif)[[1]], name=motif)
    listpcm[[motif]] <-  round(100 * listpfm[[motif]]@mat)  
  
      # plot sequence logo
      tiff(file.path(projectfolder, "TFBS", "Search_Selected_Motifs", paste0(motif, "_SequenceLogo.tiff")), 
           width = 7016 , height = 4960, res=600, compression = "lzw")
      seqLogo(listpfm[[motif]]@mat)
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
              file=file.path(projectfolder, "TFBS", "Search_Selected_Motifs", paste0("Motifcount_", projectname, ".txt")))
  
  } # end of if(!is.null(SearchSelMotifs))
  
  

# Detaching libraries not needed any more
for (lib in scriptpackages) {detach_package(lib)}

return(report)

} # end of TFsearch




