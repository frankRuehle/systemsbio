#' Calculate genomic coordinates for amino acid positions of protein domains 
#' 
#' Use genomic exon coordinates and Uniprot protein domain data to determine genomic 
#' coordinates of protein domains for plotting in genomic context.
#' 
#' Genomic exon coordinates for the desired transcript are downloaded from biomaRt. If a gene is selected,
#' the canonical transcript is determined as transcript with longest coding sequence.
#' Protein data is expected in \code{gff}-format from Uniprot giving the amino acid positions for each domain. 
#' The total protein length is read from the 2nd comment line of the \code{gff}-file.
#' The exon bp coordinates from the selected transcript are used to calculate the cooresponding bp coordinates
#' for each protein domain based on their amino acid positions.
#' The strored result table may be edited later to modify suggested plotting parameter. The name of the 
#' used gene transcript is added to the output filename.
#' 
#' @param ID character with either Ensembl transcript ID, Ensembl gene ID or gene symbol.
#' @param biomaRt biomaRt object to obtain exon ccordinates.
#' @param uniprot_domains.gff character with file path to Uniprot gff export. 
#' In Uniprot, select desired features e.g. \code{PTM/Processing} and \code{Family & Domains} 
#' and export the basket to gff-file format.
#' @param suffix.outputFilname character to be added at output filename additionally to the used transcript name.  
#' 
#' @return dataframe with exon data downloaded from \code{biomaRt}. 
#' The result domain table is stored as side effect in the filepath given by \code{uniprot_domains.gff}.
#' 
#' @author Frank Ruehle
#' 
#' @export 



makeDomainsFromExons <- function(ID, biomaRt, uniprot_domains.gff, suffix.outputFilname =".txt") {


  # load required libraries
  pcksBioc <- c("Biobase", "biomaRt")
  pks2detach <- attach_package(pcksBioc)
  
  
  
     if(grepl("ENSG", ID)) {filter.ensembl <- "ensembl_gene_id"} else {
         if(grepl("ENST", ID)) {filter.ensembl <- "ensembl_transcript_id"} else {
                                     filter.ensembl <-  "hgnc_symbol"}} # or "external_gene_id" in hg38
  

  attr.exons <- c("ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end",  
                 "rank", "genomic_coding_start", "genomic_coding_end", 
                 "cds_start", "cds_end",  # exon coordinates
                 "ensembl_gene_id", "ensembl_transcript_id", "external_gene_id", "strand")  
  
  geneexons <- getBM(attributes = attr.exons, 
                 filters = filter.ensembl,
                 values = ID, 
                 mart = biomaRt)
  
  

  # choose canonical transcript (if no transcript ID given)
  # canonical transcript: For human, the canonical transcript for a gene is set according to the following hierarchy:
  # 1. Longest CCDS translation with no stop codons. 
  # 2. If no (1), choose the longest Ensembl/Havana merged translation with no stop codons. 
  # 3. If no (2), choose the longest translation with no stop codons. 
  # 4. If no translation, choose the longest non-protein-coding transcript.
  
  # remove entries on chromosome patches
  geneexons <- geneexons[!grepl("(chr)|(patch)", geneexons$chromosome_name, ignore.case = T),]
  # remove entries without coding sequence
  geneexons <- geneexons[!is.na(geneexons$cds_start), ]
  
  # choose canonical transcript (longest coding region)
  transcript <- tapply(geneexons$cds_end, geneexons$ensembl_transcript_id, FUN = max) # named character vector with transcript length
  transcript <- names(transcript[transcript==max(transcript)])[1] # name of longest transcript
  
  geneexons <- geneexons[geneexons$ensembl_transcript_id == transcript, ] # select exons from selected transcript
  
  gene_strand <- unique(geneexons$strand)
  cat(paste("Transcript selected:", transcript, ", Strand:", gene_strand, "\n"))
  
  # invert exon order if minus strand (not used. Better keep real exon order.)
  #geneexons <- geneexons[order(geneexons$exon_chrom_start, decreasing= grepl("-", gene_strand)), ]
  
  # translate exon BP coordinates to aminoacid (AA) coordinates 
  geneexons <- geneexons[order(geneexons$cds_start, decreasing= F), , drop=F]
  geneexons$exon_cds_length <- abs(geneexons$cds_end - geneexons$cds_start) +1 # length coding region
  
  # initialise first row of AA coordinates incl. phase at end of coding sequence (modulo)
  geneexons$AA_start[1] <- geneexons$cds_start[1]
  geneexons$AA_end[1] <- geneexons$AA_start[1] -1 + round(geneexons$exon_cds_length[1] /3)
  geneexons$modulo[1] <- geneexons$exon_cds_length[1] %% 3 
  
  # determine exact AA coordinates
  if(nrow(geneexons)>=2) {
      for (i in 2:nrow(geneexons)) {
        geneexons$AA_start[i] <- geneexons$AA_end[i-1] +1
        geneexons$AA_end[i] <- geneexons$AA_start[i] -1 + floor((geneexons$exon_cds_length[i] + geneexons$modulo[i-1]) /3)
        geneexons$modulo[i] <- (geneexons$exon_cds_length[i] + geneexons$modulo[i-1]) %% 3    
      }
  }  
 
  
  
  # read domain data exported from Uniprot
  prot <- read.table(uniprot_domains.gff, sep="\t", header=F, stringsAsFactors = F)
  prot[, 10] <- NULL
  names(prot) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

  # read total protein length from second comment line of uniprot_domains.gff
  prot_length <- read.table(uniprot_domains.gff, header=F,  stringsAsFactors = F, comment.char = "", skip=1, nrows=1)
  prot$protein_length <- prot_length[, ncol(prot_length)]
  
  prot$feature_length <- prot$end - prot$start

  # find exon interval corresponding with protein domains
  index.interv.start <- findInterval(prot$start, c(geneexons$AA_start, geneexons$AA_end[length(geneexons$AA_end)]) ) 
  index.interv.end   <- findInterval(prot$end,   c(geneexons$AA_start, geneexons$AA_end[length(geneexons$AA_end)]) ) 
  
  if (grepl("-", unique(geneexons$strand))) { 
    # minus strand
    prot$BP_start <- geneexons$genomic_coding_end[index.interv.start] - 3* (prot$start - geneexons$AA_start[index.interv.start])
    prot$BP_end   <- geneexons$genomic_coding_end[index.interv.end]   - 3* (prot$end   - geneexons$AA_start[index.interv.end])
    
  } else {
    # plus strand
      prot$BP_start <- geneexons$genomic_coding_start[index.interv.start] + 3* (prot$start - geneexons$AA_start[index.interv.start])
      prot$BP_end   <- geneexons$genomic_coding_start[index.interv.end]   + 3* (prot$end   - geneexons$AA_start[index.interv.end])
      }
      
  # initialise additional columns
  prot$plotFeature <- TRUE
  prot$domain_name_plot <- prot$feature
  prot$symbol_plot <- "ellipse"
  prot$domain_height_extension <- 1
  prot$domain_color <- "black"
  prot$label_pos <- 3 # Values of 1, 2, 3 and 4 indicate positions below, left, above and right of the specified coordinates.
  prot$assignArrows2Gene <- TRUE
  #prot$domain_width_extension <- prot$feature_length # not used any more
  
  # write result table
  uniprot_domains.gff <- sub(".gff", "", uniprot_domains.gff)
  output.filename <- paste0(uniprot_domains.gff, "_", transcript, suffix.outputFilname)
  cat(paste("Write result table to", output.filename))
  write.table(prot, file=output.filename, sep="\t", quote = F, row.names = F) 
  cat("\nEdit table: select features/domains, symbol ('rectangle', 'ellipse'), color, label_pos (3=top, 1=bottom)) as requested for plotting.\n")
  
  return(geneexons)  
}  
