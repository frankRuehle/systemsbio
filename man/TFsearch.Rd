\name{TFsearch}
\alias{TFsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transcription factor binding sites enrichment
}
\description{
Function takes input genes, optionally looks up all promotor sequences by (unique) entrezID, or looks up sequences defined by coordinates, optionally filters 'sequences' for given categories and calls motifEnrichment() from 'PWMEnrich'-package to identify enriched Transcription factor binding motivs. Optionally tables generated for in 'SearchSelMotifs' preselected Motivs looked up in 'sequences'.
  }
\usage{
TFsearch(sequences, newheader = NULL, annoColumn = NULL, name.organism = "hsapiens",
  projectfolder = "GEX/TFBS", projectname = "",
  applyFilter = FALSE, filtercat1 = "adj.P.Val", filtercat1.decreasing = FALSE,
  filtercat1.function = abs, filtercat1.threshold = 0.05,
  filtercat2 = "logFC", filtercat2.decreasing = TRUE,
  filtercat2.function = abs, filtercat2.threshold = log2(1.5),
  PromLookup = TRUE, Entrez.col = "ENTREZID", PromSeqUpstreamTSS = 2000, PromSeqDownstreamTSS = 200,
  SearchSelMotifs = NULL, motif.min.score = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{
dataframe (or character with file path) with column containing Entrez IDs or sequence coordinates. Required columns as defined in 'PromLookup' (see below).
}
  \item{newheader}{
NULL if 'sequences' already supplied with header. Character vector with new header otherwise.
}
  \item{annoColumn}{
character or vector of characters. Column name(s) of 'sequences'-object with sequence annotation to maintain
}
  \item{name.organism}{
currently human data only (hg19)
}
  \item{projectfolder}{
output directory
}
  \item{projectname}{
character merged to output name
}
  \item{applyFilter}{
(boolean) If TRUE, sequences are filtered for applied categories and thresholds. Filter Values converted to ABSOLUTE values. Optional Filtering criteria (Ignored if applyFilter==FALSE).
  }
  \item{filtercat1}{
column name of first category to filter 'sequences' (e.g. p-values).
}
  \item{filtercat1.decreasing}{
(boolean) direction to order filtercat1
}
  \item{filtercat1.function}{
select transforming function for filter category1 (no quotes). e.g. abs for absolute values, identity for no transformation
}
  \item{filtercat1.threshold}{
Threshold for filtercat1 or 'top123' for top Hits
}
  \item{filtercat2}{
column name of second category to filter 'sequences' (e.g. effect size).
}
  \item{filtercat2.decreasing}{
(boolean) direction to order filtercat2
}
  \item{filtercat2.function}{
select transforming function for filter category2 (no quotes). E.g. abs for foldchanges
}
  \item{filtercat2.threshold}{
Threshold for filtercat2 or 'top123' for top Hits
}
  \item{PromLookup}{
if TRUE, all promotor sequences corresponding to genes in 'sequences' are downloaded. Therefore a column with EntrezIDs is requried (column name given in 'Entrez.col')! if FALSE, Sequences are downloaded according to given coordinates in 'sequences'. Therefore columns for chromosome, start, stop and strand information required! Additional meta columns allowed.
}
  \item{Entrez.col}{
character with name of column containing Entrez IDs
}
  \item{PromSeqUpstreamTSS}{
definition of promotor regions to download respective to TSS
}
  \item{PromSeqDownstreamTSS}{
definition of promotor regions to download respective to TSS
}
  \item{SearchSelMotifs}{
Character Vector of selected motives to search in 'sequences' separately to enrichment or NULL.
}
  \item{motif.min.score}{
minimum score to match motif pwm to target sequence (ignored if SearchSelMotifs = NULL)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
groupReport of motifEnrichment results.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Frank Ruehle}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
MotifDb
PWMEnrich
PWMEnrich.Hsapiens.background
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
