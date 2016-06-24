\name{diffLimma}
\alias{diffLimma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential expression analysis with limma}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diffLimma(GEXMTSet, matchvar = NULL,
  comparisons = pipepar[["groupcomparisonsGEX"]],
  p.value.threshold = pipepar$threshold_p, adjust.method = "BH",
  FC.threshold = log2(pipepar$threshold_FC),
  projectfolder = file.path(pipepar[["outdir"]], "GEX"),
  projectname = "", Symbol.column = "SYMBOL",
  sampleColumn = pipepar$sampleGEX, groupColumn = pipepar$groupGEX,
  useWeights = TRUE,
  geneAnno2keep = NULL,
  venn_comparisons = pipepar[["groupcomparsionsGEX"]],
  maxHM = 50, HMcexRow = 1.1, HMcexCol = 1.1,
  HMincludeRelevantSamplesOnly = TRUE,
  color.palette = "heat.colors",
  FC.heatmap.comparisons = NULL,
  FC.heatmap.geneselection = c("allprobes", "intersect", "union"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GEXMTSet}{
ExpressionSet
}
  \item{matchvar}{
(character) if paired design, name of column to be used to group samples. NULL for unpaired design.
}
  \item{comparisons}{
character vector with group comparisons in format "groupA-groupB" or nested comparisons in format "(groupA-groupB)-(groupC-groupD)".
  }
  \item{p.value.threshold}{
numeric p-value threshold
}
  \item{adjust.method}{
adjustment method for multiple testing ("none", "BH", "BY" and "holm")
}
  \item{FC.threshold}{
numeric foldchange threshold. If data is log2-transformed, mind to transform threshold too.
}
  \item{projectfolder}{
character with directory for output files (will be generated if not exisiting).
}
  \item{projectname}{
optional character prefix for output file names.
}
  \item{Symbol.column}{
character with column name of Gene Symbols in 'DEgenes.unfilt' and 'GEXMTSet' or NULL.
}
  \item{sampleColumn}{
character with column name of Sample names in pData(GEXMTSet)
}
  \item{groupColumn}{
character with column name of group names in pData(GEXMTSet). Group names must match comparisons!
}
  \item{useWeights}{
(boolean) if TRUE, arrays are weighted for their quality.
}
  \item{geneAnno2keep}{
(character vector) annotation columns in feature data to keep in output
}
  \item{venn_comparisons}{
character vector or (named) list of character vectors with group comparisons to be included in Venn Diagram (max 5 comparisons per Venn diagramm)
}
  \item{maxHM}{
(numeric) max number of top diff. regulated elements printed in Heatmap
}
  \item{HMcexRow}{
labelsize for rownames in Heatmap
}
  \item{HMcexCol}{
labelsize for colnames in Heatmap
}
  \item{HMincludeRelevantSamplesOnly}{
(boolean) if TRUE include only Samples in heatmap, which are in groups of the respective group comparison.
  }
  \item{color.palette}{
select color palette for heatmaps
}
  \item{FC.heatmap.comparisons}{
character vector or (named) list of character vectors with group comparsions to be included in foldchange heatmaps. One FC heatmap is generated for each character vector. If NULL no FC heatmaps are generated.
  }
  \item{FC.heatmap.geneselection}{
character vector containing elements of c("all", "intersect", "union"). Which probes shall be used for foldchange heatmaps? "allprobes": no restriction, "intersect": only probes are used which are differentially expressed in all respective group comparisons. "union": only probes are used which are differentially expressed in any of the respective group comparisons.Probes are prioritised by F-Test p-value calculated for all group comparisons as well as by minimum adjusted p-value of respective group comparisons. Obsolete if 'FC.heatmap.comparisons' is NULL.
}
}
\details{
Function generates design matrix, fit and differential expression results for dedicated group comparisons of input object. Analysis is performed either for one group comparison in paired design or for any number of unpaired group comparisons. The group designations in the groupColumn must match the contrasts given in comparisons ("groupA-groupB"). P-value and foldchange thresholds may be applied to filter result data. Heatmaps with sample signal intensities are generated for top differentially expressed genes for each group comparison given in 'comparisons'. Additionally heatmaps indicating foldchanges are generated for each set of group comparisons given in 'FC.heatmap.comparisons'. Selection of probes with respect to these group comparisons is characterised in 'FC.heatmap.geneselection'. If the resulting number of probes exceeds 'maxHM', probes are prioritized by either F-Test p-value or by minimum p-value of selected group comparisons.}
\value{
List of dataframes with unfiltered differential expression data for each comparison (sorted by p-value). List elements are named by the respective group comparison. An F-Test and a Venn Diagramm of all comparisons as well as filtered and unfiltered result tables and heatmaps for each group comparison are stored as side-effects in the project folder.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Frank Ruehle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
