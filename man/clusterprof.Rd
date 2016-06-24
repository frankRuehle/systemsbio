\name{clusterprof}
\alias{clusterprof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper for enrichment analysis using clusterProfiler
}
\description{
Function uses genelist and optionally background gene list as input to perform enrichment analysis using the 'clusterProfiler'-package.
}
\usage{
clusterprof(genes, newheader = NULL,
  backgroundlist = NULL, newheaderBackground = NULL,
  projectfolder = pipepar[["outdirGEX"]], projectname = "",
  enrichmentCat = c("GO", "KEGG", "Reactome", "DO"),
  maxInputGenes = 100,
  id.type = "SYMBOL", id.column = "SYMBOL",
  sortcolumn = "adj.P.Val", sortdecreasing = FALSE,
  sortcolumn.threshold = 0.05, fun.transf.incr.vales = function(x) {-log10(x)},
  FCcolumn = "logFC", threshold_FC = log2(pipepar[["threshold_FC"]]),
  org = "human", pAdjustMethod = "BH", enrich.p.valueCutoff = 0.05,
  enrich.q.valueCutoff = 0.05, nPerm = 1000, minGSSize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
vector with gene names to analyse or dataframe or character with path to dataframe. Supplied dataframe needs to include a column with identifiers specified in 'id.column'.
  }
  \item{newheader}{
optional character vector with new header information for 'genes' dataframe. Only relevant if 'genes' is a dataframe (or character string with filepath to a table) with wrong or missing header. NULL otherwise.
}
  \item{backgroundlist}{
optional background list for enrichment analysis. Can either be a vector with gene IDs or dataframe or character with path to dataframe. If given, an 'id.column' needed as in 'genes'. If "genome", all ENTREZ IDs from the annotation package of the respective org (denoted in 'org') are used as background. If NULL, full ID list from 'genes' is used as background.
  }
  \item{newheaderBackground}{
optional character vector with new header information for 'backgroundlist' dataframe.
}
  \item{projectfolder}{
character with directory for output files (will be generated if not exisiting).
}
  \item{projectname}{
optional character prefix for output file names.
}
  \item{enrichmentCat}{
character vector with categories to be enriched (GO: gene ontology (MF, BP, CC), KEGG: KEGG pathways, Reactome: Reactome pathways, DO: Disease ontology).
}
  \item{maxInputGenes}{
(numeric) max number of top diff regulated elements used for enrichment analysis.(or NULL).
}
  \item{id.type}{
character with identifier type from annotation package (e.g. "ENTREZID", "SYMBOL", "UNIGENE") IDs Will be converted to EntrezIDs prior to enrichment analysis.
  }
  \item{id.column}{
character with column name for identifier variable in 'genes'.
}
  \item{sortcolumn}{
character with column name of quantitative data in 'genes' used for ordering. If Null, ranking of genes is omitted and GSEA not possible.
  }
  \item{sortdecreasing}{
(boolean) order parameter for hierarchy of values in 'sortcolumn'. FALSE for increasing values (e.g. p-values), TRUE for decreasing values (e.g. fold changes). If FALSE, values in sortcolumn will be transformed prior to GSEA.
  }
  \item{sortcolumn.threshold}{
numeric threshold for 'sortcolumn' to be included in overepresentation analysis. If sortdecreasing=F, value < sortcolumn.threshold else value > sortcolumn.threshold.
  }
  \item{fun.transf.incr.vales}{
function definition for transforming values in 'sortcolumn' prior to GSEA. Is applied if sortdecreasing=F.
  }
  \item{FCcolumn}{
(character) optional column name of foldchanges in 'genes' if 'sortcolumn' is used elsewhere. Used only for cnetplot of enrichment results. Omitted if NULL
  }
  \item{threshold_FC}{
(numeric) Fold change threshold for filtering (threshold interpreted for log transformed foldchange values!). Only relevant for overrepresentation analysis if an unfiltered gene list is given in 'genes' to allow for parallel GSEA.
  }
  \item{org}{
character with name of organism ("human", "mouse", "rat").
}
  \item{pAdjustMethod}{
method for adjusting for multiple testing. One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
}
  \item{enrich.p.valueCutoff}{
numeric p-value threshold for returned enrichment terms.
}
  \item{enrich.q.valueCutoff}{
numeric q-value threshold for returned enrichment terms.
}
  \item{nPerm}{
permutation numbers
}
  \item{minGSSize}{
 minimal size of genes annotated by Ontology term for testing.
 }
}
\details{
 By default, overrepresentation analysis and gene set enrichment analysis (GSEA) is performed for all categories given in 'enrichmentCat'. The Gene list can be supplied as vector with ids given in 'id.type', dataframe with IDs of type 'id.type' given in column 'id.column' or character string with path of data file. The IDs are converted to ENTREZ IDs (if necessary) prior to enrichment using the annotation package for the species denoted in 'org'. Optionally, quanitative data can be included in 'sortcolumn' for sorting and filtering (using sortcolumn.threshold) the data. If quantitative data are no fold changes, fold changes may be given additionally in 'FCcolumn'. They are used only for cnetplot of enrichment results only. If Null the data in sortcolumn is also used for the cnetplots. Be aware that the legend of the cnetplots will be "Fold Change" anyway! If no quantitative data is provided, sorting and GSEA will be skipped. Optionally, a background list can be provided for enrichment analysis.
 }
\value{
List of enrichment-objects defined in 'DOSE'-packagage ('enrichResult'-object for overrepresentation analysis, 'gseaResult'-objects for gene set enrichment analysis). Enrichment tables and plots are stored in the project folder as side effects.
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
clusterProfiler
DOSE
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
