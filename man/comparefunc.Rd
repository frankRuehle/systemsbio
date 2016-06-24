\name{comparefunc}
\alias{comparefunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functional analysis of gene sets
%%  ~~function to do ... ~~
}
\description{Use dataframe (or list of dataframes) with gene sets of interest for functional analysis such as over-representation analysis or GSEA
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comparefunc(DEgenes.unfilt,
            newheader=NULL,
            comparisons=names(DEgenes.unfilt),
            projectfolder=file.path(pipepar[["outdir"]], "GEX"),
            projectname="",
            pvalue = "P.Value",
            adj.pvalue = "adj.P.Val",
            p.value.threshold = pipepar$threshold_p,
            effectsize = "logFC",
            effect.threshold = log2(pipepar$threshold_FC),
            Symbol.column = "SYMBOL",
            Entrez.column = "ENTREZID",

            # Volcano plot parameter:
            volcanoP = adj.pvalue,

            # Enrichment analysis parameter:
            use.clusterProfiler = TRUE,
            use.davidAccount = FALSE,
            use.TFsearch = TRUE,
            maxInputGenes = 100,
            enrichmentCat = c("GO", "KEGG", "Reactome", "DO"),
            org = "human",

            # Further Subroutines
            use.pwOmics = FALSE   # not implemented yet
            )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DEgenes.unfilt}{
dataframe or character with path to dataframe or list of dataframes with unfiltered data from diff expression / methylation analysis (e.g. results from groupsLimma-function).Required columns with p-value, adjusted p-value and effect size are defined below.}

\item{newheader}{
optional character vector with new header information for 'DEgenes.unfilt' dataframe. Only relevant if 'DEgenes.unfilt' is a dataframe or character string with filepath to a table with wrong or missing header. NULL otherwise.
}

\item{comparisons}{
character vector with respective names of group comparisons in format "groupA-groupB". As default, 'DEgenes.unfilt' is expected to be list of dataframes with the respective group comparison as name of each list element. Otherwise give here the group comparisons matching the data in 'DEgenes.unfilt'.
}

\item{projectfolder}{
character with directory for output files (will be generated if not exisiting).
}

\item{projectname}{
optional character prefix for output file names.
}

\item{pvalue}{
character with column name of p values in 'DEgenes.unfilt'
}

\item{adj.pvalue}{
character with column name with p value adjusted for multiple testing in 'DEgenes.unfilt'
}

\item{p.value.threshold}{
numeric p-value threshold (is applied to p-values corrected for multiple testing in column 'adj.pvalue')
}

\item{effectsize}{
character with column name of effectsize in 'DEgenes.unfilt'
}

\item{effect.threshold}{
numeric foldchange threshold
}

\item{Symbol.column}{
character with column name of Gene Symbols in 'DEgenes.unfilt'. If NULL, SYMBOLs are derived from 'Entrez.column' using the annotation package for 'org'.
}

\item{Entrez.column}{
character with column name of ENTREZ IDs in 'DEgenes.unfilt'. If NULL, ENTREZ IDs are derived from 'Symbol.column' using the annotation package for 'org'. Either 'Symbol.column' or 'Entrez.column' must be specified.
}

\item{volcanoP}{
specify if adjusted (adj.pvalue) or unadjusted (pvalue) p-values shall be used for volcano plots.
}

\item{use.clusterProfiler}{
(boolean) use clusterProfiler for enrichment analysis.
}

\item{use.davidAccount}{
email account for DAVID's Web Service connectivity or DAVIDWebService object. If NULL or FALSE, DAVID enrichment analysis is skipped.
}

\item{use.TFsearch}{
(boolean) use TFsearch.R for TF motive enrichment (human only)
}

\item{maxInputGenes}{
(numeric) max number of top diff regulated elements used for overrepresentation analysis.
}

\item{enrichmentCat}{
character vector with DAVID categories to enrich for. E.g. c("GO", "KEGG", "Reactome", "DO") for clusterProfiler or c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", "UCSC_TFBS") for DAVID. Combinations possible if clusterProfiler and DAVID are used.
}

\item{org}{
character with species name ("human", "mouse", "rat").
}

\item{use.pwOmics}{
(boolean) use pwOmics package for analysis (human only) (not implemented yet)
}



}
\details{
The function uses data.frames as input with differentially regulated elements (e.g. genes or methylation sites) resulting from a specified group comparison. Data can be delivered as single dataframe, as list of dataframes named by group comparisons or as character containing a filepath to be loaded. Either column names for Symbols or Entrez IDs must be specified. For each group comparison, volcano plots are generated and stored in the designated output folder. Optional enrichment analsis (overrepresentation or gene set enrichment) is performed on ENTREZ IDs for specified categories and/or transcription factor binding sites using the clusterProf(), DavidEnrich() and TFsearch() subroutines.
}

\value{
No value returned. Output files are stored as side-effects.
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

%## The function is currently defined as
%function (x)
%{
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
