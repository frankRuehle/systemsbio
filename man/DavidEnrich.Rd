\name{DavidEnrich}
\alias{DavidEnrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DAVID enrichment analysis using the RDAVIDWebService
}
\description{
Wrapper function which uses a genelist and optionally background gene list as input to perform enrichment analysis using the RDAVIDWebService.
}
\usage{
DavidEnrich(genes, newheader = NULL,
  backgroundlist = NULL, newheaderBackground = NULL,
  davidAccount = NULL,
  projectfolder = "GEX/DAVID_Enrichment", projectname = "",
  maxDavid = 100, davidCat = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", "UCSC_TFBS"),
  id.type = "ENTREZ_GENE_ID", id.column = "ENTREZID",
  sortcolumn = "adj.P.Val", sortdecreasing = FALSE, sortcolumn.threshold = 0.05,
  org = "Homo sapiens", threshold4enrich = 0.05, mincount = 2L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
vector with gene names to analyse or dataframe or character with path to dataframe. Supplied dataframe needs to include a column with identifiers specified in 'id.column'.
}
  \item{newheader}{
optional character vector with new header information for 'genes' dataframe. Only relevant if 'genes' is a dataframe or character string with filepath to a table with wrong or missing header. NULL otherwise.
  }
  \item{backgroundlist}{
optional dataframe containing background list for enrichment analysis. If NULL, species-specific background from DAVID is used instead.
  }
  \item{newheaderBackground}{
optional character vector with new header information for 'backgroundlist' dataframe.
}
  \item{davidAccount}{
email account for DAVID's Web Service connectivity or DAVIDWebService object. Registration for RDAVIDWebService is required at "https://david.ncifcrf.gov/webservice/register.htm".
  }
  \item{projectfolder}{
character with directory for output files (will be generated if not exisiting).
}
  \item{projectname}{
optional character prefix for output file names.
}
  \item{maxDavid}{
(numeric) max number of top diff regulated elements used for enrichment analysis.
}
  \item{davidCat}{
character vector with DAVID categories to enrich for.
}
  \item{id.type}{
character with identifier DAVID uses to enrich for (e.g. "ENTREZ_GENE_ID" or "OFFICIAL_GENE_SYMBOL")
}
  \item{id.column}{
character with respective column name for id.type in DEgenes.unfilt'.
}
  \item{sortcolumn}{
character with column name to be sorted for selecting top entries.
}
  \item{sortdecreasing}{
(boolean) drecreasing parameter for order of 'sortcolumn'.
}
  \item{sortcolumn.threshold}{
numeric threshold for 'sortcolumn' to be included in overepresentation analysis. If sortdecreasing=F, value < sortcolumn.threshold else value > sortcolumn.threshold
  }
  \item{org}{
character with species name DAVID uses for enrichment (e.g."Homo sapiens" or "Mus musculus").
}
  \item{threshold4enrich}{
numeric p-value threshold for enrichment terms.
}
  \item{mincount}{
numeric minimum gene count in enrichment terms.
}
}
\details{
Gene list can be supplied as vector with ids given in id.type, dataframe with or character string with path of data file. If a vector is given, sorting parameter will be ignored. Enrichment categories are defined in 'davidCat' and species is defined in 'org'. If no connection to RDAVIDWebService can be established, NULL is returned.
}
\value{
Functional annotation chart as R object. Additionally, a functional annotation chart is stored in the project folder as side effect.
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
RDAVIDWebService
GSEABase
KEGG.db
KEGGgraph
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
