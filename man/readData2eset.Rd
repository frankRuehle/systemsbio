\name{readData2eset}
\alias{readData2eset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make ExpressionSet object from custom data
}
\description{
Expression data, phenotype data and optional feature data are used to make an ExpressionSet object.
}
\usage{
readData2eset(exprsData, phenoData, featureData = NULL,
  ProbeID = "PROBE_ID", sampleColumn = "Sample_Name", groupColumn = "Sample_Group",
  experimentData = MIAME(),
  exprchip = pipepar$ArrayAnnotation.GEX,
  transform = log2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprsData}{
dataframe or character with filepath to expression data to be loaded. If no rownames are included, a column containing ProbeIDs is expected and used as feature names.
  }
  \item{phenoData}{
dataframe or character with filepath to sample sheet / phenotype data to be loaded.
}
  \item{featureData}{
dataframe or character with filepath to feature data to be loaded. If NULL, non-numeric data columns from 'exprsData' are used as feature data (if any).
}
  \item{ProbeID}{
character with name of the column in 'exprsData' and 'featureData' that contains identifiers used to uniquely identify each probe.
}
  \item{sampleColumn}{
character with name of sample column in 'phenoData'.
}
  \item{groupColumn}{
character with name of group column in 'phenoData'.
}
  \item{experimentData}{
MIAME object with optional experiment description.
}
  \item{exprchip}{
character string specifying chip annotation package (e.g. "Humanv4"). character() if not available.
}
  \item{transform}{
optional function definition to transform expression data (e.g. log2 for logarithm of base 2). If NULL, no data transformation is performed.
  }
}
\details{
Expression data, phenotype data and optional feature data are used to make an ExpressionSet object. Only numeric columns in 'exprsData' are used. If no feature data is given, non-numeric data columns from 'exprsData' are used as feature data (if any). If dataframes don't contain rownames as identifier, the column 'ProbeID' is used as rownames in 'exprsData' and 'featureData' and column 'sampleColumn' in 'phenoData'. Rownames of 'exprsData' must match with rownames of 'featureData' (if given). Column names of 'exprsData' must match with rownames of 'phenoData.' Characters "/" and "-" in sample names or group names are replaced by ".". If either gene symbols or EntrezIDs are missing in feature data, this annotation is added using the corresponding annotation package (species name is derived from 'exprchip'). The (optional) Function given in 'transform' is used to transform the expression data.
}
\value{
Annotated and (optionally) log2-transformed ExpressionSet object
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
