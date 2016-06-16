\name{QC_expressionset}
\alias{QC_expressionset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quality Control procedures for ExpressionSet data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QC_expressionset(eset, projectfolder = pipepar[["outdirGEX"]], projectname = "",
  groupColumn = pipepar[["groupGEX"]], phDendro = groupColumn,
  flashClustMethod = "average", cex.dendroLabels = 0.6, thresholdZ.k = -5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
ExpressionSet or ExpressionSetIllumina object
}
  \item{projectfolder}{
character vector with name of desired output folder (will be generated if not exisiting).
}
  \item{projectname}{
character vector containing prefix for names of output files.
}
  \item{groupColumn}{
character vector with name of group column in phenotype data.
}
  \item{phDendro}{
character vector with phenotypes of eset object to be displayed in sample Dendrogram
}
  \item{flashClustMethod}{
character vector with the agglomeration method to be used in WGCNA package. Can be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".}
  \item{cex.dendroLabels}{
numeric. Controls size of labels sample in dendrogram.
}
  \item{thresholdZ.k}{
numeric. Threshold for outlier prediction WGCNA package (outlier samples are highlighted in red in sample dendrogram).
}
}
\details{
Function implements several approaches for quality control of expression data, which are stored in a newly generated 'QC'-folder.
  - Illumina internal control plots, box-plots and MA-plots (if applicable).
  - Quality control as implemented in the arrayQualityMetrics-package (subfolder 'quality_metrics')
  - Quality control routines implemented in WGCNA-package. This comprises a sample dendrogram with included outlier prediction as well as boolean vectors for sample and gene quality calculated by the goodSamplesGenes()-function.
}
\value{
The returned value is a list of 3
  \item{goodGenes}{boolean vector indicating genes passing QC-criteria as implemented in WGCNA.}
  \item{goodSamples}{boolean vector indicating samples passing QC-criteria as implemented in WGCNA.}
  \item{allOK}{boolean. If TRUE all samples and all genes pass QC-criteria as implemented in WGCNA.}
Several plots and files are generated as side-effects and stored in the 'QC'-folder.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Frank Ruehle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
arrayqualitymetrics
WGCNA
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
