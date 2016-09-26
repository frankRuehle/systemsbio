\name{wrapPCAgoprom}
\alias{wrapPCAgoprom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper for PCA routines from 'pcaGoPromoter'-package.
}
\description{
2-dim and 3-dim PCA plots are generated for desired samples in the given ExpressionSet (expca). Tables of PC-associated probes and transcription factors and GO terms enriched in top associated probes are generated for any number of principal components in positive and negative orientation. All output data is stored in 'projectfolder'.
  }
\usage{
wrapPCAgoprom(expca, groupsoi = NULL,
  groupby = "Sample_Group", sample.name.column = "Sample_Name", samples2exclude = NULL,
  projectfolder = file.path("GEX", "pcaGoPromoter"), projectname = NULL,
  inputType = "geneSymbol", org = "Hs",
  PCs4table = 2, PCs2plot = c(1, 2, 3), probes2enrich = 0.025)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expca}{
ExpressionSet object or a table with expression data with variables (probes) in rows and observations in columns (samples).
}
  \item{groupsoi}{
character vector with sample groups of interest to be included in PCA. Respective samples are taken from phenoData of expca. groupnames must match entries in column given in 'groupby'
  }
  \item{groupby}{
Column name of phenoData of expca used for group names
}
  \item{sample.name.column}{
Column name of phenoData of expca used for sample names
}
  \item{samples2exclude}{
Character vector for optionally exclusion of individual samples. Used as regular expression for lookup of samples.
  }
  \item{projectfolder}{
character with directory for output files (will be generated if not exisiting).
}
  \item{projectname}{
optional character prefix for output file names.
}
  \item{inputType}{
a character vector description of the input type. Must be Affymetrix chip type, "geneSymbol" or "entrezID".
  }
  \item{org}{
a character vector specifying the organism. Either "Hs" (homo sapiens), "Mm" (mus musculus) or "Rn" (rattus norwegicus).
}
  \item{PCs4table}{
numeric or numeric vector. Indicates number of PCs (numeric) or distinct PCs (numeric vector) for which result tables of enriched transcription factor binding sites and GO-terms are calculated.
  }
  \item{PCs2plot}{
numeric or numeric vector. Indicates number of PCs (numeric) or distinct PCs (numeric vector) to use in 2-dim and 3-dim PCA plots. For 2-dim PCA plots all possible pairs of PCs are plotted. Additionally, a 3D plot is generated with the first 3 PCs in PCs2plot. Note that pca informative plot (containing TFBS and GO annotation on the axes) is restricted to first two PCs only!
  }
  \item{probes2enrich}{
numeric. Number of top PC-associated probes to look for enriched TFBS and GO terms. A value <= 1 is interpreted as fraction of total number of probes.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Several plots and files are generated as side-effects and stored are in the designated projectfolder. The returned value is a list of 3 lists.
  \item{loadsperPC }{Top associated probes for every PC in pos and neg direction}
  \item{TFtables }{dataframes containing enriched TFBS for every PC in pos and neg direction (over- and underrepresented)}
  \item{GOtreeOutput }{dataframes containing enriched GO terms for every PC in pos and neg direction}
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
pcaGoPromoter
pcaGoPromoter.Hs.hg19
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
