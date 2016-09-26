\name{read_Illu_Array}
\alias{read_Illu_Array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reading expression data from Illumina arrays
}
\description{
Unnormalized expression data from an Illumina GenomeStudio project is read into a ExpressionSetIllumina object.
}
\usage{
read_Illu_Array(dataFile, qcFile, sampleSheet,
  ProbeID = "PROBE_ID", skip = 0, controlID = "ProbeID",
  qc.skip = 0, qc.columns = list(exprs = "AVG_Signal", Detection = "Detection Pval"),
  sampleColumn = "Sample_Name", groupColumn = "Sample_Group",
  exprchip = "HumanHT-12 v4",
  org = "human",
  covarfile = NULL, covarsampleID = "ID",
  matchvar = NULL,
  method_norm = "none",
  transform = "log2",
  fields2Add = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataFile}{
character with filepath to SampleProbeProfile from GenomeStudio (tab-delimited txt-file) containing columns with AVG_Signal, Detection Pval, BEAD_STDERR und Avg_NBEADS}
  \item{qcFile}{
character with filepath to ControlProbeProfile from GenomeStudio (tab-delimited txt-file) containing columns with AVG_Signal und Detection Pval
  }
  \item{sampleSheet}{
character with filepath to Sample sheet from GenomeStudio project (csv-file). The header is expected in row 8.
}
  \item{ProbeID}{
character string with name of the column in dataFile that contains identifiers used to uniquely identify each probe
}
  \item{skip}{
number of header lines to skip at the top of dataFile.
}
  \item{controlID}{
character string specifying the column in qcFile that contains the identifiers that uniquely identify each control probe
  }
  \item{qc.skip}{
number of header lines to skip at the top of qcFile
}
  \item{qc.columns}{
list defining the column headings in qcFile which correspond to the matrices stored in the QCInfo slot of the final ExpressionSetIllumina object
  }
  \item{sampleColumn}{
Name of sample column in sample sheet
}
  \item{groupColumn}{
Name of group column in sample sheet
}
  \item{exprchip}{
character string specifying expression chip type (e.g. "HumanHT-12 v?", "MouseWG-6 v?", "MouseRef-8 v?")
}
  \item{org}{
character string specifying organism: "human", "rat", "mouse".
}
  \item{covarfile}{
character with filepath if a covariates file is given. NULL otherwise.
}
  \item{covarsampleID}{
character with name of sample column in (optional) covariates file.
}
  \item{matchvar}{
NULL if unpaired study design. In  paired study design, character with variable name indicating pairing of samples. This column must be contained either in sampleSheet or covariates file.
  }
  \item{method_norm}{
character with normalisation method. Options are "quantile", "qspline", "vsn", "rankInvariant", "median" and "none"
}
  \item{transform}{
character with data transformation method. Options are "none", "log2", "neqc", "rsn" and "vst".
}
  \item{fields2Add}{
character vector with names of Illumina mappings to add to feature data (chip type dedicated in 'exprchip').
}
}
\details{
If a covariate file is given, covariates are included in the phenotype data of the object. If a paired sample design was chosen, the variable indicating sample pairing must be included either in the sample sheet file or in the covariates file. In default settings, the expression data is log2-transformed but not normalized. The feature data of the ExpressionSetIllumina object is supplemented by annotation data given in the respective Illumina annotation package (if available "REANNOTATED"-columns used). The respective chip-type is given in 'exprchip'.
  }
\value{
Annotated (and log-transformed) ExpressionSetIllumina object
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
