% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEX_read_Illu_expr_array.R
\name{read_Illu_expr_array}
\alias{read_Illu_expr_array}
\title{Reading expression data from Illumina arrays}
\usage{
read_Illu_expr_array(
  dataFile,
  qcFile,
  sampleSheet,
  ProbeID = "ProbeID",
  skip = 0,
  controlID = "ProbeID",
  qc.skip = 0,
  qc.columns = list(exprs = "AVG_Signal", Detection = "Detection Pval"),
  sampleColumn = "Sample_Name",
  groupColumn = "Sample_Group",
  exprchip = NULL,
  org = NULL,
  covarfile = NULL,
  covarsampleID = "ID",
  method_norm = "none",
  transform = "log2",
  fields2Add = NULL
)
}
\arguments{
\item{dataFile}{character with filepath to SampleProbeProfile from GenomeStudio (tab-delimited txt-file) 
containing columns with \code{AVG_Signal}, \code{Detection Pval}, \code{BEAD_STDERR} und \code{Avg_NBEADS}.
Additional columns are ignored except for \code{"PROBE_ID"}, \code{"SYMBOL"} and \code{"ENTREZ_GENE_ID"}.}

\item{qcFile}{character with filepath to \code{ControlProbeProfile} from GenomeStudio (tab-delimited txt-file)
containing columns with \code{AVG_Signal} und \code{Detection Pval}}

\item{sampleSheet}{character with filepath to Sample sheet from GenomeStudio project (csv-file). The header is expected in row 8.}

\item{ProbeID}{character string with name of the column in dataFile that contains identifiers used to uniquely identify each probe}

\item{skip}{number of header lines to skip at the top of \code{dataFile}.}

\item{controlID}{character string specifying the column in qcFile that contains the identifiers that 
uniquely identify each control probe. If the \code{ControlGeneProfile} from GenomeStudio is used, 
you will need to set \code{controlID="TargetID"}}

\item{qc.skip}{number of header lines to skip at the top of qcFile}

\item{qc.columns}{list defining the column headings in qcFile which correspond to the matrices stored in the 
QCInfo slot of the final ExpressionSetIllumina object}

\item{sampleColumn}{Name of sample column in sample sheet}

\item{groupColumn}{Name of group column in sample sheet}

\item{exprchip}{character string specifying expression chip type (e.g. \code{"HumanHT-12 v4"}, \code{"MouseWG-6 v2"}, \code{"MouseRef-8 v2"}, \code{"RatRef-12 v1"}.}

\item{org}{character string specifying organism: \code{"human"}, \code{"rat"}, \code{"mouse"}.}

\item{covarfile}{character with filepath if a covariates file is given. NULL otherwise.}

\item{covarsampleID}{character with name of sample column in (optional) covariates file.}

\item{method_norm}{character with normalisation method. Options are \code{"quantile"}, \code{"qspline"}, 
\code{"vsn"}, \code{"rankInvariant"}, \code{"median"} and \code{"none"}}

\item{transform}{character with data transformation method. Options are \code{"none"}, \code{"log2"}, \code{"neqc"}, \code{"rsn"} and \code{"vst"}.}

\item{fields2Add}{character vector with names of Illumina mappings to add to feature data (chip type dedicated in \code{exprchip}).
Remark: works only if feature data is \code{"PROBE_ID"}, but this identifier may be non-unique and may therefore throw an error when loading.}
}
\value{
Annotated (and log-transformed) ExpressionSetIllumina object
}
\description{
Create an ExpressionSetIllumina object from Illumina GenomeStudio output files.
}
\details{
Unnormalized expression data from an Illumina GenomeStudio project is read into a \code{ExpressionSetIllumina} object.
If a covariate file is given, covariates are included into the phenotype data of the object. 
If a paired sample design was chosen, the variable indicating sample pairing must be included either in
the sample sheet file or in the covariates file. In default settings, the expression data is log2-transformed
but not normalized.
The feature data of the ExpressionSetIllumina object is supplemented by annotation data given in 
the respective Illumina annotation package (if available \code{REANNOTATED}-columns used). 
The respective chip-type is given in \code{exprchip}.
}
\author{
Frank Ruehle
}
