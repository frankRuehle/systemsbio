\name{setparam}
\alias{setparam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set parameters used in systemsbio pipeline }
\description{
This functions optionally organizes parameter for analysis settings frequently used within analysis pipeline and returns a list with relevant parameters, which may be used as standard input for further functions. The denoted working directory is set and required analysis and annotation packages are installed and/or loaded.}
\usage{
setparam(workingdir, outdir = "Analysis", projectname = "",
  org = "human", threshold_p = 0.05, threshold_fdr = 0.05, threshold_FC = 1.5,
  exprchip = NULL, exprchip.manufacturer = "Illumina", outdirGEX = file.path(outdir, "GEX"),
  sampleGEX = "Sample_Name", groupGEX = "Sample_Group", matchvarGEX = NULL,
  covarGEXSampleID = "ID", groupcomparisonsGEX = c("case-control"),
  methPlatform = NULL, methPlatform.manufacturer = "Illumina", methannotation.file = NULL,
  outdirMT = file.path(outdir, "MT"), threshold_beta = 0.1, refgroupMeth = "control", typeMethvar = "categorical",
  sampleMT = "Sample_Name", groupMT = "Sample_Group", matchvarMeth = "none",
  covarMTSampleID = "IID", groupcomparsionsMT = groupcomparisonsGEX,
  snpchip = NULL, snpchip.manufacturer = "Illumina", outdirGT = file.path(outdir, "GT"), covarGTSampleID = "IID",
  updateBiocPackages = FALSE, source.bioc = "http://bioconductor.org/biocLite.R")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{workingdir}{
character with working directory.
}
  \item{outdir}{
character with output directory within working directory for result files (will be created if not existing).
}
  \item{projectname}{
character used as suffix for output file names.
}
  \item{org}{
character with species name. Currently "human", "rat" or "mouse". Corresponding annotation packages will be loaded.
}
  \item{threshold_p}{
numeric p-value threshold.
}
  \item{threshold_fdr}{
numeric threshold for FDR (false discovery rate) corrected p-value.
}
  \item{threshold_FC}{
numeric foldchange threshold.
}

# GEX parameter (only relevant if exprchip != NULL)
  \item{exprchip}{
character with expression chip type, e.g. "HumanHT-12 v4", "MouseWG-6 v2", "MouseRef-8 v1".
}
  \item{exprchip.manufacturer}{
character with chip manufacturer, e.g. "Illumina" or "Affymetrix". Will be used to load corresponding annotation package. NULL if no supplier specified.
  }
  \item{outdirGEX}{
character with output directory for gene expression results.
}
  \item{sampleGEX}{
character with column name containing sample names in sample data (e.g. in Illumina sample sheet).
}
  \item{groupGEX}{
character with column name containing group names in sample data (e.g. in Illumina sample sheet).
}
  \item{matchvarGEX}{
In case of matched study design, character with column name indicating corresponding samples. If NULL, unmatched study design assumed.
  }
  \item{covarGEXSampleID}{
character with column name of sample names in optional covariate file.
}
  \item{groupcomparisonsGEX}{
character vector with desired group comparisons to analyse in format "groupA-groupB" or nested comparisons in format "(groupA-groupB)-(groupC-groupD)". Group names must correspond to data given by 'groupGEX'. Analysis can be performed either for any number of unpaired group comparisons or for one group comparison in paired design.}

  # Methylation parameter (only relevant if methPlatform != NULL)
  \item{methPlatform}{
character with methylation array platform, e.g. "450k-UCSC" for UCSC CpG Islands for 450k array probes, "450k-HMM" for HMM CpG Islands for 450k array probes, "27k" for UCSC CpG Islands for 27k array probes.
  }
  \item{methPlatform.manufacturer}{
character with chip manufacturer, e.g. "Illumina".
}
  \item{methannotation.file}{
(if no methPlatform specified) character with path to custom tab-text file for annotation with header "SiteID  Chr  Loc  Gene  Island".
  }
  \item{outdirMT}{
character with output directory for gene methylation results
}
  \item{threshold_beta}{
numeric threshold for beta value. beta is defined as: methylated signal / (methylated + unmethylated signal).
}
  \item{refgroupMeth}{
character with reference group of methylation analysis.
}
  \item{typeMethvar}{
character with type of methylation variable ("categorical" or "continuous")
}
  \item{sampleMT}{
character with column name containing sample names in sample data.
}
  \item{groupMT}{
character with column name containing group names in sample data.
}
  \item{matchvarMeth}{
In case of matched study design, character with column name indicating corresponding samples. "none" if unpaired Samples
}
  \item{covarMTSampleID}{
character with column name of sample names in optional covariate file.
}
  \item{groupcomparsionsMT}{
character vector with desired group comparisons to analyse.
}

# GT Parameter (only relevant if snpchip != NULL)
  \item{snpchip}{
character with genotyping chip type, e.g."HumanOmniExpress-24", "HumanOmni5-Quad", "HumanOmni5Exome", "HumanCore-24", "HumanCoreExome-24".
}
  \item{snpchip.manufacturer}{
character with chip manufacturer, e.g. "Illumina".
}
  \item{outdirGT}{
character with output directory for genotyping results.
}
  \item{covarGTSampleID}{
character with column name of sample names in optional covariate file.
}

# other parameter
  \item{updateBiocPackages}{
boolean. if TRUE automatically updates all installed Bioconductor Packages.
}
  \item{source.bioc}{
character with bioconductor source script.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List with parameters given or derived from input arguments. A log-file showing these parameter is stored in the output directory as side-effect.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Frank Ruehle}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
