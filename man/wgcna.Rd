\name{wgcna}
\alias{wgcna}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper for 'WGCNA'-package for Whole Genome Coexpression Network Analysis.
}
\description{
This function makes use of the 'WGCNA'-package from Steve Horvath and  Peter Langfelder to construct whole genome coexpression networks.
  }
\usage{
wgcna(GEXMTSet,
  projectfolder = file.path(pipepar[["outdir"]], "GEX"),
  softThresholdPower = "auto",
  corType = "bicor",
  networkType = "signed",
  TOMType = "signed",
  maxBlockSize = 45000,
  TOMplot = FALSE, MDSplot = FALSE,
  phDendro = NULL, phModule = NULL,
  sampleColumn = pipepar$sampleGEX, groupColumn = pipepar$groupGEX,
  groupsets = NULL, symbolColumn = NULL,
  flashClustMethod = "average", dendroRowText = T, autoColorHeight = FALSE,
  colorHeight = 0.1, cex.dendroLabels = 0.6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GEXMTSet}{
ExpressionSet or MethylSet. If 'GEXMTSet' is character containing a filepath, the functions assumes previously stored network object to be loaded from this path. If 'GEXMTSet' is "load_default", network object is loaded from default directory "file.path(projectfolder, "WGCNA", "TOM", "networkConstruction-auto.RData")".
  }
  \item{projectfolder}{
character with directory for output files (will be generated if not exisiting).
}
  \item{softThresholdPower}{
soft-thresholding power for network construction. If "auto", function selects soft-thresholding power automatically. If Null, network construction is omitted.
  }
  \item{corType}{
character string specifying the correlation to be used. Allowed values are "pearson" and "bicor", corresponding to Pearson and bidweight midcorrelation, respectively. Missing values are handled using the pairwise.complete.obs option.
  }
  \item{networkType}{
character with network type. Allowed values are "unsigned", "signed", "signed hybrid". "unsigned" means negative correlation of genes are treated the same as positive correlation. In an "signed" network, negatively correlated genes will not be put into one module, but will be treated as not correlated.
  }
  \item{TOMType}{
character with one of "none", "unsigned", "signed". If "none", adjacency will be used for clustering. If "unsigned", the standard TOM will be used (more generally, TOM function will receive the adjacency as input). If "signed", TOM will keep track of the sign of correlations between neighbours.
  }
  \item{maxBlockSize}{
integer giving maximum block size for module detection. If the number of genes in 'GEXMTSet' exceeds maxBlockSize, genes will be pre-clustered into blocks whose size should not exceed maxBlockSize (maxBlockSize must not exceed 46340). It's intended to use as big block sizes as possible, but mind that big blocksizes will heavily impact memory usage.
  }
  \item{TOMplot}{
boolean. If TRUE make Topological Overlap Matrix (TOM) plot (also known as connectivity plot) of the network connections. Light color represents low topological overlap and progressively darker red color represents higher overlap. Modules correspond to red squares along the diagonal.
  }
  \item{MDSplot}{
boolean. If TRUE make Multidimensional scaling plot (MDS) to visualize pairwise relationships specified by a dissimilarity matrix. Each row of the dissimilarity matrix is visualized by a point in a Euclidean space. Each dot (gene) is colored by the module assignment.
  }
  \item{phDendro}{
character vector with phenotypes of 'GEXMTSet' object to be displayed in sample dendrogram.
}
  \item{phModule}{
character vector with phenotypes to correlate module eigengenes with in heatmap.
}
  \item{sampleColumn}{
character with column name of Sample names in pheno data of 'GEXMTSet'.
}
  \item{groupColumn}{
character with column name of group names in pheno data of 'GEXMTSet'.
}
  \item{groupsets}{
character vector with names of group sets in format "groupA-groupB". Groups summarized in parentheses "(groupA-groupB)" are coded as ONE group. They are used for correlation of module eigengenes with corresponding samples of selected groupsets. Mind that eigengenes are calculated using all samples, while correlation is calculated for samples of denoted groupsets only.Group names must match names in 'groupColumn'. Omitted if NULL.
  }
  \item{symbolColumn}{
character with name of feature identifier in feature data of 'GEXMTSet'.
}
  \item{flashClustMethod}{
character with agglomeration method used for hierarchical clustering in 'flashClust'-package. Either "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
  }
  \item{dendroRowText}{
boolean. If TRUE, phenotype names are plotted beneath the sample dendrogram.
}
  \item{autoColorHeight}{
boolean. If TRUE, the height of the color area below the dendrogram is adjusted automatically for the number of phenotypes.
  }
  \item{colorHeight}{
numeric specifying the height of the color area under dendrogram as a fraction of the height of the dendrogram area. Only effective when autoColorHeight above is FALSE.
  }
  \item{cex.dendroLabels}{
numeric with character expansion factor for dendrogram (sample) labels.
}
  \item{\dots}{
further arguments to be passed to blockwiseModules()-function of 'WGCNA'-package.}
}
\details{
If no value is given in 'softThresholdPower', the function analyses scale free topology for multiple soft thresholding powers to help choosing the appropriate value for obtaining an approximately scale free network topology. For each power the scale free topology fit index is calculated and returned along with other information on connectivity. If 'softThresholdPower' is set to 'auto' and the function finds an appropriate value, it is directly used for network construction. Network construction is performed in block-wise manner with respect to 'maxBlockSize'. Genes are clustered using average linkage hierarchical clustering and coexpressed gene modules are identified in the resulting dendrogram by the Dynamic Hybrid tree cut. Modules whose module eigengenes (MEs) are highly correlated are merged.
The function calclulates the following parameter:

- kME: INTRAmodular connectivity for finding intramodular hubs. Also known as module membership measure (MM). Correlation of the gene with the corresponding module eigengene. kME close to 1 means that the gene is a hub gene.

- GS: gene significance: correlation of the gene with a phenotype.

- Module-trait relationship: correlation of a module eigengene with a phenotype.

Phenotypes are taken from phenotype data of 'GEXMTSet' as specified in 'phModule'. Furthermore, membership of samples in groups which are defined in 'groupsets' are also used as phenotypes (e.g. two groups from a differential gene expression experiment). When correlation with group membership is calculated, only those samples are included which belong to the denoted groupset (mind that gene modules were calculated using expression data from all samples). All correlation coefficients are calculated using pearson correlation for numeric phenotypes and intraclass correlation (ICC) for categorical phenotypes (e.g. groupsets with more than two groups). Categorical variables with only two levels are coded numerically.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

value depends on parameter 'softThresholdPower'. If a softThresholdPower is given or could be chosen automatically, value is a list with the following components:
\item{colors }{a vector of color or numeric module labels for all genes.}
\item{unmergedColors }{a vector of color or numeric module labels for all genes before module merging.}
\item{MEs }{a data frame containing module eigengenes of the found modules (given by colors).}
\item{goodSamples }{numeric vector giving indices of good samples, that is samples that do not have too many missing entries.}
\item{goodGenes }{numeric vector giving indices of good genes, that is genes that do not have too many missing entries.}
\item{dendrograms }{a list whose components conatain hierarchical clustering dendrograms of genes in each block.}
\item{TOMFiles }{character vector (one string per block), giving the file names in which blockwise topological overlaps were saved.}
\item{blockGenes }{a list whose components give the indices of genes in each block.}
\item{blocks }{a vector of length equal number of genes giving the block label for each gene. Note that block labels are not necessarilly sorted in the order in which the blocks were processed }
\item{MEsOK }{logical indicating whether the module eigengenes were calculated without errors.}

If 'softThresholdPower' is NULL or could not be chosen automatically, value is a list with the following components:
\item{powerEstimate }{estimate of an appropriate soft-thresholding power: the lowest power for which the scale free topology fit R^2 exceeds RsquaredCut. If R^2 is below RsquaredCut for all powers, NA is returned.}
\item{fitIndices }{data frame containing the fit indices for scale free topology. The columns contain the soft-thresholding power, adjusted R^2 for the linear fit, the linear coefficient, adjusted R^2 for a more complicated fit models, mean connectivity, median connectivity and maximum connectivity.}

Side-effects:
Diagrams for SoftThreshold power, gene and sample dendrograms generated by hierarchical clustering with phenotypes given in 'phDendro' or 'phModule' printed underneath as well as correlation heatmaps are plotted into the projectfolder. Additionally, tables with module eigengenes and correlation results of eigengenes with phenotypes and groupsets are generated. Scatterplots are generated with module membership and gene significance for each phenotype/groupset and the 8 top associated modules.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Frank Ruehle}
\note{
The procedure is divided in several steps:
  Step1: Selection of an appropriate softThresholdPower for network construction.

  Step2: Automatic network construction and module detection.

  Step3: Plot sample dendrogram and gene dendrogramm with phenotype information.
      Calculate gene significance for traits: GS.datTraits(i) = |cor(gene,Trait)|.
      GSPvalue[i] = corPvalueStudent(GS.datTraits[i], nSamples).

  Step4: Correlation of modules with phenotypes (traits).
      moduleTraitCor =  cor(MEs, Trait).
      moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples).
      make 'groupsetMat' as phenotype matrix from group memberships for groups denoted in 'groupsets'.
      moduleGroupsetCor = cor(MEs, groupsetMat).
      moduleGroupsetPvalue = corPvalueStudent(moduleGroupsetCor, nSamplesInGroups).
      Heatmaps are generated for correlation results of phenotypes and groupsets.

  Step5: Intramodular analysis - Find hub genes in modules.
      datKME: INTRAmodular connectivity for finding intramodular hubs. Also known as module membership measure (MM).
      MMPvalue = corPvalueStudent(datKME, nSamples).
      Calculage gene significance for group memberships: geneGroupsetCor = cor(gene,groupsetMat).
      geneGroupsetPvalue = corPvalueStudent(geneGroupsetCor, nSamplesInGroups).
      Scatterplots are generated with module membership and gene significance for each phenotype/groupset and the 8 top associated modules.

  Step6: Generate output tables.
      networkDatOutput0 = data.frame(featureGEXMTSet, moduleColors, GS.datTraits, GSPvalue, geneGroupsetCor, geneGroupsetPvalue).
      networkDatOutputMM = data.frame(networkDatOutput0, datKME[,modOrder], MMPvalue[,modOrder]).

  Step7: Visualization of networks within R.
      TOMplot.
      MDSplot.
  }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
