% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HLP_genelabels_plot.R
\name{determineNonOverlapGenelayers}
\alias{determineNonOverlapGenelayers}
\title{Non-overlapping gene layers for plotting}
\usage{
determineNonOverlapGenelayers(
  gr,
  ystart = 0,
  ysize = 1,
  minDistance = 0.1,
  gene_name_column = NULL,
  prot_domain_column = NULL,
  scale.factor = 1,
  units = "user",
  sortGeneStart = FALSE
)
}
\arguments{
\item{gr}{GRanges object containing gene coordinates.}

\item{ystart}{numeric start y-coordinate for first layer.}

\item{ysize}{numeric added to each new layer.}

\item{minDistance}{numeric minimum distance for two genes to be considered as non-overlapping.
If \code{minDistance} < 1, it is interpreted as ratio of the width of the full \code{gr} object.}

\item{gene_name_column}{character giving column name of gene names in \code{gr}. 
If not NULL, overlapping gene names are considered additionally to gene coordinates.}

\item{prot_domain_column}{character giving column name in \code{gr} with logical 
values indicating genes to be plotted with protein domains. For those, needed space in layer
below the respective gene layer is considered when determining non-overlapping gene layers.
Omitted if NULL.}

\item{scale.factor}{numeric scale factor as used in \code{regionalplot.genelabels}.
Required to determine character extension of gene names (if \code{gene_name_column} is not NULL).}

\item{units}{character indicating units for \code{strwidth}-function, either \code{"user"}, \code{"inches"} 
or \code{"figure"}. Measurements in \code{"user"} units are only available after \code{plot.new} has been called.
Required to determine character extension of gene names (if \code{gene_name_column} is not NULL).}

\item{sortGeneStart}{logical. If TRUE, gene levels are sorted according to minimum start coordinate.
Useful to keep free space for legend in a plot.}
}
\value{
GRanges object with column \code{ypos} added to meta data.
}
\description{
Determine gene layers for plotting of genes in regional plot without overlap
}
\details{
This is a help function for regional plots. For a given set of genes, the functions
determines separate layers necessary to plot the genes without overlapping.
}
\author{
Frank Ruehle
}
