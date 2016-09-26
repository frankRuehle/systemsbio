\name{filterGeneLists}
\alias{filterGeneLists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter gene table for specified categories
}
\description{
Filters dataframe or character directing to dataframe for up to two categories (e.g. p-value and foldchange). Columns for filtercriteria must be included. Values in filtercriteria may be transformed e.g. to ABSOLUTE values by filtercat.function = abs, which is needed for log foldchanges. No tranformation if filtercat.function == identity".
}
\usage{
filterGeneLists(genes, newheader = NULL,
  filtercat1 = "adj.P.Val", filtercat1.decreasing = FALSE,
  filtercat1.function = abs, filtercat1.threshold = 0.05,
  filtercat2 = "logFC", filtercat2.decreasing = TRUE,
  filtercat2.function = abs,
  filtercat2.threshold = log2(1.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
dataframe or character with path directing to table with genelist and columns to filter for.
}
  \item{newheader}{
NULL if 'genes' already supplied with header. Character vector with new header otherwise.
}
  \item{filtercat1}{
column name of first category to filter 'genes' (e.g. p-values). Skipped if NULL.
}
  \item{filtercat1.decreasing}{
(boolean) direction to order and filter filtercat1
 }
  \item{filtercat1.function}{
select transforming function for filter category1 (no quotes). e.g. abs for absolute values, identity for no transformation
}
  \item{filtercat1.threshold}{
Threshold for filtercat1 or 'top123' for top Hits
}
  \item{filtercat2}{
column name of second category to filter 'genes' (e.g. effect size). Skipped if NULL.
}
  \item{filtercat2.decreasing}{
(boolean) direction to order and filter filtercat2
}
  \item{filtercat2.function}{
select transforming function for filter category2 (no quotes). E.g. abs for foldchanges
}
  \item{filtercat2.threshold}{
Threshold for filtercat2 or 'top123' for top Hits
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
dataframe filtered for desired criteria and sorted for filtercat1.
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
