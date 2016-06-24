\name{norm_Illu}
\alias{norm_Illu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalisation and filtering of Expression data using beadarray
}
\description{
Expression data is normalised and/or transformed by algorithms dedicated in 'method_norm' and/or 'transform'. Probes assigned a `Bad' or `No match' quality score after normalisiation are removed from the dataset. Lastly, the object class is changed from 'ExpressionSetIllumina' to 'ExpressionSet' to avoid later incompabilities with the 'limma'-package.}
\usage{
norm_Illu(eset, method_norm = "quantile", transform = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
ExpressionSet or ExpressionSetIllumina
}
  \item{method_norm}{
character with normalisation method. Options are "quantile", "qspline", "vsn", "rankInvariant", "median" and "none"
}
  \item{transform}{
character with data transformation method. Options are "none", "log2", "neqc", "rsn" and "vst".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ExpressionSet object
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
beadarray
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
