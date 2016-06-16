\name{basicAnno}
\alias{basicAnno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate gene symbols and EntrezIDs using an organism-specific annotation package.
}
\description{
If either gene symbols or EntrezIDs are missing in 'data', this annotation is added using the corresponding annotation package of the organism specified in 'org'.
  }
\usage{
basicAnno(data, Symbol.column = "SYMBOL", Entrez.column = "ENTREZID", org = "human")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
dataframe or character with filepath to data to be loaded.
}
  \item{Symbol.column}{
character with column name of Gene Symbols in 'data'. If NULL, SYMBOLs are derived from 'Entrez.column' using the annotation package for 'org'.
  }
  \item{Entrez.column}{
character with column name of ENTREZ IDs in 'data'. If NULL, ENTREZ IDs are derived from 'Symbol.column' using the annotation package for 'org'. Either 'Symbol.column' or 'Entrez.column' must be specified.
  }
  \item{org}{
character with species name ("human", "mouse", "rat").
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
dataframe with added annotation column.
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
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
