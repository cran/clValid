\name{readAnnotationFile}
\Rdversion{1.1}
\alias{readAnnotationFile}
\title{
  Read in biological annotation files from external file
}
\description{
  This function reads in a biological annotation file detailing the
  functional classes for a given set of genes from a comma separated file.
}
\usage{
readAnnotationFile(filename)
}
\arguments{
  \item{filename}{The name of the external file.  The values must be
    comma separated.  
  }
}
\details{
  The required format is comma separated, with
  the first column indicating the biological functional category, and
  the remaining columns containing the gene identifiers for those genes
  belonging to that category. 
}
\value{
  Returns a list where each item in the list is the set of genes
  belonging to a particular functional class. This can be converted to a
  TRUE/FALSE matrix using the \code{\link{annotationListToMatrix}}
  function, though it is not necessary to do this prior to using
  \code{\link{clValid}}. 
}
\author{
Guy Brock, Vasyl Pihur
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{annotationListToMatrix}}, \code{\link{clValid}}
}
\examples{
## For example use, see the package vignette
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

