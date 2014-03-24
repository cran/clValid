\name{annotationListToMatrix}
\alias{annotationListToMatrix}
\title{ Change annotation list to matrix}
\description{
  Change biological functional annotation from list to TRUE / FALSE matrix
}
\usage{
annotationListToMatrix(annotation, genenames)
}
\arguments{
  \item{annotation}{ functional annotation of genes, as a list}
  \item{genenames}{ vector of genenames}
}
\details{
  Converts biological annotation from ``list'' to ``matrix'' format.  In
  ``list'' format, each item in the list is a vector giving genes
  belonging to a particular biological class.  In ``matrix'' format,
  each column is a logical vector indicating which genes belong to the
  biological class.  Both \code{\link{BHI}} and
  \code{\link{BSI}} use the TRUE / FALSE matrix format for inputting
  predetermined biological functional classes.
}
\value{
  A logical matrix, where each element indicates whether the gene in row
  i belongs to the biological functional class of column j
}
\author{Guy Brock}
\note{Special thanks to Rainer Machne, who initially suggested this change.}

\seealso{\code{\link{BHI}}, \code{\link{BSI}}}
\examples{

data(mouse)
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
fc <- tapply(rownames(express),mouse$FC[1:25], c)
fc <- fc[-match( c("EST","Unknown"), names(fc))]
fc <- annotationListToMatrix(fc, rownames(express))

## see package vignette for example use when reading from Excel file

}

\keyword{cluster}

