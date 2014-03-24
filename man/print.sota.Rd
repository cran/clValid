\name{print.sota}
\alias{print.sota}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print Function for a SOTA Object}
\description{
 A default print method for a SOTA object.
}
\usage{
\method{print}{sota}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a SOTA object as returned by the \code{\link{sota}} function}
  \item{\dots}{Additional arguments to pass to \code{\link{print}}}
}
\value{
    The \code{print} function does not return anything. It simply displays in the console
    window general information about the paritioning (cluster ID, Size, and Diversity Score),
    as well as cluster centroids (average profiles within each cluster) and the distance that has
    been used.
}
\references{Herrero, J., Valencia,
A, and Dopazo, J. (2005). A hierarchical unsupervised growing neural
network for clustering gene expression patterns. Bioinformatics, 17, 126-136.}

\author{Vasyl Pihur, Guy Brock, Susmita Datta, Somnath Datta}

\seealso{\code{\link{sota}}, \code{\link{print.sota}} }
\examples{
data(mouse)
express <- mouse[,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID

sotaCl <- sota(as.matrix(express), 4)
names(sotaCl)
sotaCl
plot(sotaCl)
plot(sotaCl, cl=2)

}

\keyword{cluster}
\keyword{print}