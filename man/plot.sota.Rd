\name{plot.sota}
\alias{plot.sota}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Function for a SOTA Object}
\description{
  'plot.sota' is used to obtain a visual representation of profiles within each individual cluster.
  Corresponding cluster average profiles are also available. By default, plots for all clusters
  are displayed side by side.
}
\usage{
\method{plot}{sota}(x, cl = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{SOTA object, an object returned by function \code{\link{sota}}.}
  \item{cl}{\code{cl} specifies which cluster is to be plotted by setting it to the cluster ID. By default,
    \code{cl} is equal to 0 and the function plots all clusters side by
    side.}
  \item{\dots}{Additional arguments to pass to \code{\link{plot}}.}
  
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
\keyword{hplot}