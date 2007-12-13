\name{connectivity}
\alias{connectivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Connectivity Measure}
\description{
  Calculates the connectivity validation measure for a given cluster partitioning.
}
\usage{
connectivity(distance = NULL, clusters, Data = NULL, neighbSize = 10, method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{The distance matrix (as a matrix object) of the
    clustered observations.  Required if \code{Data} is NULL.}
  \item{clusters}{An integer vector indicating the cluster partitioning}
  \item{Data}{The data matrix of the clustered observations. Required if
    \code{distance} is NULL.}  
  \item{neighbSize}{The size of the neighborhood}
  \item{method}{The metric used to determine the distance
    matrix.  Not used if \code{distance} is provided.}
}
\details{
  The connectivity indicates the degree of connectedness of the
  clusters, as determined by the k-nearest neighbors.  The
  \code{neighbSize} argument specifies the number of neighbors to use.
  The connectivity has a value between 0 and infinity and should be minimized.
  For details see the package vignette.
}
\value{
  Returns the connectivity measure as a numeric value.
}

\references{

  Handl, J., Knowles, K., and Kell, D. (2005).
  Computational cluster validation in post-genomic data analysis.
  Bioinformatics 21(15): 3201-3212.

}  

\author{Guy Brock, Vasyl Pihur, Susmita Datta, Somnath Datta}
\note{
  The main function for cluster validation is \code{\link{clValid}}, and
  users should call this function directly if possible.
}
\seealso{
  For a description of the function 'clValid' see \code{\link{clValid}}.
  
  For a description of the class 'clValid' and all available methods see
  \code{\link{clValidObj}} or \code{\link{clValid-class}}.

  For additional help on the other validation measures see
  \code{\link{dunn}},
  \code{\link{stability}}, 
  \code{\link{BHI}}, and
  \code{\link{BSI}}.
}


\examples{
data(mouse)
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
## hierarchical clustering
Dist <- dist(express,method="euclidean")
clusterObj <- hclust(Dist, method="average")
nc <- 2 ## number of clusters      
cluster <- cutree(clusterObj,nc)
connectivity(Dist, cluster)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
