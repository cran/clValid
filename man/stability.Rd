\name{stability}
\alias{stability}
\alias{APN}
\alias{AD}
\alias{ADM}
\alias{FOM}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stability Measures}
\description{
  Calculates the stability measures the average proportion of non-overlap (APN),
  the average distance (AD), the average distance between means (ADM),
  and the figure of merit (FOM).
}
\usage{
stability(mat, Dist=NULL, del, cluster, clusterDel, method="euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{The data matrix of the clustered observations}  
  \item{Dist}{The distance matrix (as a matrix or dist object) of the
    clustered observations.  If NULL then \code{method} is used with
    \code{mat} to determine the distance matrix.}
  \item{del}{An integer indicating which column was removed}
  \item{cluster}{An integer vector indicating the cluster partitioning
    based on all the data}  
  \item{clusterDel}{An integer vector indicating the cluster
    partitioning based on the data with column \code{del} removed.}
  \item{method}{The metric used to determine the distance
    matrix.  Not used if \code{distance} is provided.}  
}
\details{
  The stability measures evaluate the stability of a clustering result by comparing it
  with the clusters obtained by removing one column at a time.
  These measures include the average proportion of non-overlap (APN),
  the average distance (AD), the average distance between means (ADM),
  and the figure of merit (FOM).  The APN, AD, and ADM are all based on the
  cross-classification table of the original clustering with the
  clustering based on the removal of one column.  The APN measures the
  average proportion of observations not placed in the same cluster under both
  cases, while the AD measures the average distance between observations placed
  in the same cluster under both cases and the ADM measures the average
  distance between cluster centers for observations placed in the same cluster
  under both cases.  The FOM measures the average intra-cluster variance
  of the deleted column, where the clustering is based on the remaining
  (undeleted) columns.  In all cases the average is taken over all the
  deleted columns, and all measures should be minimized.
  For details see the package vignette.

  NOTE: The \code{stability} function only calculates these measures for
  the particular column specified by \code{del} removed.  To get the
  overall scores, the user must average the measures
  corresponding to each removed column.
}
\value{
  Returns a numeric vector with the APN, AD, ADM, and FOM measures
  corresponding to the particular column that was
  removed.  
}
\references{
  Datta, S. and Datta, S. (2003).
  Comparisons and validation of statistical clustering techniques for microarray gene expression data.
  Bioinformatics 19(4): 459-466.  
}  

\author{Guy Brock, Vasyl Pihur, Susmita Datta, Somnath Datta}
\note{
  The main function for cluster validation is \code{\link{clValid}}, and
  users should call this function directly if possible.

  To get the overall values, the stability measures
  corresponding to each removed column should be averaged (see the
  examples below).  
}
\seealso{
  For a description of the function 'clValid' see \code{\link{clValid}}.
  
  For a description of the class 'clValid' and all available methods see
  \code{\link{clValidObj}} or \code{\link{clValid-class}}.

  For additional help on the other validation measures see
  \code{\link{connectivity}},   \code{\link{dunn}},
  \code{\link{BSI}}, and 
  \code{\link{BHI}}.
}

\examples{

data(mouse)
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
## hierarchical clustering
Dist <- dist(express,method="euclidean")
clusterObj <- hclust(Dist, method="average")
nc <- 4 ## number of clusters      
cluster <- cutree(clusterObj,nc)

stab <- matrix(0,nrow=ncol(express),ncol=4)
colnames(stab) <- c("APN","AD","ADM","FOM")

## Need loop over all removed samples
for (del in 1:ncol(express)) {
  matDel <- express[,-del]               
  DistDel <- dist(matDel,method="euclidean")
  clusterObjDel <- hclust(DistDel, method="average")
  clusterDel <- cutree(clusterObjDel,nc)
  stab[del,] <- stability(express, Dist, del, cluster, clusterDel)
}
colMeans(stab)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
