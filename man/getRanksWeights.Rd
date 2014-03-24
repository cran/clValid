\name{getRanksWeights}
\alias{getRanksWeights}
\title{Extract ranks and weights from clValid object}
\description{
  Creates matrix of ranks and weights from \code{\linkS4class{clValid}}
  object, to use as input for rank aggregation using
  \code{\link[RankAggreg]{RankAggreg}} in package \pkg{RankAggreg}}
\usage{
getRanksWeights(clVObj, measures = measNames(clVObj), nClust =
                nClusters(clVObj), clAlgs = clusterMethods(clVObj))
}
\arguments{
  \item{clVObj}{a clValid object}
  \item{measures}{the cluster validation measures to use for rank aggregation}
  \item{nClust}{the number of clusters to evaluate}
  \item{clAlgs}{the clustering algorithms to evaluate}
}
\details{
  This function extracts cluster validation measures from a
  \code{\linkS4class{clValid}} object, and creates a matrix of rankings
  where each row contains a list of clustering algorithms which are
  ranked according to the validation measure for that row.  The function
  also returns the cluster validation measures as a matrix of weights,
  for use with weighted rank aggregation in the function
  \code{\link[RankAggreg]{RankAggreg}}.  Any combination of validation
  measures, numbers of clusters, and clustering algorithms can be
  selected by the user.  Number of clusters and clustering algorithms
  are appended into a single name.
}
\value{
  A list with components
  \item{ranks}{Matrix with rankings for each validation measure in each row}
  \item{weights}{Matrix of weights, corresponding to the cluster
    validation measures, which are used for weighted rank aggregation}
}
\references{

  Brock, G., Pihur, V., Datta, S. and Datta, S. (2008)
  clValid: An R Package for Cluster Validation
  Journal of Statistical Software 25(4)
  \url{http://www.jstatsoft.org/v25/i04}
  
  Pihur, V., Datta, S. and Datta, S. (2009)
  RankAggreg, an R package for weighted rank aggregation
  BMC Bioinformatics 10:62  
  \url{http://www.biomedcentral.com/1471-2105/10/62}
  
}
\author{Guy Brock}
\seealso{\code{\linkS4class{clValid}}, \code{\link[RankAggreg]{RankAggreg}}}
\examples{
data(mouse)
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
clv <- clValid(express, 4:6, clMethods=c("hierarchical","kmeans","pam"), 
                  validation=c("internal","stability"))
res <- getRanksWeights(clv)
if(require("RankAggreg")) {
  CEWS <- RankAggreg(x=res$ranks, k=5, weights=res$weights, seed=123, verbose=FALSE)
  CEWS
}
}

\keyword{cluster}
