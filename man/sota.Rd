\name{sota}
\alias{sota}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Self-organizing Tree Algorithm (SOTA)}
\description{
  Computes a Self-organizing Tree Algorithm (SOTA) clustering of a dataset returning a SOTA object.
}
\usage{
sota(data, maxCycles, maxEpochs = 1000, distance = "euclidean", wcell = 0.01, 
     pcell = 0.005, scell = 0.001, delta = 1e-04, neighb.level = 0, 
     maxDiversity = 0.9, unrest.growth = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix or data frame. Cannot have a profile ID as the first column.}
  \item{maxCycles}{integer value representing the maximum number of iterations allowed. The resulting number
  of clusters returned by \code{sota} is maxCycles+1 unless \code{unrest.growth} is set to FALSE and the
  \code{maxDiversity} criteria is satisfied prior to reaching the maximum number of iterations.}
  \item{maxEpochs}{integer value indicating the maximum number of training epochs allowed per cycle. By default,
  \code{maxEpochs} is set to 1000.}
  \item{distance}{character string used to represent the metric to be used for calculating
  dissimilarities between profiles. 'euclidean' is the default, with 'correlation' being another option.}
  \item{wcell}{value specifying the winning cell migration weight. The default is 0.01.}
  \item{pcell}{value specifying the parent cell migration weight. The default is 0.005.}
  \item{scell}{value specifying the sister cell migration weight. The default is 0.001.}
  \item{delta}{value specifying the minimum epoch error improvement. This value is used as a threshold for signaling
  the start of a new cycle. It is set to 1e-04 by default.}
  \item{neighb.level}{integer value used to indicate which cells are candidates to accept new profiles. This number
  specifies the number of levels up the tree the algorithm moves in the search of candidate cells for the redistribution
  of profiles. The default is 0.}
  \item{maxDiversity}{value representing a maximum variability allowed within a cluster. 0.9 is the default value.}
  \item{unrest.growth}{logical flag: if TRUE then the algorithm will run \code{maxCycles} iterations regardless of
  whether the \code{maxDiversity} criteria is satisfied or not and \code{maxCycles}+1 clusters will be produced; if FALSE
  then the algorithm can potentially stop before reaching the \code{maxCycles} based on the current state of cluster
  diversities. A smaller than usual number of clusters will be obtained. The default value is TRUE.}
  \item{\dots}{Any other arguments.}
}

\details{
The Self-Organizing Tree Algorithm (SOTA) is an unsupervised neural network with a binary tree topology. It combines
the advantages of both hierarchical clustering and Self-Organizing Maps (SOM). The algorithm picks a node with
the largest Diversity and splits it into two nodes, called Cells. This process can be stopped at any level, assuring a fixed number of
hard clusters. This behavior is achieved with setting the \code{unrest.growth} parameter to TRUE. Growth of the
tree can be stopped based on other criteria, like the allowed maximum Diversity within the cluster and so on.

Further details regarding the inner workings of the algorithm can be found in the paper listed in the Reference
section.
}
\value{
  \item{data }{data matrix used for clustering}
  \item{c.tree }{complete tree in a matrix format. Node ID, its Ancestor, and whether it's a terminal node (cell)
  are listed in the first three columns. Node profiles are shown in the remaining columns.}
  \item{tree }{incomplete tree in a matrix format listing only the terminal nodes (cells).
  Node ID, its Ancestor, and 1's for a cell indicator
  are listed in the first three columns. Node profiles are shown in the remaining columns.}
  \item{clust }{integer vector whose length is equal to the number of profiles in a data matrix indicating
  the cluster assingments for each profile in the original order.}
  \item{totals }{integer vector specifying the cluster sizes.}
  \item{dist }{character string indicating a distance function used in the clustering process.}
  \item{diversity }{vector specifying final cluster diverisities.}
}
\references{Herrero, J., Valencia,
A, and Dopazo, J. (2005). A hierarchical unsupervised growing neural
network for clustering gene expression patterns. Bioinformatics, 17, 126-136.}

\author{Vasyl Pihur, Guy Brock, Susmita Datta, Somnath Datta}

\seealso{\code{\link{plot.sota}}, \code{\link{print.sota}} }
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