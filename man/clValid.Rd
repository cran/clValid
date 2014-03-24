\name{clValid}
\alias{clValid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Validate Cluster Results}
\description{
  \code{clValid} reports validation measures for clustering
  results.  The function returns an object of class
  \code{"\linkS4class{clValid}"}, which
  contains the clustering results in addition to the validation
  measures.  The validation measures fall into three general categories:
  "internal", "stability", and "biological".
}
\usage{
clValid(obj, nClust, clMethods = "hierarchical", validation =
"stability", maxitems = 600, metric = "euclidean", method = "average",
neighbSize = 10, annotation = NULL, GOcategory = "all",
goTermFreq=0.05, dropEvidence=NULL, verbose=FALSE, ...)
}


\arguments{
  \item{obj}{Either a numeric matrix, a data frame, or an \code{\link[Biobase]{ExpressionSet}}
    object.  Data frames must contain all numeric columns.  In all
    cases, the rows are the items to be clustered (e.g., genes),
    and the columns are the samples.}
  \item{nClust}{A numeric vector giving the numbers of clusters
    to be evaluated.  e.g., 4:6 would evaluate the number of clusters
    ranging from 4 to 6.}
  \item{clMethods}{A character vector giving the clustering methods.
    Available options are "hierarchical", "kmeans", "diana", "fanny",
    "som", "model", "sota", "pam", "clara", and "agnes", with multiple
    choices allowed.}
  \item{validation}{A character vector giving the type of
    validation measures to use.  Available options are "internal",
    "stability", and "biological", with multiple choices allowed.  }
  \item{maxitems}{The maximum number of items (rows in matrix) which can be clustered.}
  \item{metric}{The metric used to determine the distance
    matrix.  Possible choices are "euclidean", "correlation", and "manhattan".}
  \item{method}{For hierarchical clustering (\code{hclust} and \code{agnes}), the
    agglomeration method used.  Available choices are "ward", "single",
    "complete", and "average".}
  \item{neighbSize}{For internal validation, an integer giving the neighborhood size used for the
    connectivity measure.}
  \item{annotation}{For biological validation,
    either a character string naming the Bioconductor annotation
    package for mapping genes to GO categories, or a list with the names of the functional classes
    and the observations belonging to each class.}
  \item{GOcategory}{For biological validation, gives which GO
    categories to use for biological validation.  Can be one of "BP",
    "MF", "CC", or "all".}
  \item{goTermFreq}{For the
    BSI, what threshold frequency of GO terms to use for functional
    annotation.}
  \item{dropEvidence}{Which GO evidence codes to omit.  Either NULL or a
    character vector, see 'Details' below.}
  \item{verbose}{Logical - if TRUE will produce detailed output on the
    progress of cluster validation.}
  \item{\dots}{Additional arguments to pass to the clustering functions.}
}


\details{

  This function calculates validation measures for a given set of
  clustering algorithms and number of clusters.  A variety of clustering
  algorithms are available, including hierarchical, self-organizing maps
  (SOM), K-means, self-organizing tree algorithm (SOTA), and model-based.
  The available validation measures fall into the three general
  categories of "internal", "stability", and "biological".  A brief
  description of each measure is given below, for further details refer
  to the package vignette and the references.

  \describe{
    \item{\bold{Internal measures:}}{The internal
  measures include the connectivity, and Silhouette Width, and Dunn
  Index.  The connectivity indicates the degree of connectedness of the
  clusters, as determined by the k-nearest neighbors.  The
  \code{neighbSize} argument specifies the number of neighbors to use.
  The connectivity has a value between 0 and infinity and should be minimized.
  Both the Silhouette Width and the Dunn Index combine measures of
  compactness and separation of the clusters.  The Silhouette Width is
  the average of each observation's Silhouette value.  The Silhouette
  value measures the degree of confidence in a particular clustering
  assignment and lies
  in the interval [-1,1], with well-clustered observations having values
  near 1 and poorly clustered observations having values near -1.  See
  the \code{\link[cluster]{silhouette}} function in package \pkg{cluster} for
  more details.  The
  Dunn Index is the ratio between the smallest distance between
  observations not in the same cluster to the largest intra-cluster
  distance.  It has a value between 0 and infinity and should be maximized.}

\item{\bold{Stability measures:}}{The stability measures are a special version of internal measures
  which evaluate the stability of a clustering result by comparing it
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
  deleted columns, and all measures should be minimized.  }


\item{\bold{Biological measures:}}{There are two biological validation measures, the biological homogeneity index (BHI) and
  biological stability index (BSI).  The observations are typically taken to
  represent a `gene' (e.g., ORF, SAGE
  tag, affy ID).  The BHI measures the average
  proportion of gene pairs that are clustered together which have
  matching biological functional classes.
  The BSI is similar to the other stability
  measures, but 
  inspects the consistency of
  clustering for genes with similar biological functionality. Each sample is removed one
  at a time, and the cluster membership for genes with similar functional annotation is
  compared with the cluster membership using all available samples.

  For biological validation, the user has two options.
  The first option is to explicity specify the
  functional clustering of the genes via either a named list or logical
  matrix. In ``list'' format, each item in the list is a vector giving genes
  belonging to a particular biological class.  In ``matrix'' format,
  each column is a logical vector indicating which genes belong to the
  biological class. \code{clValid} will convert the biological
  annotation to matrix format internally if initially given in list format.

  The second option is to specify the appropriate
  annotation package from Bioconductor (\url{http://www.bioconductor.org})
  and GO terms to determine the functional classes of
  the genes.  To use the second option requires the \pkg{Biobase},
  \pkg{annotate}, and \pkg{GO} packages from Bioconductor, in addition
  to the annotation package for the particular data type.
  If the annotation package cannot be loaded, \code{clValid} will
  attempt to automatically download the package from
  \url{www.bioconductor.org} (using the
  \file{biocLite.R} installation script).

  The \code{GOcategory} options are "MF", "BP", "CC", or "all",
  corresponding to molecular function, biological process, cellular
  component, and all of the ontologies.

  The \code{dropEvidence} argument indicates which GO evidence codes to
  omit.  For example, "IEA" is a relatively weak association based only
  on electronic information, and users may wish to omit this evidence
  when determining the functional annotation classes.
  
}}



}
\value{

  \code{clValid} returns an object of class
  \code{"\linkS4class{clValid}"}.  See the help file for the class description.

}
\references{

  Brock, G., Pihur, V., Datta, S. and Datta, S. (2008)
  clValid: An R Package for Cluster Validation
  Journal of Statistical Software 25(4)
  \url{http://www.jstatsoft.org/v25/i04}
  
  Datta, S. and Datta, S. (2003)
  Comparisons and validation of statistical clustering techniques for microarray gene expression data.
  Bioinformatics 19(4): 459-466

  Datta, S. and Datta, S. (2006)
  Methods for evaluating clustering algorithms for gene expression data
  using a reference set of functional classes.
  BMC Bioinformatics 7:397
  \url{http://www.biomedcentral.com/1471-2105/7/397}

  Handl, J., Knowles, K., and Kell, D. (2005)
  Computational cluster validation in post-genomic data analysis.
  Bioinformatics 21(15): 3201-3212
  
}

\author{Guy Brock, Vasyl Pihur, Susmita Datta, Somnath Datta}
\note{
The only package which is automatically attached is \pkg{cluster}.  To
use the clustering methods \code{som} and \code{Mclust} you will need
to load the packages \pkg{kohonen} and \pkg{mclust}, respectively.

Unless the the list of genes corresponding to functional classes is prespecified,
to perform biological clustering validation
will require the \pkg{Biobase}, \pkg{annotate} and \pkg{GO} packages from
Bioconductor, and in addition the annotation package for your
particular data type.  Please see \url{http://www.bioconductor.org} for installation instructions.

Further details of the validation measures and instructions in use can
be found in the package vignette.
}

\seealso{
  For a description of the class 'clValid' and all available methods see \code{\link{clValidObj}}
  or \code{\link{clValid-class}}.

  For help on the clustering methods see \code{\link{hclust}} and
  \code{\link{kmeans}}  in package \pkg{stats}, 
  \code{\link[cluster]{agnes}}, \code{\link[cluster]{clara}},  \code{\link[cluster]{diana}},
  \code{\link[cluster]{fanny}}, and \code{\link[cluster]{pam}} in package \pkg{cluster},
  \code{\link[kohonen]{som}} in package \pkg{kohonen}, \code{\link[mclust]{Mclust}}
  in package \pkg{mclust}, and \code{\link{sota}} (in this package).

  For additional help on the validation measures see
  \code{\link{connectivity}},   \code{\link{dunn}},
  \code{\link{stability}}, 
  \code{\link{BHI}}, and 
  \code{\link{BSI}}. 

  % links to other validation packages??
  
}

\examples{


data(mouse)

## internal validation
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
intern <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                  validation="internal")

## view results
summary(intern)
optimalScores(intern)
plot(intern)

## stability measures
stab <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                validation="stability")
optimalScores(stab)
plot(stab)

## biological measures
## first way - functional classes predetermined
fc <- tapply(rownames(express),mouse$FC[1:25], c)
fc <- fc[-match( c("EST","Unknown"), names(fc))]
bio <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
               validation="biological", annotation=fc)
optimalScores(bio)
plot(bio)

## second way - using Bioconductor
if(require("Biobase") && require("annotate") && require("GO.db") && require("moe430a.db")) {
  bio2 <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                  validation="biological",annotation="moe430a.db",GOcategory="all")
  optimalScores(bio2)
  plot(bio2)
}

}

  

\keyword{cluster}
