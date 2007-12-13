\name{BHI}
\alias{BHI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biological Homogeneity Index}
\description{
  Calculates the biological homogeneity index (BHI) for a given
  statistical clustering partition and biological annotation.
}
\usage{
BHI(statClust, annotation, names = NULL, category = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statClust}{An integer vector indicating the statistical cluster partitioning}
  \item{annotation}{Either a character string naming the Bioconductor
    annotation package for mapping genes to GO categories, or a list
    with the names of the functional classes
    and the observations belonging to each class}
  \item{names}{A vector of labels to associate with the 'genes', to be
    used in conjunction with the Bioconductor annotation package.  Not
    needed if \code{annotation} is a list providing the functional classes. }
  \item{category}{Indicates the GO 
    categories to use for biological validation.  Can be one of "BP",
    "MF", "CC", or "all".}
}
\details{
  The BHI measures how homogeneous the clusters are biologically.  The
  measure checks whether genes placed
  in the same statistical cluster also belong to the same functional
  classes.  The BHI is in the range [0,1], with
  larger values corresponding to more biologically homogeneous clusters.
  For details see the package vignette.
}
\value{
  Returns the BHI measure as a numeric value.
}
\references{

  Datta, S. and Datta, S. (2006).
  Methods for evaluating clustering algorithms for gene expression data
  using a reference set of functional classes.
  BMC Bioinformatics 7:397.
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
  \code{\link{connectivity}},   \code{\link{dunn}},
  \code{\link{stability}}, and
  \code{\link{BSI}}.
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

## first way - functional classes predetermined
fc <- tapply(rownames(express),mouse$FC[1:25], c)
fc <- fc[-match( c("EST","Unknown"), names(fc))]
BHI(cluster, fc)

## second way - using Bioconductor
if(require("Biobase") && require("annotate") && require("GO") &&
require("moe430a")) {
  BHI(cluster, annotation="moe430a", names=rownames(express), category="all")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

