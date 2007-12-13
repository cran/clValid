\name{BSI}
\alias{BSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biological Stability Index}
\description{
  Calculates the biological stability index (BSI) for a given
  statistical clustering partition and biological annotation.
}

\usage{
BSI(statClust, statClustDel, annotation, names = NULL, category = "all", goTermFreq = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statClust}{An integer vector indicating the statistical cluster partitioning}  
  \item{statClustDel}{An integer vector indicating the statistical cluster
    partitioning based on one column removed}
  \item{annotation}{Either a character string naming the Bioconductor
    annotation package for mapping genes to GO categories, or a list
    with the names of the functional classes
    and the observations belonging to each class.}
  \item{names}{An optional vector of names for the observations}
  \item{category}{Indicates the GO 
    categories to use for biological validation.  Can be one of "BP",
    "MF", "CC", or "all".}
  \item{goTermFreq}{What threshold frequency of GO terms to use for
    functional annotation.}  
}
\details{
  The BSI  inspects the consistency of clustering for genes with
  similar biological functionality.  Each sample is removed, 
  and the cluster membership for genes with similar functional annotation
  is compared with the cluster membership using all available samples.
  The BSI is in the range [0,1], with larger values corresponding to more stable
  clusters of the functionally annotated genes.
  For details see the package vignette.

  NOTE: The \code{BSI} function only calculates these measures for
  one particular column removed.  To get the
  overall scores, the user must average the measures
  corresponding to each removed column.
}
\value{
  Returns the BSI value corresponding to the particular column that was
  removed.  
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

  To get the overall BSI value, the BSI values
  corresponding to each removed column should be averaged (see the
  examples below).
}
\seealso{
  For a description of the function 'clValid' see \code{\link{clValid}}.
  
  For a description of the class 'clValid' and all available methods see
  \code{\link{clValidObj}} or \code{\link{clValid-class}}.

  For additional help on the other validation measures see
  \code{\link{connectivity}},   \code{\link{dunn}},
  \code{\link{stability}}, and
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

## first way - functional classes predetermined
fc <- tapply(rownames(express),mouse$FC[1:25], c)
fc <- fc[-match( c("EST","Unknown"), names(fc))]
bsi <- numeric(ncol(express))
## Need loop over all removed samples
for (del in 1:ncol(express)) {
  matDel <- express[,-del]               
  DistDel <- dist(matDel,method="euclidean")
  clusterObjDel <- hclust(DistDel, method="average")
  clusterDel <- cutree(clusterObjDel,nc)
  bsi[del] <- BSI(cluster, clusterDel, fc)
}
mean(bsi)

## second way - using Bioconductor
if(require("Biobase") && require("annotate") && require("GO") &&
   require("moe430a")) {
  bsi <- numeric(ncol(express))
  for (del in 1:ncol(express)) {
    matDel <- express[,-del]               
    DistDel <- dist(matDel,method="euclidean")
    clusterObjDel <- hclust(DistDel, method="average")
    clusterDel <- cutree(clusterObjDel,nc)
    bsi[del] <- BSI(cluster, clusterDel, annotation="moe430a",
                    names=rownames(express), category="all")
  }
  mean(bsi)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
