\name{clValid-class}
\docType{class}
\alias{clValidObj}
\alias{clValid-class}

\alias{clusterMethods,clValid-method}
\alias{clusters,clValid-method}
\alias{measNames,clValid-method}
\alias{measures,clValid-method}
\alias{nClusters,clValid-method}
\alias{optimalScores,clValid-method}

\alias{clusterMethods}
\alias{clusters}
\alias{measNames}
\alias{measures}
\alias{nClusters}
\alias{optimalScores}

\alias{plot,clValid,missing-method}
\alias{print,clValid-method}
\alias{show,clValid-method}
\alias{summary,clValid-method}

\title{Class "clValid"}
\description{The class \code{"clValid"} contains the clustering results
  and validation measures from the accompanying call to the function
  \code{\link{clValid}}.
}
\section{Objects from the Class}{
  Objects can be created using the function \code{\link{clValid}}.
}

\section{Slots}{
  \describe{
    \item{\code{clusterObjs}:}{Object of class \code{"list"}. A list
      containing the results from the clustering methods.}
    \item{\code{measures}:}{Object of class \code{"array"}. A
      3-dimensional array which contains the
      validation measures for the clustering results.  The first
      dimension indicates the validation measures, the second the number of
      clusters, and the third the clustering methods.}
    \item{\code{measNames}:}{Object of class \code{"character"}. The
      names of the validation measures.}
    \item{\code{clMethods}:}{Object of class \code{"character"}. A
      character vector giving the clustering methods.}
    \item{\code{labels}:}{Object of class \code{"character"}. A
      character vector giving the item (gene) labels.}    
    \item{\code{nClust}:}{Object of class \code{"numeric"}. A numeric
      vector giving the numbers of clusters
      which were evaluated.}
    \item{\code{validation}:}{Object of class \code{"character"}. A character vector giving the type of
      validation measures used, consisting of some combination of
      "internal", "stability", or "biological".}
    \item{\code{metric}:}{Object of class \code{"character"}. The metric used to determine the distance
      matrix.}
    \item{\code{method}:}{Object of class \code{"character"}. For
      hierarchical clustering, the agglomeration method used.}
    \item{\code{neighbSize}:}{Object of class \code{"numeric"}. For internal validation, the neighborhood size used for the
      connectivity measure.}
    \item{\code{annotation}:}{Object of class \code{"character or array
	or list"}.
      Either a character string naming the Bioconductor annotation
      package for mapping genes to GO categories, or a list with the names of the functional classes
      and the observations belonging to each class.}
    \item{\code{GOcategory}:}{Object of class \code{"character"}. For biological validation, gives which GO
      categories to use for biological validation.  Can be one of "BP",
      "MF", "CC", or "all"}
    \item{\code{goTermFreq}:}{Object of class \code{"numeric"}. For the
      BSI, what threshold frequency of GO terms to use for functional annotation.}
    \item{\code{call}:}{Object of class \code{"call"}.  Gives the call
      to \code{\link{clValid}} used to create the \code{clValid} object.}
  }
}

\section{Methods}{
  \describe{
    \item{clusterMethods}{\code{signature(object = "clValid")}: Returns the
      names of the clustering methods. }
    \item{clusters}{\code{signature(object = "clValid")}: Returns the
      results from the clustering methods.
      
      Additional arguments:
      \describe{
	\item{\code{method = clMethods(object)}}{The clustering method(s)
	to extract.}
      }
    }
    \item{measNames}{\code{signature(object = "clValid")}: Returns the
      names of the validation measures.}
    \item{measures}{\code{signature(object = "clValid")}: Returns
      the validation measures.
      
      Additional arguments:
      \describe{
	\item{\code{measures = measNames(object)}}{ The validation measure(s)
	to extract.}
      }
    }
    \item{nClusters}{\code{signature(object = "clValid")}: Returns the
      numbers of clusters evaluated.}
    \item{optimalScores}{\code{signature(object = "clValid")}:
      Returns the optimal value for each validation measure, along with
      the corresponding clustering method and number of clusters.
      
      Additional arguments:
      \describe{
	\item{\code{measures = measNames(object)}}{ The validation measure(s)
	to extract.}
    }}
    \item{plot}{\code{signature(x = "clValid", y = "missing")}: Plots
      the validation measures.
      
      Additional arguments:
      \describe{
	\item{\code{measures=measNames(x)}}{ The validation measures to
	plot.}
	\item{\code{legend=TRUE}}{ If TRUE provides a legend.}
	\item{\code{legendLoc="topright"}}{ The location of the legend.}
	\item{\code{main=NULL}}{ Title of graph.}
	\item{\code{pch=NULL}}{ Plotting characters to use.}
	\item{\code{type="b"}}{ Type of plot.}
	\item{\code{ask=prod(par("mfcol")) < length(measures) &&
	    dev.interactive()}}{ Logical.  If \code{TRUE} the user is
	prompted before each plot.}
      }
    }
    \item{print}{\code{signature(x = "clValid")}: Print method for class
      \code{clValid}.}
    \item{show}{\code{signature(object = "clValid")}: Same as print.}
    \item{summary}{\code{signature(object = "clValid")}: Summary method
      for class \code{clValid}.
      
      Additional arguments:
      \describe{
	\item{\code{digits = max(3,getOption("digits")-3)}}{ The number
	of significant digits to use.}
      }
    }
}}

\references{

  Brock, G., Pihur, V., Datta, S. and Datta, S. (2008)
  clValid: An R Package for Cluster Validation
  Journal of Statistical Software 25(4)
  \url{http://www.jstatsoft.org/v25/i04}
  
  Datta, S. and Datta, S. (2003)
  Comparisons and validation of statistical clustering techniques for microarray gene expression data.
  Bioinformatics 19(4): 459-466.

  Datta, S. and Datta, S. (2006)
  Methods for evaluating clustering algorithms for gene expression data
  using a reference set of functional classes.
  BMC Bioinformatics 7:397.
  \url{http://www.biomedcentral.com/1471-2105/7/397}
  
  Handl, J., Knowles, K., and Kell, D. (2005)
  Computational cluster validation in post-genomic data analysis.
  Bioinformatics 21(15): 3201-3212.

}
\author{Guy Brock, Vasyl Pihur, Susmita Datta, Somnath Datta}
\note{See the vignette for an illustration of the class.}
\seealso{
  For a description of the function 'clValid' see \code{\link{clValid}}.

  For help on the clustering methods see \code{\link{hclust}} and
  \code{\link{kmeans}} in package \pkg{stats},
  \code{\link{kmeans}}  in package \pkg{stats}, 
  \code{\link[cluster]{agnes}}, \code{\link[cluster]{clara}},  \code{\link[cluster]{diana}},
  \code{\link[cluster]{fanny}}, and \code{\link[cluster]{pam}} in package \pkg{cluster},
  \code{\link[kohonen]{som}} in package \pkg{kohonen}, \code{\link[mclust]{Mclust}}  
  in package \pkg{mclust}, and \code{\link{sota}}.

  For additional help on the validation measures see
  \code{\link{connectivity}},   \code{\link{dunn}},
  \code{\link{stability}}, 
  \code{\link{BHI}}, and 
  \code{\link{BSI}}. 

}


\examples{

## to delete
library(clValid)


data(mouse)

## internal validation
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
intern <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                  validation="internal")
slotNames(intern)

## view results
intern
summary(intern)
optimalScores(intern)
plot(intern)

## Extract objects from slots
measures(intern)
hierClust <- clusters(intern,"hierarchical")
plot(hierClust)
measNames(intern)
nClusters(intern)

}

\keyword{classes}
\keyword{cluster}
