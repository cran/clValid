\name{mouse}
\alias{mouse}
\docType{data}
\title{Mouse Mesenchymal Cells}
\description{
Data from an Affymetrix microarray experiment (moe430a) comparing comparing gene expression of mesenchymal
cells from two distinct lineages, neural crest and mesoderm derived. The
dataset consists of 147 genes and ESTs which were determined to be significantly
differentially expressed between the two cell lineages, with at least a 1.5 fold increase
or decrease in expression. There are three samples for each of the neural crest
and mesoderm derived cells.
}
\usage{data(mouse)}
\format{
  A data frame with 147 observations on the following 8 variables.
  \describe{
    \item{\code{ID}}{The Affymetric GeneChip ID, from the moe430a chip}
    \item{\code{M1}}{Mesoderm derived cell sample} 
    \item{\code{M2}}{Mesoderm derived cell sample} 
    \item{\code{M3}}{Mesoderm derived cell sample} 
    \item{\code{NC1}}{Nueral crest derived cell sample}
    \item{\code{NC2}}{Nueral crest derived cell sample}
    \item{\code{NC3}}{Nueral crest derived cell sample}
    \item{\code{FC}}{The functional class of each gene/EST}
  }
}

\source{
V. Bhattacherjee, P. Mukhopadhyay, S. Singh, C. Johnson, J. T. Philipose, C. P.
Warner, R. M. Greene, and M. M. Pisano. Neural crest and mesoderm lineagedependent
gene expression in orofacial development. Differentiation, 2007.
}

\examples{
data(mouse)

## table of fuctional classifications
table(mouse$FC)

## hierarchical clustering of expression values
express <- mouse[,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID
hc <- hclust(dist(express))
plot(hc)

}
\keyword{datasets}
