% Generated by roxygen2 (4.0.2.9000): do not edit by hand
\name{hier_clust}
\alias{hier_clust}
\title{hier_clust}
\usage{
hier_clust(tdm, clusters = 5, sparsity = 0.98, method = "ward.D")
}
\arguments{
\item{tdm}{Required. Term-document matrix.}

\item{clusters}{Optional. Number of clusters based on distance between terms, defaults to 5.}

\item{sparsity}{Optional. Maximum allowed sparsity. Between 0 and 1, defaults to 0.98.}

\item{method}{Optional. Clustering method, defaults to ward.D (see details) - other: single linkage, complete linkage, average linkage, median and centroid}
}
\value{
Clusters, plot fit.
}
\description{
Hierarchical cluster analysis presented in dendrogram.
}
\details{
Transforms tdm, removes sparse terms and clusters words. Default method is ward.D which clusters terms according to the varianced.
}
\examples{
\dontrun{

library(tm)

data(crude)

tdm <- TermDocumentMatrix(crude, control=list(minWordLength=1))

fit <- hier_clust(tdm, sparsity = 0.97)

plot(fit)

}
}
\seealso{
\code{\link{pam_k}} to cluster around medoids
\code{\link[tm]{TermDocumentMatrix}} TermDocumentMatrix from corpus
\code{\link[stats]{hclust}} Clustering methods
}

