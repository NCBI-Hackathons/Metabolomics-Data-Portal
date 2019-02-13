\name{reach}
\alias{reach}
\title{
Calculate largest reach centrality
}
\description{
Calculate largest reach centrality
}
\usage{
reach(graph, weights=E(graph)$weight, mode=c("all", "in", "out"))
}
\arguments{

  \item{graph}{an \code{\link[igraph]{igraph}} object}
  \item{mode}{mode of the centrality}
  \item{weights}{If the edges in the graph have weight, then by default, the weight is used to calculate the length of the shortest path. Set it to NULL to supress the weight.}

}
\details{
The largest reach centrality measures how far a node can send or receive the information in the network.
It is defined as the largest length of the shortest path from all the other nodes in the network.
}
\examples{
# There is no example
NULL

}
