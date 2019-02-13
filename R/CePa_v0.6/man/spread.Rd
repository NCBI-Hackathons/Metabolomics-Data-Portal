\name{spread}
\alias{spread}
\title{
Calculate radiality centrality
}
\description{
Calculate radiality centrality
}
\usage{
spread(graph, mode = c("all", "in", "out"),
    weights = E(graph)$weight, f = function(x) 1/x)
}
\arguments{

  \item{graph}{an \code{\link[igraph]{igraph}} object}
  \item{mode}{mode of the centrality}
  \item{weights}{If edges in the graph have weight, then by default, the weight is used to calculate the length of the shortest path. Set it to NULL to supress the weight}
  \item{f}{function for the weaken rate}

}
\details{
The spread centrality measures how wide the node can send or receive the information in the network.
Like the water wave, the effect would be weakened with the increase of the distance to other nodes.
}
\examples{
# There is no example
NULL

}
