\name{radiality}
\alias{radiality}
\title{
Calculate radiality centrality
}
\description{
Calculate radiality centrality
}
\usage{
radiality(graph, mode = c("all", "in", "out"))
}
\arguments{

  \item{graph}{an \code{\link[igraph]{igraph}} object}
  \item{mode}{mode of the centrality}

}
\details{
The radiality is defined as \code{sum(d_G + 1 - d(v, w))/(n - 1)}. where \code{d(w, v)} is the
length of the shortest path from node \code{w} to node \code{v}, \code{d_G} is the diameter of the network,
n is the size of the network.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(igraph)
pathway = barabasi.game(200)
radiality(pathway)
}
