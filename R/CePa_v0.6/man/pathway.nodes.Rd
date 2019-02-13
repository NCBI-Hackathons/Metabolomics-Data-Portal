\name{pathway.nodes}
\alias{pathway.nodes}
\title{
names of the pathway nodes
}
\description{
names of the pathway nodes
}
\usage{
pathway.nodes(pathway)
}
\arguments{

  \item{pathway}{an \code{\link[igraph]{igraph}} object}

}
\details{
If nodes in the pathway have names, then it returns a vector of nodes
names. If nodes in the pathway have no name, it just returns the index
of nodes (start from 1, after igraph version 0.6).
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
interaction = rbind(c("a", "b"),
                    c("a", "c"))
g = generate.pathway(interaction)
pathway.nodes(g)
}
