\name{generate.pathway}
\alias{generate.pathway}
\title{
Generate igraph object from edge list
}
\description{
Generate igraph object from edge list
}
\usage{
generate.pathway(el)
}
\arguments{

  \item{el}{edge list, matrix with two columns. The first column is the input node and the second column is the output node.}

}
\details{
The function is a wrapper of \code{\link[igraph]{graph.edgelist}} and it generates
a directed graph.

In the function, repeated edged for two nodes will be eliminated.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{cepa}}, \code{\link[igraph]{graph.edgelist}}
}
\examples{
edgelist = rbind(c("a", "b"), c("a", "b"), c("a", "c"))
g = generate.pathway(edgelist)
}
