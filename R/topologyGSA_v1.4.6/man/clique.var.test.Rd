\name{clique.var.test}
\alias{clique.var.test}
\title{Testing the concentration matrices of all cliques of a pathway.}

\description{
Test the equality of concentration matrices in two experimental conditions for all cliques of a pathway
}
\usage{
clique.var.test(y1,y2,dag,alpha)
}
\arguments{
\item{y1}{a matrix  with n1 individuals (rows) in the first experimental condition and p genes (columns).}
\item{y2}{a matrix  with n2 individuals (rows) in the second experimental condition and p genes (columns).
The genes in the two experimental conditions must be the same.}
\item{dag}{graphNEL object, a directed acyclic graph (DAG) corresponding to the pathway of interest. See package \code{gRbase} for more details.}
\item{alpha}{significance level of the test.}
}

\details{
The function tests  the equality of the concentration matrices of each
clique of a pathway in two experimental conditions. The graph of a
pathway is first converted into a DAG, then the moral graph is obtained
and if the latter graph is decomposable then the test is performed on
all its cliques. If the moral graph is not decomposable, its
triangulated version is obtained and then the test is performed on all
its cliques.

The expression data may contain some genes differing
from those in the pathway: in such case the function automatically takes
the intersection between the two gene sets.

This function requires \code{gRBase} and \code{qpgraph} packages.
}
\value{
\item{p.value}{the observed value of the significance level.}
\item{cliques}{the list of cliques where the test is performed.}
\item{var.equal}{logical, a vector with a TRUE for a significant clique and
  a FALSE otherwise.}
\item{lambda.var}{the observed value of the statistic.}
\item{graph}{the triangulated and moral graphs.}
}
\author{M. Sofia Massa, Gabriele Sales}

\references{Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set analysis exploiting the topology of a pathway. BMC Systems Biology, 4:121
\url{http://www.biomedcentral.com/1752-0509/4/121}
}
\seealso{\code{\link{pathway.var.test}}, \code{\link{pathway.mean.test}},
\code{\link{clique.mean.test}},
}
\examples{
data(examples)

clique.var.test(y1, y2, dag_bcell, 0.05)
}
