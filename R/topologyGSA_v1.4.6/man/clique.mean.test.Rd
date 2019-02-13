\name{clique.mean.test}
\alias{clique.mean.test}
\title{Comparison of all the cliques of a pathway.}
\description{
Test the equality of the means in two experimental conditions for all the cliques of a pathway, conditioning on the results of the test on the concentration matrices
}
\usage{
clique.mean.test(y1,y2,dag,alpha,perm.num=1000,paired=FALSE)
}

\arguments{

\item{y1}{a matrix  with n1 individuals (rows) in the first experimental condition and p genes (columns). }
\item{y2}{a matrix  with n2 individuals (rows) in the second experimental condition and p genes (columns).
The genes in the two experimental conditions must be the same.}
\item{dag}{ graphNEL object, a directed acyclic graph (DAG) corresponding to the pathway of interest. See package \code{gRbase} for more details.}
\item{alpha}{ significance level of the test.}
\item{perm.num}{number of permutation used by the test.}
\item{paired}{logical flag. If \code{TRUE} statistical tests for paired samples is calculated and the test on the variances is not performed.}
}
\details{
  The function tests  the equality of the means of each clique of a
  pathway in two experimental conditions, conditioning on the results of
  \code{clique.var.test}.
  The graph of a pathway is first converted into a DAG, then the moral
  graph is obtained and if the latter graph is decomposable then the
  test is performed on all its cliques. If the moral graph is not
  decomposable, its triangulated version is obtained and then the test
  is performed on all its cliques.

  The expression data may contain some genes differing
  from those in the pathway: in such case the function automatically takes
  the intersection between the two gene sets.

  This function requires \code{gRBase} and \code{qpgraph} packages.
}
\value{
  \item{p.value}{p-value of the test on the mean for each clique.}
  \item{p.value.var}{p-value of the test on the variance for each clique. If paired is \code{TRUE}, this is not returned.}
  \item{cliques}{the list of the corresponding cliques.}
  \item{var.equal}{logical, a vector with a TRUE for a significant clique and a FALSE otherwise. If paired is \code{TRUE}, this is not returned.}
  \item{t.value}{the observed values of the statistic to test the mean
    for each clique.}
  \item{lambda.value}{the observed values of the statistic to test the
    variance for each clique. If paired is \code{TRUE}, this is not returned.}
  \item{graph}{the triangularized and moralized graph.}
}

\author{M. Sofia Massa, Gabriele Sales}

\references{
Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set analysis exploiting the topology of a pathway. BMC Systems Biology, 4:121
\url{http://www.biomedcentral.com/1752-0509/4/121}
}
\seealso{
\code{\link{pathway.var.test}},\code{\link{pathway.mean.test}},
\code{\link{clique.var.test}},
}
\examples{
data(examples)

clique.mean.test(y1, y2, dag_bcell, 0.05)
}
