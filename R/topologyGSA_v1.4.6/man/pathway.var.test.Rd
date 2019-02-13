\name{pathway.var.test}
\alias{pathway.var.test}
\title{Testing the concentration matrices of a pathway.}

\description{Test the equality of concentration matrices in two experimental conditions for a pathway
}
\usage{
pathway.var.test(y1,y2,dag,alpha,variance=FALSE,s1=NULL,s2=NULL)
}

\arguments{
\item{y1}{a matrix with n1 individuals (rows) in the first experimental condition and p genes (columns).}
\item{y2}{a matrix with n2 individuals (rows) in the second experimental condition and p genes (columns).
The genes in the two experimental conditions must be the same.}
\item{dag}{graphNEL object, directed acyclic graph (DAG) corresponding to the pathway of interest. See package \code{gRbase} for more details.}
\item{alpha}{significance level of the test.}
\item{variance}{logical flag. If \code{TRUE} the estimates of the
  covariance matrices are included in the result.}
\item{s1}{y1 covariance matrix estimation.}
\item{s2}{y2 covariance matrix estimation.}
}
\details{The graph of a pathway is first converted into a DAG and then
  into a moral graph. The data is modelled with two Gaussian graphical
  models with zero mean and graph provided by the moral graph. The
  function tests the equality of the two concentration matrices
  (inverse of the covariance matrices).

  The expression data may contain some genes differing
  from those in the pathway: in such case the function automatically takes
  the intersection between the two gene sets.

  A necessary condition for the existence of the covariance estimates is that
  the number of statistical units (samples) is greater than the number of
  variables. If this is not the case, penalized techniques for estimating
  \eqn{\hat{\Sigma}_{1}^{-1}}{s1^(-1)} and \eqn{\hat{\Sigma}_{2}^{-1}}{s2^(-1)}
  have to be employed, that are currently not provided by the package. In theses
  cases, one can perform penalized estimation of
  \eqn{\hat{\Sigma}_{1}^{-1}}{s1^(-1)} and \eqn{\hat{\Sigma}_{2}^{-1}}{s2^(-1)}
  outside \pkg{topologyGSA}, and then provide such estimates as input arguments
  to the function \code{pathway.var.test} to compute the value of the test for
  homogeneity. In this case, computation of the p-value deserves attention,
  as standard results on the asymptotic distribution of the test statistic may
  no longer be valid. Therefore, computation of the p-value has to be dealt
  with by the user.

  This function requires \code{gRBase} and \code{qpgraph} packages.
}
\value{
\item{lambda.value}{the observed value of the test statistic.}
\item{qchisq.value}{the 1 - alpha quantile value of the null
  distribution of the test statistic on the variance.}
\item{p.value}{the observed value of the significance level of the test.}
\item{cli.moral}{a list containing the cliques of the moral graph.}
\item{var.equal}{logical flag. If TRUE variances are heteroschedastic.}
\item{graph}{a graphNEL object representing the moral graph.}
\item{df}{the degrees of freedom of the null distribution.}
\item{s1}{if variance=TRUE, the estimate of y1 covariance.}
\item{s2}{if variance=TRUE, the estimate of y2 covariance.}
}

\author{M. Sofia Massa, Gabriele Sales}

\references{Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set
  analysis exploiting the topology of a pathway. BMC Systems Biology,
  4:121
  \url{http://www.biomedcentral.com/1752-0509/4/121}

Lauritzen, S.L. (1996). Graphical models. Clarendon Press, Oxford.
}
\seealso{\code{\link{pathway.mean.test}}, \code{\link{clique.var.test}},
\code{\link{clique.mean.test}}.
}
\examples{
data(examples)

pathway.var.test(y1, y2, dag_bcell, 0.05)
}
