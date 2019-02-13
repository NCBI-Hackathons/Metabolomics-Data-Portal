\name{pathway.mean.test}
\alias{pathway.mean.test}
\title{Testing the means of a pathway.}

\description{Test the equality of means in two experimental conditions for a pathway,
conditioning on the result of the test on the concentration matrices
}
\usage{pathway.mean.test(y1,y2,dag,alpha,perm.num=10000,variance=TRUE,paired=FALSE)}

\arguments{
\item{y1}{a matrix with n1 individuals (rows) in the first experimental condition and p genes (columns).}
\item{y2}{a matrix with n2 individuals (rows) in the second experimental condition and p genes (columns).
The genes in the two experimental conditions must be the same.}
\item{dag}{graphNEL object, a directed acyclic graph (DAG) corresponding to the pathway of interest. See package \code{gRbase} for more details.}
\item{alpha}{significance level of the test.}
\item{variance}{logical flag. If \code{TRUE} estimated variances s1  and s2 are reported.}
\item{perm.num}{number of permutations of the test.}
\item{paired}{logical flag. If \code{TRUE} Hotelling test for paired samples is calculated and the test on the variances is not performed.}
}
\details{The function tests  the equality of the means of  two
  experimental conditions, conditioning on the result of
  \code{\link{pathway.var.test}}.

  The expression data may contain some genes differing
  from those in the pathway: in such case the function automatically takes
  the intersection between the two gene sets.

  This function requires \code{gRBase} and \code{qpgraph} packages.
}

\value{
  \item{p.value}{the p-value of the test on the mean.}
  \item{p.value.var}{the p-value of the test on the variance. If paired is \code{TRUE}, this is not returned.}
  \item{lambda.value}{the observed value of the test statistic for the
    variance. If paired is \code{TRUE}, this is not returned.}
  \item{qchisq.value}{the 1 - alpha quantile value of the null
    distribution of the test statistic on the variance. If paired is \code{TRUE}, this is not returned.}
  \item{cli.moral}{a list containing the cliques of the moral graph.}
  \item{var.equal}{logical flag. If \code{TRUE} variances are heteroschedastic. If paired is \code{TRUE}, this is not returned.}
  \item{graph}{a graphNEL object representing the moral graph.}
  \item{df.var}{the degrees of freedom of the null distribution for the
    variance test. If paired is \code{TRUE}, this is not returned.}
  \item{df.mean}{the degrees of freedom of the null distribution for the
    mean test. If paired is \code{TRUE}, this is not returned.}
  \item{t.value}{the observed value of the test statistic for the mean.}
}


\author{M. Sofia Massa,  Gabriele Sales}
\references{Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set analysis exploiting the topology of a pathway. BMC Systems Biology, 4:121
\url{http://www.biomedcentral.com/1752-0509/4/121}
}
\seealso{\code{\link{pathway.var.test}}, \code{\link{clique.var.test}},
\code{\link{clique.mean.test}},
}
\examples{
data(examples)

pathway.mean.test(y1, y2, dag_bcell, 0.05, 100)
}
