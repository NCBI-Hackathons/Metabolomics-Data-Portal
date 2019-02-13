\name{cepa.all.parallel}
\alias{cepa.all.parallel}
\title{
use CePa package through parallel computing
}
\description{
use CePa package through parallel computing
}
\usage{
cepa.all.parallel(dif = NULL, bk = NULL, mat = NULL, label = NULL,
    pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)),
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2)
}
\arguments{

  \item{dif}{differential gene list}
  \item{bk}{background gene list. If background gene list are not specified, use whole human genes}
  \item{mat}{expression matrix in which rows are genes and columns are samples}
  \item{label}{a \code{\link{sampleLabel}} object identify the design of the microarray experiment}
  \item{pc}{a \code{pathway.catalogue} object storing information of pathways}
  \item{cen}{centrality measuments, it can ce a string, or a function}
  \item{cen.name}{centrality measurement names. By default it is parsed from \code{cen} argument}
  \item{nlevel}{node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs". Also self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.}
  \item{plevel}{pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank". Also, self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.}
  \item{iter}{number of simulations}
  \item{ncores}{number of cores for parallel computing}

}
\details{
The function divides the pathway list into several parts and each part is sent to a core for 
parallel computing.

The package for parallel computing is \code{snow}.

Note: there may be warnings saying connections not closed. In fact I have closed
connections after the parallel computing is done. I don't know why this
happens. Maybe you breaked the computing ahead manually. However it does not matter 
unless you have obsessive compulsive disorder.
}
\value{
A \code{\link{cepa.all}} class object
}
\references{
Gu Z, Liu J, Cao K, Zhang J, Wang J. Centrality-based pathway enrichment: a systematic 
approach for finding significant pathways dominated by key genes. BMC Syst Biol. 2012 Jun 6;6(1):56.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
cepa.all
}
\examples{
\dontrun{
data(PID.db)
# ORA extension
data(gene.list)
res.ora = cepa.all.parallel(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, ncores = 4)
# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53_symbol.gct")
label = read.cls("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", 
    treatment="MUT", control="WT")
res.gsa = cepa.all.parallel(mat = eset, label = label, pc = PID.db$NCI, ncores = 4)
}
}
