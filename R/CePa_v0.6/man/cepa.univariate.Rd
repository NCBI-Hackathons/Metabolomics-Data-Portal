\name{cepa.univariate}
\alias{cepa.univariate}
\title{
Apply centrality-extended GSA on a single pathway
}
\description{
Apply centrality-extended GSA on a single pathway
}
\usage{
cepa.univariate(mat, label, pc, pathway = NULL, id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
    iter = 1000, nlevel = "tvalue_abs", plevel = "mean",
    node.level.from.expr = NULL, node.level.t.value = NULL,
    r.node.level.from.expr = NULL)
}
\arguments{

  \item{mat}{expression matrix in which rows are genes and columns are samples}
  \item{label}{a \code{\link{sampleLabel}} object identify the design of the microarray experiment}
  \item{pc}{a \code{pathway.catalogue} object storing information of pathways}
  \item{pathway}{\code{\link[igraph]{igraph}} object or edge list}
  \item{id}{identify the number of the pathway in the catalogue}
  \item{cen}{centrality measuments, it can ce a string, or function has been quote}
  \item{cen.name}{centrality measurement names}
  \item{nlevel}{node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs". Also self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.}
  \item{plevel}{pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank". Also, self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.}
  \item{node.level.from.expr}{for simplicity of computing}
  \item{node.level.t.value}{for simplicity of computing}
  \item{r.node.level.from.expr}{for simplicity of computing}
  \item{iter}{number of simulations}

}
\details{
The function is always called by \code{\link{cepa.univariate.all}}. But you can still
use it if you realy want to analysis just one pathway under one centrality.
}
\value{
A \code{\link{cepa}} class object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{

data(PID.db)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.univariate(mat = eset, label = label, pc = PID.db$NCI, id = 2)
}
}
