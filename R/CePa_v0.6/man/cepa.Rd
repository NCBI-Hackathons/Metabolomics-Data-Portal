\name{cepa}
\alias{cepa}
\title{
Apply CePa algorithm on a single pathway
}
\description{
Apply CePa algorithm on a single pathway
}
\usage{
cepa(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, pathway = NULL,
    id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen))
    else if(mode(cen) == "name") deparse(cen) else cen,
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000)
}
\arguments{

  \item{dif}{differential gene list}
  \item{bk}{background gene list. If background gene list are not specified, use whole human genes}
  \item{mat}{expression matrix in which rows are genes and columns are samples}
  \item{label}{a \code{\link{sampleLabel}} object identify the design of the microarray experiment}
  \item{pc}{a \code{pathway.catalogue} object storing information of pathways}
  \item{pathway}{an \code{\link[igraph]{igraphtest}} object or edge list}
  \item{id}{identify which pathway should be analysis in the pathway catalogue}
  \item{cen}{centrality measuments, it can ce a string, or function has been quote}
  \item{cen.name}{centrality measurement names. This argument should be set if the \code{cen} is a function.}
  \item{nlevel}{node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs". Also self-defined functions are allowed, see \code{\link{cepa.univariate}} for detail.}
  \item{plevel}{pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank". Also, self-defined functions are allowed, see \code{\link{cepa.univariate}} for detail.}
  \item{iter}{number of simulations}

}
\details{
The function is a wrapper of \code{\link{cepa.ora}} and \code{\link{cepa.univariate}}.
Selection of which function depends on the arguments specified.

If \code{dif}, \code{bk}, \code{pc}, \code{pathway}, \code{id}, \code{cen}, \code{cen.name} and \code{iter}
are specified, the arguments are passed to \code{\link{cepa.ora}}. The centrality-extension 
of over-representation analysis (ORA) will be applied on the list of differential genes.

If \code{mat}, \code{label}, \code{pc}, \code{pathway}, \code{id}, \code{cen}, \code{cen.name}, \code{nlevel},
\code{plevel} and \code{iter} are specified, the arguments are passed to \code{\link{cepa.univariate}}.
The centrality-extension of gene-set analysis (GSA) will be applied on the whole gene expressions.

This function is always called by \code{\link{cepa.all}}. But you can still use it
if you want to analysis a single pathway under a specific centrality.
}
\value{
A \code{\link{cepa}} class object
}
\seealso{
\code{\link{cepa.all}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{

data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will take about 45 min
res.gsa = cepa(mat = eset, label = label, pc = PID.db$NCI, id = 2)
}
}
