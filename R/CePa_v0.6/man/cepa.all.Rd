\name{cepa.all}
\alias{cepa.all}
\title{
Apply CePa algorithm on a list of pathways under multiple centralities
}
\description{
Apply CePa algorithm on a list of pathways under multiple centralities
}
\usage{
cepa.all(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)),
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000)
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

}
\section{Item}{
cen.name
}
\section{Code}{
cen
}
\section{Nlevel}{
node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs". Also self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.
}
\section{Plevel}{
pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank". Also, self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.
}
\section{Iter}{
number of simulations
}
\section{Item}{
cen.name
}
\section{Code}{
cen
}
\section{Nlevel}{
node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs". Also self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.
}
\section{Plevel}{
pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank". Also, self-defined functions are allowed, see \code{\link{cepa.univariate.all}} for detail.
}
\section{Iter}{
number of simulations
}
\details{
All the calculation can be achieved by this function. The function is wrapper of 
both ORA extension and GSA extension. It chooses corresponding procedure according 
to the arguments specified. If the arguments contain gene lists, then the calculation 
is sent to functions doing ORA extension. While if the arguments contain an expression 
matrix and a phenotype label, the GSA extension is evoked.

The function is a wrapper of \code{\link{cepa.ora.all}} and \code{\link{cepa.univariate.all}}.

This is the core function of the package. User can refer to the vignette to find
how to use it (\code{vignette("CePa")}).

If \code{dif}, \code{bk}, \code{pc}, \code{cen}, \code{cen.name} and \code{iter}
are specified, the arguments are passed to \code{cepa.ora.all}. The centrality-extension 
of over-representation analysis (ORA) will be applied on the list of differential genes.

If \code{mat}, \code{label}, \code{pc}, \code{cen}, \code{cen.name}, \code{nlevel},
\code{plevel} and \code{iter} are specified, the arguments are passed to \code{cepa.univariate.all}.
The centrality-extension of gene-set analysis (GSA) will be applied on the whole gene expressions.

There is a parallel version of the function: \code{\link{cepa.all.parallel}}.
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
\code{\link{cepa}}, \code{\link{cepa.ora.all}}, \code{\link{cepa.univariate.all}}, \code{\link{cepa.all.parallel}}
}
\examples{
\dontrun{

data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53_symbol.gct")
label = read.cls("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
}
}
