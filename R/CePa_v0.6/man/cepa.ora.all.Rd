\name{cepa.ora.all}
\alias{cepa.ora.all}
\title{
Apply centrality-extented ORA on a list of pathways
}
\description{
Apply centrality-extented ORA on a list of pathways
}
\usage{
cepa.ora.all(dif, pc, bk = NULL, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)),
    iter = 1000)
}
\arguments{

  \item{dif}{differential gene list}
  \item{pc}{a \code{pathway.catalogue} class object}
  \item{bk}{background gene list. If background gene list are not specified, use whole human genes}
  \item{cen}{centrality measuments, it can ce a string, or a function}
  \item{cen.name}{centrality measurement names. By default it is parsed from \code{cen} argument}
  \item{iter}{number of simulations}

}
\details{
The traditional over-representation analysis (ORA) to find significant pathways 
uses a 2x2 contingency table to test the independency of genes belonging to a 
functional category and these genes being differentially expressed, usually by 
Fisher's exact test. The ORA only consider the number of genes and the function
extend traditional ORA with network centralities.

The differential gene list and the background gene list should be indicated
with the same identifiers (e.g. gene symbol or refseq ID). All genes in
the differential gene list should exist in the background gene list. If users 
use the \code{\link{PID.db}} data, all genes should be formatted in gene symbol.

If the centrality measurement is set as a string, only pre-defined "equal.weight",
"in.degree", "out.degree", "degree", "betweenness", "in.reach", "out.reach",
"reach", "in.spread", "out.spread" and "spread" are allowed. More centrality
measurements can be used by setting it as a function (such as closeness,
cluster coefficient). In the function, we recommand users choose
at least two centrality measurements. The default centralities are "equal.weight",
"in.degree", "out.degree", "betweenness", "in.reach" and "out.reach".

However, in most circumstance, the function is called by \code{\link{cepa.all}}.
}
\value{
A \code{\link{cepa.all}} class object
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
res.ora = cepa.ora.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
}
}
