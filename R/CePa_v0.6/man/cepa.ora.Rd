\name{cepa.ora}
\alias{cepa.ora}
\title{
Apply centrality-extended ORA on a single pathway
}
\description{
Apply centrality-extended ORA on a single pathway
}
\usage{
cepa.ora(dif, pc, bk = NULL, pathway = NULL, id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen))
    else if(mode(cen) == "name") deparse(cen)
    else cen,
    iter = 1000)
}
\arguments{

  \item{dif}{differential gene list}
  \item{pc}{a \code{pathway.catalogue} class object}
  \item{bk}{background gene list. If background gene list are not specified, use whole human genes}
  \item{pathway}{\code{\link[igraph]{igraphtest}} object or edge list}
  \item{id}{identify which pathway in the catalogue}
  \item{cen}{centrality measuments, it can ce a string, function, or function that has been quoted}
  \item{cen.name}{centrality measurement names. This argument should be set if the \code{cen} is a function.}
  \item{iter}{number of simulations}

}
\details{
The function is always called by \code{\link{cepa.ora.all}}. But you can still
use it if you realy want to analysis just one pathway under one centrality.
}
\value{
A \code{cepa} class object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{
data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)
}
}
