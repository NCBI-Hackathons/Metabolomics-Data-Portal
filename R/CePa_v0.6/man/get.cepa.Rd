\name{get.cepa}
\alias{get.cepa}
\title{
get single cepa object from cepa.all object
}
\description{
get single cepa object from cepa.all object
}
\usage{
get.cepa(x, id = NULL, cen = 1)
}
\arguments{

  \item{x}{a \code{\link{cepa.all}} object}
  \item{id}{index or the name of the pathway}
  \item{cen}{index or the name of the centrality}

}
\details{
The `cepa.all object contains the result for pathways under several centrality
measurements. In \code{\link{cepa.all}} object, each pathway under a specific centrality
is a single \code{\link{cepa}} object. The \code{\link{get.cepa}} function is used to get the \code{\link{cepa}}
object from the \code{\link{cepa.all}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{cepa}}, \code{\link{cepa.all}}
}
\examples{
\dontrun{
data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
ora = get.cepa(res.ora, id = 5, cen = 3)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
gsa = get.cepa(res.gsa, id = 5, cen = 3)
}
}
