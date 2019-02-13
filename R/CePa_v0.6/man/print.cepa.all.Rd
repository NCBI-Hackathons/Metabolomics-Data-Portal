\name{print.cepa.all}
\alias{print.cepa.all}
\title{
print the cepa.all object
}
\description{
print the cepa.all object
}
\usage{
\method{print}{cepa.all}(x, ...)
}
\arguments{

  \item{x}{a \code{\link{cepa.all}} object}
  \item{...}{other arguments}

}
\details{
The function print the number of significant pathways under various
centrality measures at p-value <= 0.01.
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
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
res.ora

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
res.gsa
}
}
