\name{print.cepa}
\alias{print.cepa}
\title{
\preformatted{
  print the cepa object  }
}
\description{
\preformatted{
  print the cepa object  }
}
\usage{
\method{print}{cepa}(x, ...)
}
\arguments{

  \item{x}{a \code{\link{cepa}} object}
  \item{...}{other arguments}

}
\details{
The function print procedure of the analysis, the centrality and the p-value for the pathway.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{cepa}}
}
\examples{
\dontrun{

data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)
res.ora

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa(mat = eset, label = label, pc = PID.db$NCI, id = 2)
res.gsa
}
}
