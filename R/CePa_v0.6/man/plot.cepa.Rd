\name{plot.cepa}
\alias{plot.cepa}
\title{
Plot the cepa object
}
\description{
Plot the cepa object
}
\usage{
\method{plot}{cepa}(x, type = c("graph", "null"), ...)
}
\arguments{

  \item{x}{a \code{\link{cepa}} object}
  \item{type}{identify the type for the plot}
  \item{...}{arguments passed to \code{\link{plotGraph}}}

}
\details{
The function is wrapper of \code{\link{plotGraph}} and \code{\link{plotNull}}.
If type is specified to "graph", the graph of the network will be plotted (see \code{\link{plotGraph}} for details).
If type is specified to "null", the null distribution of the pathway score
in the pathway will be plotted (see \code{\link{plotNull}} for details).
}
\value{
if type is set to "graph", the function will return a \code{\link[igraph]{igraph}} object or a \code{graphML} object of the pathway. Else it is NULL.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{cepa}}, \code{\link{plotNull}}, \code{\link{plotGraph}}
}
\examples{
\dontrun{

data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)
plot(res.ora)
plot(res.ora, type = "null")

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa(mat = eset, label = label, pc = PID.db$NCI, id = 2)
plot(res.gsa, type = "null")
}
}
