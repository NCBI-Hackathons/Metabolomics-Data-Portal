\name{plotNull}
\alias{plotNull}
\title{
Plot the null distribution of the pathway score
}
\description{
Plot the null distribution of the pathway score
}
\usage{
plotNull(x)
}
\arguments{

  \item{x}{a \code{\link{cepa}} object}

}
\details{
There are two figures in the plotting.

\itemize{
  \item A) Distribution of node score in the pathway under simulation. Since a pathway contains a list of nodes. The distribution of node score for the pathway in each simulation is measures by maximum value,  the 75th quantile, median value and minimum value. The distribution of node score for  the pathway in the real data is highlighted.
}

\itemize{
  \item B) Histogram of simulated pathway scores. 
}

The function is always called through \code{\link{plot.cepa.all}} and \code{\link{plot.cepa}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{cepa}}, \code{\link{plot.cepa}}
}
\examples{
\dontrun{
data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
ora = get.cepa(res.ora, id = 5, cen = 3)
plotNull(ora)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
gsa = get.cepa(res.gsa, id = 5, cen = 3)
plotNull(gsa)
}
}
