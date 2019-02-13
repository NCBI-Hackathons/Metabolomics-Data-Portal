\name{p.table}
\alias{p.table}
\title{
Table of p-values of pathways
}
\description{
Table of p-values of pathways
}
\usage{
p.table(x, adj.method = NA, cutoff = ifelse(adj.method == "none", 0.01, 0.05))
}
\arguments{

  \item{x}{a \code{\link{cepa.all}} object}
  \item{adj.method}{methods in \code{\link[stats]{p.adjust}}, available methods are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}
  \item{cutoff}{cutoff for significance}

}
\details{
Since the p-values for each pathway are calculated for several centralities, the
whole p-values are represented as a table.

Also it can extract significant pathways only.
}
\value{
A data matrix where rows are pathways and columns are centralities.
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
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
p.table(res.ora)
p.table(res.ora, adj.method = "BH")

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
p.table(res.gsa)
}
}
