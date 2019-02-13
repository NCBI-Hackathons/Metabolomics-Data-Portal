\name{report}
\alias{report}
\title{
Generate report for CePa analysis
}
\description{
Generate report for CePa analysis
}
\usage{
report(x, adj.method = "none", cutoff = ifelse(adj.method == "none", 0.01, 0.05),
    template.file = system.file(package = "CePa", "extdata", "cepa.template"),
    only.sig = TRUE, dir.path = NULL, ...)
}
\arguments{

  \item{x}{a \code{\link{cepa.all}} object}
  \item{adj.method}{methods in \code{\link[stats]{p.adjust}}, available methods are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}
  \item{cutoff}{cutoff for significance}
  \item{template.file}{path of the template file}
  \item{only.sig}{whether to generate detailed report for every pathway. If it is set to FALSE, the page for every pathway under every centrality would be generated (there would be so many images!).}
  \item{dir.path}{dir name}
  \item{...}{other arguments}

}
\details{
The report is in HTML format that you can view it in you web browser. Networks
for pathways can be visualized interactively (by using Cytoscape Web, in which 
you can drag the network, zoom in and zoom out the network). To load Flash Player
successful in you browser, you may need to set the Flash security settings on your
machine. See \url{http://cytoscapeweb.cytoscape.org/tutorial} and change the settings
via \url{http://www.macromedia.com/support/documentation/en/flashplayer/help/settings_manager04.html} .

The report would locate at the current working directory. View the report
by clicking \code{index.html} in the report directory.

There is also another popular method qvalue to adjust p-values. Turn to \code{\link{plot.cepa.all}}
to find out how to use qvalue.
}
\section{Source}{
\url{http://cytoscapeweb.cytoscape.org/}}
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
report(res.ora)
}
}
