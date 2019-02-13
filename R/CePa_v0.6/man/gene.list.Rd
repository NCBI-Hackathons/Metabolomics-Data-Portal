\name{gene.list}
\docType{data}
\alias{gene.list}
\title{
Differential gene list and background gene list
}
\description{
Differential gene list and background gene list
}
\usage{
data(gene.list)
}
\details{
Differential gene list and background gene list was extracted from
microarray data from GEO database. The accession number for the data set
is GSE22058. The t-test was applied to find differentially expressed genes.
Top 2000 genes were selected as the gene list.
}
\value{
A list containing two componets:

\describe{
  \item{\code{bk}}{background gene list, gene symbol}
  \item{\code{dif}}{differentially expressed gene list, gene symbol}
}
}
\section{Source}{
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22058}}
\examples{
data(gene.list)
names(gene.list)
}
