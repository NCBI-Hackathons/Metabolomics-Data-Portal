\name{plot.pathway.catalogue}
\alias{plot.pathway.catalogue}
\title{
plot pathway.catalogue object
}
\description{
plot pathway.catalogue object
}
\usage{
\method{plot}{pathway.catalogue}(x, ...)
}
\arguments{

  \item{x}{a \code{pathway.catalogue} object}
  \item{...}{other arguments}

}
\details{
There are three fugures: A) Distribution of the number of member genes in each node; 
B) Distribution of the number of nodes in which a single gene resides; C) Relationship 
between node count and gene count in biological pathways.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{set.pathway.catalogue}}
}
\examples{
data(PID.db)
NCI = PID.db$NCI
plot(NCI)
}
