\name{print.pathway.catalogue}
\alias{print.pathway.catalogue}
\title{
print pathway.catalogue object
}
\description{
print pathway.catalogue object
}
\usage{
\method{print}{pathway.catalogue}(x, ...)
}
\arguments{

  \item{x}{a \code{pathway.catalogue} object}
  \item{...}{other arguments}

}
\details{
Simply print the number of pathways in the catalogue
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
NCI
}
