\name{sampleLabel}
\alias{sampleLabel}
\title{
Generate data structure of sample labels
}
\description{
Generate data structure of sample labels
}
\usage{
sampleLabel(label, treatment, control)
}
\arguments{

  \item{label}{sample label vector}
  \item{treatment}{treatment label}
  \item{control}{control label}

}
\details{
Since sample label will not be modified in the analysis, this function is used to
integrate all the label information in one single data object.
}
\value{
A \code{sampleLabel} class object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
sampleLabel(c("A", "B", "B", "A", "A", "A", "B", "B"), treatment = "A", control = "B")
}
