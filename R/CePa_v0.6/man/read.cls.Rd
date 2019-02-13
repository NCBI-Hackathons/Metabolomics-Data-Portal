\name{read.cls}
\alias{read.cls}
\title{
Read CLS file which stores the phenotype data
}
\description{
Read CLS file which stores the phenotype data
}
\usage{
read.cls(file, treatment, control)
}
\arguments{

  \item{file}{cls file path}
  \item{treatment}{string of treatment label in cls file}
  \item{control}{string of control label in cls file}

}
\details{
The CLS file format defines the phenotype data of microarray experiments.
The first line is the number of samples, number of classes and the third number always be 1.
These three numbers are seperated by spaces or tabs.
The second line begins with #. The next two strings usually are the label of the phenotype. 
The third line is the label of each samples where same label represents 
the same class.

The first and the second line is ignored by this function 
and class labels are taken from the factor of the vector parsed from the third line.
}
\value{
A \code{sampleLabel} class object
}
\section{Source}{
\url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{read.gct}}, \code{\link{sampleLabel}}
}
\examples{
\dontrun{
# P53.cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
label = read.cls("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", 
    treatment="MUT", control="WT")
}
}
