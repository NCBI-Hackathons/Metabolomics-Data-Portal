\name{read.gct}
\alias{read.gct}
\title{
Read GCT format file which stores the expression values
}
\description{
Read GCT format file which stores the expression values
}
\usage{
read.gct(file)
}
\arguments{

  \item{file}{gct file path}

}
\details{
The GCT format is a tab delimited file format that stores the expression value matrix.
The first line of the file is the version number which always be #1.2.
The second line is the number of the size of genes and samples, 
seperated by space, usually for the initiation of reading the expression matrix.
The third line contains a list of identifiers for the samples associated 
with each of the columns in the remainder of the file.
From the fourth line will be the expression value of each gene.

GCT file is used together with CLS file.
}
\value{
A matrix of the expression values, with rows correponding to genes and cols to
samples.
}
\section{Source}{
\url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
# expression data stored in a gct format file
# P53_symbol.gct can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53_symbol.gct")
head(eset)
}
}
