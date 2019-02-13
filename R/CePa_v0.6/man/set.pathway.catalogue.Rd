\name{set.pathway.catalogue}
\alias{set.pathway.catalogue}
\title{
store pathway data and pre-processing
}
\description{
store pathway data and pre-processing
}
\usage{
set.pathway.catalogue(pathList, interactionList, mapping,
    min.node = 5, max.node = 500, min.gene = min.node, max.gene = max.node, ...)
}
\arguments{

  \item{pathList}{list of pathways}
  \item{interactionList}{list of interactions}
  \item{mapping}{a data frame or matrix providing mappings from gene id to pathway node id.  The first column is node id and the second column is gene id.}
  \item{min.node}{minimum number of connected nodes in each pathway}
  \item{max.node}{maximum number of connected nodes in each pathway}
  \item{min.gene}{minimum number of genes in each pathway}
  \item{max.gene}{maximum number of genes in each pathway}
  \item{...}{other arguments, should have names, these data will be stored as a list member in the returned value from the function}

}
\details{
The pathway data will not be changed in the analysis, so the pathway data is integrated
in one single data object by this function. Also, the function will do a little 
preprocess of the pathway data.

Basicly, a pathway contains a list of interactions. The pathList argument
is a list where elements in the list is the vector of interaction IDs 
in the pathway. The interactions in the pathway can be got from a interaction
list pool represented as interactionList argument. The interactionList argument
stores the total interaction list in the pathway catalogue. It represents
as a three columns data frame or matrix where the first column is the interaction id,
the second column is the input node id and the third column is the output node id.

The mapping data frame provide the mapping from node id to gene id. The 
first column is the node id and the second column is the gene id.

Besides the pathList, interactionList and mapping arguments, more arguments can
be added to the function. These new data will be stored as the member of the list
that returned by the function. E.g., in the \code{\link{PID.db}} data, each catalogue
is a \code{pathway.catalogue} object. Except the pathList, interactionList and mapping arguments,
there are also node.name, node.type and version arguments.

The summary can be visualized by \code{\link{plot.pathway.catalogue}}.
}
\value{
A \code{pathway.catalogue} class object
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
catalogue = set.pathway.catalogue(pathList = PID.db$NCI$pathList[1:20],
            interactionList = PID.db$NCI$intertionList, mapping = PID.db$NCI$mapping)
}
}
