\name{PID.db}
\docType{data}
\alias{PID.db}
\title{
pathway catalogues from Pathway Interaction Database(PID)
}
\description{
pathway catalogues from Pathway Interaction Database(PID)
}
\usage{
data(PID.db)
}
\details{
The pathway data is parsed from XML format file provided by PID FTP site.
The perl code for parsing pathway data can be found at author's website
(\url{http://mcube.nju.edu.cn/jwang/lab/soft/cepa/index.html} ).

There are four pathway catalogues which are NCI_Nature, BioCarta, KEGG
and Reactome.

Each catalogue contains at least three members: pathway list (pathList), 
interaction list (interactionList) and mappings from node id to gene id (mapping). 
The pathway list contains a list of pathways in which each pathway is represented
by a list of interaction id. The interactions can be queried from the interaction list
by the interaction id. The interaction list data is represented as a data frame
in which the first column is the interaction id, the second column is the input
node id and the third column is the output node id. In real biological pathways,
a node in the pathway can be proteins, complex, families and none-gene nodes, so the 
mapping from node ids to gene ids is also provided. It must be noted that in this
package, gene symbol is selected as the primary gene id, so if users apply the PID.db
data, they should pay attension to the gene ids they transform.

Besides, in each catalogue, there also a node.name, node.type and version data.
The node.name provides the custumed name for each node. The node.type provides the 
type for each node (e.g. it is a complex or compound). The version provides the version
of the catalogue data.

Data has been updated to the lastest version by the day the package released (at 2012_07_19 09:34::20).

Note only part of pathways in the XML file are listed on the PID website. Also, 
we have set the minimum and maximum connected nodes when extracting pathways 
from PID, so not all the pathways listed on the PID website are in PID.db.
}
\value{
A list containing four component:

\describe{
  \item{NCI}{NCI_Nature-curated pathway catalogue}
  \item{BioCarta}{BioCarta pathway catalogue}
  \item{KEGG}{KEGG pathway catalogue}
  \item{Reactome}{Reactome pathway catalogue}
}

Each pathway catalogue is a \code{pathway.catalogue} class object. Each pathway
catalogue can be used directly in \code{\link{cepa.all}} and \code{\link{cepa}}
}
\section{Source}{
\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/NCI-Nature_Curated.xml.gz}

\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/BioCarta.xml.gz}

\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/Reactome.xml.gz}

\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/KEGG.xml.gz}}
\examples{
data(PID.db)
names(PID.db)
PID.db$NCI
plot(PID.db$NCI)
}
