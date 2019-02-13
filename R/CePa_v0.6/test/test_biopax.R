
for(script in dir("R")) {
	source(paste("R/", script, sep = ""))
}

library(rBiopaxParser)
library(igraph)
library(Rgraphviz)

nci = import_biopax("test/NCI-Nature_Curated.bp2.owl")
biocarta = import_biopax("test/BioCarta.bp2.owl")
reactome = import_biopax("test/Reactome.bp2.owl")
kegg = import_biopax("test/KEGG.bp2.owl")
