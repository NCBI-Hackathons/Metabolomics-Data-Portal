# == title
# names of the pathway nodes
#
# == param
# -pathway an `igraph::igraph` object
#
# == details
# If nodes in the pathway have names, then it returns a vector of nodes
# names. If nodes in the pathway have no name, it just returns the index
# of nodes (start from 1, after igraph version 0.6).
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# interaction = rbind(c("a", "b"),
#                     c("a", "c"))
# g = generate.pathway(interaction)
# pathway.nodes(g)

# pathway = barabasi.game(200)
# pathway.nodes(pathway)
pathway.nodes = function(pathway) {
    if(class(pathway) != "igraph") {
        stop("Wrong class for pathway.\n")
    }
    
    n = vcount(pathway)
    name = get.vertex.attribute(pathway, "name")
    
    if(length(name)) {
        return(name)
    }
    else {
        return(1:n)
    }
}

# == title
# Generate igraph object from edge list
#
# == param
# -el edge list, matrix with two columns. The first column is the input node
#      and the second column is the output node.
#
# == details
# The function is a wrapper of `igraph::graph.edgelist` and it generates
# a directed graph.
#    
# In the function, repeated edged for two nodes will be eliminated.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `cepa`, `igraph::graph.edgelist`
#
# == example
# edgelist = rbind(c("a", "b"), c("a", "b"), c("a", "c"))
# g = generate.pathway(edgelist)
generate.pathway = function(el) {
    
    if(dim(el)[2] != 2) {
        stop("Second dimension of edgelist should be 2.\n")
    }
    
    # remove duplicat connections between two nodes
    el.vector = apply(el, 1, paste, collapse="--")
    el.vector = unique(el.vector)
    
    el = t(as.matrix(as.data.frame(strsplit(el.vector, "--"))))
    
    g = graph.edgelist(el, directed = TRUE)
    
    return(g)
}


# == title
# store pathway data and pre-processing
#
# == param
# -pathList list of pathways
# -interactionList list of interactions
# -mapping a data frame or matrix providing mappings from gene id to pathway node id. 
#      The first column is node id and the second column is gene id.
# -min.node minimum number of connected nodes in each pathway
# -max.node maximum number of connected nodes in each pathway
# -min.gene minimum number of genes in each pathway
# -max.gene maximum number of genes in each pathway
# -... other arguments, should have names, these data will be stored as
#     a list member in the returned value from the function
#
# == details
# The pathway data will not be changed in the analysis, so the pathway data is integrated
# in one single data object by this function. Also, the function will do a little 
# preprocess of the pathway data.
#
# Basicly, a pathway contains a list of interactions. The pathList argument
# is a list where elements in the list is the vector of interaction IDs 
# in the pathway. The interactions in the pathway can be got from a interaction
# list pool represented as interactionList argument. The interactionList argument
# stores the total interaction list in the pathway catalogue. It represents
# as a three columns data frame or matrix where the first column is the interaction id,
# the second column is the input node id and the third column is the output node id.
#
# The mapping data frame provide the mapping from node id to gene id. The 
# first column is the node id and the second column is the gene id.
#
# Besides the pathList, interactionList and mapping arguments, more arguments can
# be added to the function. These new data will be stored as the member of the list
# that returned by the function. E.g., in the `PID.db` data, each catalogue
# is a ``pathway.catalogue`` object. Except the pathList, interactionList and mapping arguments,
# there are also node.name, node.type and version arguments.
#
# The summary can be visualized by `plot.pathway.catalogue`.
#
# == value
# A ``pathway.catalogue`` class object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `cepa.all`
#
# == example
# \dontrun{
# data(PID.db)
# catalogue = set.pathway.catalogue(pathList = PID.db$NCI$pathList[1:20],
#             interactionList = PID.db$NCI$intertionList, mapping = PID.db$NCI$mapping)
# }
set.pathway.catalogue = function(pathList, interactionList, mapping,
    min.node = 5, max.node = 500, min.gene = min.node, max.gene = max.node, ...) {
    
    # pathway should be list
    if(!is.list(pathList)) {
        stop("pathList should be a list.\n")
    }
    
    # pathway should have name
    if(is.null(names(pathList))) {
        stop("pathList should have names.\n")
    }
    
    # interactionList is a matrix
    if(length(dim(interactionList)) != 2) {
        stop("interactionList should be two dimension matrix.\n")
    }
    
    # interactionList contains three columns
    if(dim(interactionList)[2] != 3) {
        stop("interactinList should contain 3 columns.\n")
    }
    
    l = sapply(pathList, function(x) {
                   it = interactionList[interactionList[, 1] %in% x, 2:3]   # interaction list for the pathway
                   node = unique(c(it[, 1], it[, 2]))                       # nodes in the pathway
                   l.node = length(node)                                    # number of nodes
                   gene = unique(mapping[mapping[,1] %in% node, 2])         # genes in the pathway
                   l.gene = length(gene)                                    # number of genes
                   return(c(l.node, l.gene))
               })
    pathList = pathList[l[1, ] >= min.node & l[1, ] <= max.node
                        & l[2, ] >= min.gene & l[2, ] <= max.gene]
    
    if(length(pathList) == 0) {
        warning("Your pathway catalogue is empty!")
    }
    
    pc = list(pathList = pathList, 
              interactionList = interactionList,
              mapping = mapping,
              ...)        # name of the catalogue, such as KEGG, PID ...
    class(pc) = "pathway.catalogue"
    return(pc)
}

# == title
# print pathway.catalogue object
#
# == param
# -x a ``pathway.catalogue`` object
# -... other arguments
#
# == details
# Simply print the number of pathways in the catalogue
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `set.pathway.catalogue`
#
# == example
# data(PID.db)
# NCI = PID.db$NCI
# NCI
print.pathway.catalogue = function(x, ...) {
    cat("\n  The catalogue contains", length(x$pathList), "pathways.\n\n")
}

# == title
# plot pathway.catalogue object
#
# == param
# -x a ``pathway.catalogue`` object
# -... other arguments
#
# == details
# There are three fugures: A) Distribution of the number of member genes in each node; 
# B) Distribution of the number of nodes in which a single gene resides; C) Relationship 
# between node count and gene count in biological pathways. 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `set.pathway.catalogue`
#
# == example
# data(PID.db)
# NCI = PID.db$NCI
# plot(NCI)
plot.pathway.catalogue = function(x, ...) {

    r1 = numeric(0)
    r2 = numeric(0)
    i = 0
    
    # first column is the number of genes in each pathway
    # second column is the number of nodes in each pathway
    gn = matrix(0, nrow=length(x$pathList), ncol=2)
    colnames(gn) = c("node", "gene")
    for(pa in x$pathList) {
        # interaction ID list
        int = pa
        # node IDs in the pathway
        l  = x$interactionList[, 1] %in% int
        node = unique(c(x$interactionList[l, 2], x$interactionList[l, 3]))
        # mapping in this pathway
        m = x$mapping[x$mapping[, 1] %in% node, ]
        # number of nodes contain k genes
        ta = table(table(m[, 1]))
        # how many genes that a node contains
        tname = names(ta)
        tname = as.integer(tname)
        r1 = c(r1, rep(tname, ta))
        
        ta = table(table(m[, 2]))
        tname = names(ta)
        tname = as.integer(tname)
        r2 = c(r2, rep(tname, ta))
                
        i = i+1
        
        # number of nodes
        gn[i, 1] = length(node)
        # number of genes
        gn[i, 2] = length(unique(m[, 2]))
        
    }
    
    op = par(no.readonly = TRUE)
    par(mfrow=c(1,3))
    t1 = table(r1)
    plot((as.integer(names(t1))), (as.vector(t1)), pch=16, cex=0.8, xlab="Number of member genes in each node", ylab="Frequency", log="y", main="(A) Distribution of the number\nof member genes in each node")
    t2 = table(r2)
    plot((as.integer(names(t2))), (as.vector(t2)), pch=16, cex=0.8, xlab="Number of nodes in which a single gene resides", ylab="Frequency", log="y", main="(B) Distribution of the number\nof nodes in which a single gene resides")
    plot(gn[,2], gn[,1], pch=16, cex=0.8, xlim=range(gn), ylim=range(gn), xlab="Count of genes", ylab="Count of nodes", main="(C) Relationship between node count\nand gene count in pathways")
    abline(a=0,b=1)
    
    par(op)
    
}


# import_biopax = function(data) {
	
# 	if(!inherits(data, "biopax")) {
# 		message("importing biopax")
# 		suppressMessages(biopax <- readBiopax(data, verbose = FALSE))
# 	} else {
# 		biopax = data
# 	}
	
# 	if(biopax$biopaxlevel != 2) {
# 		stop("Only Biopax level 2 is supported.")
# 	}
	
# 	pathway_df = listPathways(biopax)
# 	pathway_list = sapply(pathway_df[[1]], function(pid) {
# 		suppressWarnings(graph <- pathway2RegulatoryGraph(biopax, pid, expandSubpathways = TRUE,
# 		    splitComplexMolecules = FALSE, useIDasNodenames = TRUE, verbose = FALSE))
		
# 		if(!is.null(graph)) {
# 			if(length(edges(graph)) == 0) {
# 				graph = NULL
# 			} else {
				
# 				edge = edges(graph)
# 				input = rep(names(edge), sapply(edge, length))
# 				output = unlist(edge)
# 				interaction_id = paste(pid, seq_along(output), sep = "_")
# 				graph = data.frame(interaction.id = interaction_id, input = input, output = output, stringsAsFactors = FALSE)
# 			}
# 		}
# 		return(graph)
# 	})
# 	pathway_list = pathway_list[!sapply(pathway_list, is.null)]
# 	pathList = lapply(pathway_list, function(pathway) pathway[[1]])
# 	interactionList = do.call("rbind", pathway_list)
	
# 	# nodes in pathway_list are complex ids
# 	all_nodes = c(interactionList[[2]], interactionList[[3]])
# 	all_nodes = unique(all_nodes)
		
# 	# mapping from nodes to molecules
# 	bp2 = selectInstances(biopax, id = all_nodes)
# 	l = isOfClass(bp2, "complex")
# 	complex_nodes = unique(bp2[l]$id)
# 	non_complex_nodes = all_nodes[! all_nodes %in% complex_nodes]
	
# 	nl = c(lapply(complex_nodes, function(nid) {
# 				splitComplex(biopax, nid, returnIDonly = TRUE)
# 			}),
# 			lapply(non_complex_nodes, function(nid) {
# 				nid
# 			}))
# 	names(nl) = c(complex_nodes, non_complex_nodes)
# 	node_list = data.frame(node.id = rep(names(nl), sapply(nl, length)),
# 	                       molecule.id = unlist(nl),
# 						   stringsAsFactors = FALSE)
	
# 	mapping = data.frame(node.id = node_list$node.id,
# 	                     name = sapply(node_list$molecule.id, function(nid) getInstanceProperty(biopax, nid, property = "NAME")),
# 						 class = sapply(node_list$molecule.id, function(nid) getInstanceClass(biopax, nid)),
# 						 stringsAsFactors = FALSE)
# 	mapping = mapping[mapping$class == "protein", 1:2]
	
# 	node.type = sapply(all_nodes, function(nid) getInstanceClass(biopax, nid))
# 	node.name = sapply(all_nodes, function(nid) getInstanceProperty(biopax, nid, property = "NAME"))
	
# 	res = list(pathList = pathList,
# 	           interactionList = interactionList,
# 			   mapping = mapping,
# 			   node.type = node.type,
# 			   node.name = node.name)
# 	class(res) = "pathway.catalogue"
# 	return(res)
# }
