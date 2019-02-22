
getCepaORA <- function(input) {
    
    
}
 
#' Apply centrality-extented ORA on a list of pathways
#'
#' @param dif differential gene list
#' @param bk background gene list. If background gene list are not specified, use whole human genes
#' @param pmap.path path to directory containing igraph objects of individual pathways
#' @param pathway.name character vector of pathways
#' @param cen centrality measuments, it can ce a string, or a function
#' @param cen.name centrality measurement names. By default it is parsed from ``cen`` argument
#' @param iter number of simulations; should be greater than 100
#
# == details
# The traditional over-representation analysis (ORA) to find significant pathways 
# uses a 2x2 contingency table to test the independency of genes belonging to a 
# functional category and these genes being differentially expressed, usually by 
# Fisher's exact test. The ORA only consider the number of genes and the function
# extend traditional ORA with network centralities.
#
# The differential gene list and the background gene list should be indicated
# with the same identifiers (e.g. gene symbol or refseq ID). All genes in
# the differential gene list should exist in the background gene list. If users 
# use the `PID.db` data, all genes should be formatted in gene symbol.
#
# If the centrality measurement is set as a string, only pre-defined "equal.weight",
# "in.degree", "out.degree", "degree", "betweenness", "in.reach", "out.reach",
# "reach", "in.spread", "out.spread" and "spread" are allowed. More centrality
# measurements can be used by setting it as a function (such as closeness,
# cluster coefficient). In the function, we recommand users choose
# at least two centrality measurements. The default centralities are "equal.weight",
# "in.degree", "out.degree", "betweenness", "in.reach" and "out.reach".
#
# However, in most circumstance, the function is called by `cepa.all`.
#
# == value
# A `cepa.all` class object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == example
# \dontrun{
# data(PID.db)
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa.ora.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
# }
cepa.ora.metab.all = function(dif,bk,pmap.path,pathway.name, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    iter = 1000) {
    
    dif = dif[dif %in% bk]
    
    if(length(cen) < 1) {
        stop("cen argument must be specified.\n")
    }
    
    # if cen argument is a function, the function should be quoted or substituted
    # because we need the function name
    for(ce in cen) {
        if(is.function(ce)) {
            stop("Functions cannot be used directly because we need the function name, use quote or substitute.\n")
        }
    }

    n.pathway = length(pathway.name)
    
    pathway.result = list()
    length(pathway.result) = n.pathway
    # pathway.result is like a two layer list
    pathway.result = lapply(pathway.result, function(x) {
                              y = list()
                              length(y) = length(cen.name)
                              names(y) = cen.name
                              return(y)
                            })
    names(pathway.result) = pathway.name
    
    cat("  Calculate pathway scores...\n")
    for(i in 1:n.pathway) {
        
        cat("    ", i, "/", n.pathway, ", ", pathway.name[i], "...\n", sep="")
        
#        pathway = generate.pathway(as.matrix(inter))
        pathway = getPathwayIgraph(Pathway.Name=pathway.name[i], pmap.path=pmap.path)
        
        j = 0
        # to this pathway, applying various centralities
        for(ce in cen) {
            j = j + 1
            pathway.result[[i]][[j]] = cepa.ora.metab(dif = dif, bk = bk, pathway = pathway, cen = ce, iter = iter)
            cat("      - ", ce, ": ", round(pathway.result[[i]][[j]]$p.value, 3), "\n", sep = "")
        }
    }

    class(pathway.result) = "cepa.all"
    return(pathway.result)
}


#' Apply centrality-extended ORA on a single pathway
#'
#'
#' @param dif differential gene list
#' @param bk background gene list. If background gene list are not specified, use whole human genes
#' @param pathway `igraph::igraphtest` object or edge list
#' @param cen centrality measuments, it can ce a string, function, or function that has been quoted
#' @param cen.name centrality measurement names. This argument should be set if the ``cen`` is a 
#' function.
#' @param iter number of simulations.  Should be >= 100.
#' 
#' @import igraph
#' @importFrom data.table data.table
#
# == details
# The function is always called by `cepa.ora.all`. But you can still
# use it if you realy want to analysis just one pathway under one centrality.
#
# == value
# A ``cepa`` class object
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
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)

# }
cepa.ora.metab = function(dif, bk, pathway = NULL, cen = "equal.weight",
                cen.name = if(is.function(cen)) deparse(substitute(cen)) 
                           else if(mode(cen) == "name") deparse(cen) 
                           else cen,
                iter = 1000) {
    

    # some checking of the arguments
    if(length(dif) > length(bk)) {
        stop("Length of differential genes should not be larger than the length of background genes.\n")
    }
    if(sum(dif %in% bk) != length(dif)) {
        stop("Differential genes must be all in background list.\n")
    }
    
    # you can specify a pathway object or a pathway id
    if(! is.null(pathway)) {   # a pathway is specified
       if (class(pathway) != "igraph") {    # it should be an igraph object
            stop("Since pathway is not formatted as edge list, it should be an igraph object.")
        }
    }  else {  # one of pathway and id should be set
        stop("You should specify pathway argument.")
    }
    
    if(iter < 100) {
        stop("Iterations should not be smaller than 100.\n")
    }
    
    # single centrality!
    if(length(cen) != 1) {
        stop("Length of cen must be equal to 1.\n") 
    }
    
    ### PSchange when we reconfigure our code as package ####
    weight = CePa:::centrality(pathway, cen)
    if (is.null(names(weight))) names(weight) <- pathway.nodes(pathway)
    
    
    # is our algorithm, only non-negative centrality is allowed
    if(any(weight < 0)) {
        stop("Weight should not be negative.")
    }

    add = 0
    # if there are none-zero weight values
    if(any(weight > 0)) {
        add = min(weight[weight > 0])/100
    }
    weight = weight + add
    
    # nodes in the pathway
    node = pathway.nodes(pathway)
    
    # only the mapping in the pathway
#    mapping = pc$mapping[pc$mapping[, 1] %in% node, ]
    mapping = data.frame(node=node, names=get.vertex.attribute(pathway, "label"), shape=get.vertex.attribute(pathway, "shape"))
    mapping$weight <- weight[mapping$node]
    
    dt <- data.table(mapping)
    dt <- dt[,.I[which.max(weight)],by=.(names, shape)]
    mapping <- mapping[dt$V1,]
    
       
    # get node names formatted with genes
    node.name = mapping$node
    node <- mapping$node
    weight <- mapping$weight

    member = character(0)
 
    for(i in 1:length(node)) {
        # genes that exsit in the node
        l = mapping[, 1] == node[i]

        # if find nodes with genes mapped
        if(sum(l)) {
            member = sort(unique(mapping[l, 2]))
            # mark the diff genes
            member[member %in% dif] = paste("[", member[member %in% dif], "]", sep="")
            node.name[i] = paste(member, collapse = "\n")
        }
    }
    
 

    # map dif genes to node id
    dif.node = unique(mapping[mapping[, 2] %in% dif, 1])
    is.dif.node = as.numeric(node %in% dif.node)
    
    # get weights
    node.level = is.dif.node * weight
    s = sum(node.level)

    # distribution of node level value (combined with weight)
    if(sum(is.dif.node > 0) == 0) {   # if there is no differential genes
        ds = c(0, 0, 0, 0)
    } else {
        ds = quantile(node.level, c(1, 0.75, 0.5, 0))
    }
    names(ds) = c("max", "q75", "median", "min")

    # sampling
    p.dif = length(dif) / length(bk)  # probability to be a differential gene
    s.random = numeric(iter)
    ds.random = matrix(0, iter, 4)   # descriptive statistic of the node
    colnames(ds.random) = c("max", "q75", "median", "min")

    # genes in the pathway
    gene = unique(mapping[mapping[, 1] %in% node, 2])

    for(i in 1:iter) {
        # simulated random genes
        dif.gene.random = gene[as.logical(rbinom(length(gene), 1, p.dif))]

        # then map to node id
        dif.node.random = unique(mapping[mapping[, 2] %in% dif.gene.random, 1])
        # find which node is differentially affected
        is.dif.node.random = as.numeric(node %in% dif.node.random)
        # calculate the score
        node.level.random = is.dif.node.random * weight
        s.random[i] = sum(node.level.random)
        if(sum(is.dif.node.random > 0) == 0) {
            ds.random[i, ] = c(0, 0, 0, 0)
        }
        else {
            ds.random[i, ] = quantile(node.level.random, c(1, 0.75, 0.5, 0))
        }
    }

    p.value = (sum(s.random >= s) + 1) / (iter + 1)

    dif.gene = intersect(dif, gene)
    n.dif.node = length(dif.node)
    n.dif.gene = length(dif.gene)
    n.node = length(node)
    n.gene = length(gene)

    count = c(n.dif.node, n.node, n.dif.gene, n.gene)
    names(count) = c("n.dif.node", "n.node", "n.dif.gene", "n.gene")


    res = list("score" = s,                                  # pathway score
               "score.distribution" = ds,                    # distribution of node value in the pathway
               "score.random" = s.random,                    # simulated pathway scores
               "score.distribution.random" = ds.random,      # distribution of node value in the pathway in each simulation
               "p.value" = p.value,                          # p value
               "centrality" = cen.name,                      # centrality name
               "weight" = weight,                            # weight for each node
               "node.level.t.value" = as.integer(is.dif.node),    # value for each node, exclude the centrality part
               "node.level" = node.level,                    # value for each node, exclude the centrality part
               "node.name" = node.name,                      # node names
               "pathway" = pathway,                          # pathway in igraph format
               "count" = count,
               "framework" = "ora")

    class(res) = "cepa"

    return(invisible(res))

}

