# caculate the centrality of nodes in pathways
# graph: igrpah object
# method: string or function or name
centrality = function(graph, method = "equal.weight") {
    
    if(length(method) > 1) {
        stop("Length of method must be equal to 1.\n") 
    }
    
    if(is.function(method)) {
        return(method(graph))
    } else if(mode(method) == "name") {
        method = eval(method)
        return(method(graph))
    } else if(method == "equal.weight") {
        return(rep(1, vcount(graph)))
    } else if(method == "in.degree") {
        return(degree(graph, mode="in"))
    } else if(method == "out.degree") {
        return(degree(graph, mode="out"))
    } else if(method == "degree") {
        return(degree(graph))
    } else if(method == "betweenness") {
        return(betweenness(graph))
    } else if(method == "in.reach") {
        return(reach(graph, mode="in"))
    } else if(method == "out.reach") {
        return(reach(graph, mode="out"))
    } else if(method == "reach") {
        return(reach(graph))
    } else if(method == "in.spread") {
        return(spread(graph, mode="in"))
    } else if(method == "out.spread") {
        return(spread(graph, mode="out"))
    } else if(method == "spread") {
        return(spread(graph))
    } else {
        stop("`method` should be centrality function or pre-defined centrality measurements.\n")
    }
}

# == title
# Calculate radiality centrality
#
# == param
# -graph an `igraph::igraph` object
# -mode mode of the centrality
# -weights If edges in the graph have weight, then by default, the weight
#         is used to calculate the length of the shortest path. Set it to NULL to supress
#         the weight
# -f function for the weaken rate
#
# == details
# The spread centrality measures how wide the node can send or receive the information in the network.
# Like the water wave, the effect would be weakened with the increase of the distance to other nodes.

# If the weaken function is defined as ``1/x``, then the spread centrality is calculated as
# ``sum(1/d(w, v))`` where ``d(w, v)`` is the length of the shortest path of node ``w`` and node ``v``.
# 
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `reach`, `radiality`
#
# == example
# require(igraph)
# pathway = barabasi.game(200)
# spread(pathway)
spread = function(graph, mode = c("all", "in", "out"),
                  weights = E(graph)$weight, f = function(x) 1/x) {
    mode = match.arg(mode)[1]
    
    sp = shortest.paths(graph, mode = mode, weights = weights)
    s = apply(sp, 1, function(x) {
            return(sum(f(x[x > 0])))
        })
    return(s)
}

# == title
# Calculate largest reach centrality
#
# == param
# -graph an `igraph::igraph` object
# -mode mode of the centrality
# -weights If the edges in the graph have weight, then by default, the weight
#          is used to calculate the length of the shortest path. Set it to NULL to supress
#          the weight.
#
# == details
# The largest reach centrality measures how far a node can send or receive the information in the network.
# It is defined as the largest length of the shortest path from all the other nodes in the network. 

# The largest reach centrality is calculated as ``max(d(w, v))`` where ``d(w, v)`` is the
# length of the shortest path from node ``w`` to node ``v``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `radiality`, `spread`
#
# == example
# require(igraph)
# pathway = barabasi.game(200)
# reach(pathway)
reach = function(graph, weights=E(graph)$weight, mode=c("all", "in", "out")) {
    mode = match.arg(mode)[1]
    
    sp = shortest.paths(graph, weights = weights, mode = mode)
    s = apply(sp, 1, function(x) {
            if(all(x == Inf)) {
                return(0)
            }
            else {
                return(max(x[x != Inf]))
            }
        })
    return(s)
}

# == title
# Calculate radiality centrality
#
# == param
# -graph an `igraph::igraph` object
# -mode mode of the centrality
#
# == details
# The radiality is defined as ``sum(d_G + 1 - d(v, w))/(n - 1)``. where ``d(w, v)`` is the
# length of the shortest path from node ``w`` to node ``v``, ``d_G`` is the diameter of the network,
# n is the size of the network.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(igraph)
# pathway = barabasi.game(200)
# radiality(pathway)
radiality = function(graph, mode = c("all", "in", "out")) {
    mode = match.arg(mode)[1]
    
    sp = shortest.paths(graph, mode = mode, weights = NA)
    n = vcount(graph)
    diam = diameter(graph, directed = ifelse(mode == "all", FALSE, TRUE))
    s = apply(sp, 1, function(x) {
            if(all(x == Inf)) {
                return(0)
            }
            else {
                return(sum(diam+1-x[x != Inf])/(n-1))
            }
        })
    return(s)
}

