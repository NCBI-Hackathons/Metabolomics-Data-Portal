# == title
# plot the cepa.all object
#
# == param
# -x          a `cepa.all` object
# -id         index or the name for the pathway
# -cen        index or the name for the centrality
# -type       If the aim is to plot single pathway, then this argument is to identify the kind of the plotting.
# -tool       Use which tool to visualize the graph. Choices are 'igraph' and 'Rgraphviz'
# -node.name  node.name for each node
# -node.type  node.type for each node
# -adj.method method of `stats::p.adjust`, available methods are "holm", "hochberg", 
#             "hommel", "bonferroni", "BH", "BY", "fdr", "none"
# -only.sig   whether to show all pathways. If just show significant pathways, 
#             the names for each significant pathway will be draw.
# -cutoff     cutoff for significance
# -...        other arguments
#
# == details
# This function has two applications. First, it can draw heatmaps of p-values
# of all pathways under different centrality measurements. To do it, users should set
# ``x``, ``adj.method``, ``only.sig``, ``cutoff`` arguments. 
#  
# Second, it can draw figures of single
# pathway under specific centrality measurement. Under this circumstance, 
# this function is just a wrapper of `plot.cepa`. To do it, 
# users should set ``x``, ``id``, ``cen``, ``type``, ``tool``, ``node.name`` and ``node.type`` arguments. The
# ``id`` and ``cen`` arguments are used to get single `cepa` object that sent to the 
# plot function.
# 
# It must be noted that these two kinds of arguments should not be mixed.
# 
# There is also another popular method ``qvalue`` to adjust p-values. However, errors
# may occur when adjusting some kind of p-value list by ``qvalue``.
# So ``qvalue`` was not implemented into CePa. But still users can override the default
# p.adjust to support qvalue by themselves, see the vignette.
# 
# == seealso
# `cepa.all`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# data(PID.db)
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
# plot(res.ora)
# plot(res.ora, id = 3)
# plot(res.ora, id = 3, type = "null")
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# plot(res.gsa)
# plot(res.gsa, id = 3, cen = 2)
# plot(res.gsa, id = 3, cen = 2, type = "null")
# }
plot.cepa.all = function(x, id = NULL, cen = 1, type = c("graph", "null"), tool = c("igraph", "Rgraphviz"),
                        node.name = NULL, node.type = NULL,
                        adj.method = "none", only.sig = FALSE,
                        cutoff = ifelse(adj.method == "none", 0.01, 0.05), ...) {
    
    if(class(x) != "cepa.all") {
        stop("x should be top.all object.\n")
    }
    
    if(!is.null(id)) {
        if(length(cen) > 1) {
            stop("Length of cen must be equal to 1.\n")
        }
        
        if(is.function(cen)) {
            cen = deparse(substitute(cen))
        }
        else if(mode(cen) == "name") {
            cen = deparse(cen)
        }
        
        plot(get.cepa(x, id, cen), type = type, tool = tool, node.name = node.name, node.type = node.type, ...)

    } else {    
        p.heatmap(x, adj.method = adj.method, only.sig = only.sig, cutoff = cutoff)
    }
}


# #
# # ###########################################################
# # # 
# # #     write a new p.adjust to supprot qvalue
# # #
# library(qvalue)
# p.adjust = function(p, method = c("holm", "hochberg", "hommel", "bonferroni",
#                     "BH", "BY", "fdr", "none", "qvalue"), ...) {
#     if(method == "qvalue") {
#         # qvalue has more arguments, pass them by ...
#         qvalue(p, ...)$qvalue
#     } else {
#         stats::p.adjust(p, method)
#     }
# }
# }


p.heatmap = function(x, adj.method = "none", only.sig = TRUE,
                     cutoff = ifelse(adj.method == "none", 0.01, 0.05)) {

    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    p.value = p.table(x)
    p.foo = apply(p.value, 2, p.adjust, adj.method)
    if(!is.matrix(p.foo)) {
        p.foo = matrix(p.foo, nrow = 1)
        rownames(p.foo) = rownames(p.value)
        colnames(p.foo) = colnames(p.value)
        p.value = p.foo
    }
    p.value = p.foo
    
    np = ncol(p.value)
    
    if(only.sig) {
        l = apply(p.value, 1, function(x) { sum(x <= cutoff) > 0})
        p.value = p.value[l, , drop=FALSE]
        if(sum(l) == 0) {
            plot(1,1, axes=FALSE, ann=FALSE, type="n")
            text(1,1,"No significant pathway!", cex=2)
            return(invisible(NULL))
        }
    }

    # get the order of p.fisher
    o = order(apply(-log(p.value + 1e-8), 1, mean), decreasing = TRUE)
    
    p.smallest = min(p.value) + 1e-4
    p.value[p.value == 0] = p.smallest
    p2 = -log10(p.value)
    p2 = p2[o,,drop=FALSE]
    
    # smallest p value is 0.001
    #p2[p2 > floor(-log10(p.smallest))] = floor(-log10(p.smallest))
    
    p2 = t(p2)
    
    if(only.sig) {  # do not draw pathway names
        cn = colnames(p2)
        max.length = max(nchar(cn))
        layout(rbind(c(1, 2),
                     c(3, 4),
                     c(5, 0)),
               widths=c(7,lcm(3.5)),
               heights=c(1.5,5,max.length/4))
    } else {
        layout(rbind(c(1, 2),
                     c(3, 4)),
               widths=c(7,lcm(3.5)),
               heights=c(1.5,5))
    }
    
    # 1.only draw titles
    par("mar" = c(0, 0, 0, 0))
    plot(0, 0, type="n", axes=FALSE, ann=FALSE)
    if(only.sig) {
        text(0, 0, paste("Heatmap of ", ifelse(adj.method=="none", "p-value", "FDR"),"s of pathways (only significant)", sep=""), cex=1.5)
    } else {
        text(0, 0, paste("Heatmap of ", ifelse(adj.method=="none", "p-value", "FDR"),"s of pathways", sep=""), cex=1.5)
    }
    
    fc = c(0, -log10(cutoff), floor(-log10(p.smallest))+1)
    colors = c("green", "white", "red")
    
    # 2. legend
    g = seq(0, floor(-log10(p.smallest))+1, length.out=100)
    lg = length(g)
    plot(c(0,0), c(lg, 3), type="n", xlim=c(0, lg*1.1), ylim=c(0, 3), axes=FALSE, ann=FALSE)
    text(1, 1.5, 1)
    text(lg, 1.5, 10^(-fc[3]))
    text((fc[2]-fc[1])/(fc[3]-fc[1])*(lg-1)+1, 1.5, 10^(-fc[2]))
    for (i in 1:lg) {
        rect(i-1, 0, i, 1, col=get_color(g[i], colors=colors, fc=fc), border=NA)
    } 
    
    # 3. heatmap
    par("mar" = c(0, 0, 0, 0))
    ncol = dim(p2)[2]
    nrow = dim(p2)[1]
    plot(c(0,0), c(ncol, nrow), type="n", xlim=c(0, ncol), ylim=c(0, nrow), axes=FALSE, ann=FALSE, xaxs="i")
    
    for (i in 1:nrow) {
        for (j in 1:ncol) {
            rect(j-1, i-1, j, i, col=get_color(p2[i, j], colors=colors, fc=fc), border=NA)
        }
    }
    
    # 4. rownames
    rn = rownames(p2)
    par("mar" = c(0, 0, 0, 2))
    plot(c(0,0), c(1, nrow), type="n", xlim=c(0, 1), ylim=c(0, nrow), axes=FALSE, ann=FALSE)
    for (i in 1:nrow) {
        text(0, i-0.5, rn[i], adj=c(0, 0.5), cex=ifelse(only.sig, 1.5, 1.2))
    }
    
    # 5. colnames
    if(only.sig) {
        cn = colnames(p2)
        par("mar" = c(2, 0, 0, 0))
        plot(c(0,0), c(ncol, 1), type="n", xlim=c(0, ncol), ylim=c(0, 1), axes=FALSE, ann=FALSE, xaxs="i")
        for (i in 1:ncol) {
            text(i-0.5, 1, cn[i], srt=-90, adj=c(0, 0.5), cex=1.5)
        }
    }
    
    par("mar" = c(5.1, 4.1, 4.1, 2.1))
    layout(matrix(1, 1, 1))
    return(invisible(NULL))
}


# get color code from data
get_color = function(x,
                     colors = c("green", "black", "red"),
                     fc = c(-5, 0, 5),
                     gradient = function(x)x, ...) {
                     
    col_MIN = as.vector(col2rgb(colors[1]))
    col_MEDIAN = as.vector(col2rgb(colors[2]))
    col_MAX = as.vector(col2rgb(colors[3]))
    
    fc = sign(fc)*gradient(abs(fc))
    
    color = character(length(x))
    for(i in 1:length(x)) {
        if(!is.numeric(x[i])) {
            color[i] = rgb(128,12,128, maxColorValue = 255)
            next
        }
        value = sign(x[i])*(gradient(abs(x[i])))
        if(x[i] <= fc[2]) {
            col_num = (value - fc[2])*(col_MEDIAN - col_MIN)/(fc[2] - fc[1]) + col_MEDIAN
        }
        else {
            col_num = (value - fc[2])*(col_MEDIAN - col_MAX)/(fc[2] - fc[3]) + col_MEDIAN
        }
        
        col_num = as.integer(col_num)
        col_num[col_num > 255] = 255
        col_num[col_num < 0] = 0
        
        color[i] = rgb(col_num[1], col_num[2], col_num[3], maxColorValue = 255, ...)
    }
    
    return(color)
}

# == title
# print the cepa.all object
# 
# == param
# -x   a `cepa.all` object
# -... other arguments
#
# == details
# The function print the number of significant pathways under various
# centrality measures at p-value <= 0.01.
#
# == seealso
# `cepa.all`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# data(PID.db)
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
# res.ora
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# res.gsa
# }
print.cepa.all = function(x, ...) {
    cat("\n")
    cat("number of pathways:", length(x), "\n")
    cat("\n")
    
    p.value = p.table(x)
    centrality = colnames(p.value)

    cat("Significant pathways (p.value <= 0.01):\n")
    sig = matrix(0, nrow=length(centrality), 1)
    rownames(sig) = centrality
    colnames(sig) = "Number"
    for(i in 1:length(centrality)) {
        sig[i, 1] = sum(p.value[ ,i] <= 0.01)
    }
    print(sig)
    
    cat("\n")
}

# == title
# Plot the cepa object
#
# == param
# -x    a `cepa` object
# -type identify the type for the plot
# -...   arguments passed to `plotGraph`
#
# == details
# The function is wrapper of `plotGraph` and `plotNull`.
# If type is specified to "graph", the graph of the network will be plotted (see `plotGraph` for details).
# If type is specified to "null", the null distribution of the pathway score
# in the pathway will be plotted (see `plotNull` for details).
#
# == return
# if type is set to "graph", the function will return a `igraph::igraph` object or a ``graphML`` object of the pathway. Else it is NULL.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `cepa`, `plotNull`, `plotGraph`
#
# == example
# \dontrun{
#
# data(PID.db)
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)
# plot(res.ora)
# plot(res.ora, type = "null")
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa(mat = eset, label = label, pc = PID.db$NCI, id = 2)
# plot(res.gsa, type = "null")
# }
plot.cepa = function(x, type = c("graph", "null"), ...) {
    type = type[1]
    if(type == "graph") {
        plotGraph(x, ...)
    } else if(type == "null") {
        plotNull(x)
    }
}

# == title
# Plot the null distribution of the pathway score
#
# == param
# -x a `cepa` object
#
# == details
# There are two figures in the plotting.
#
# - A) Distribution of node score in the pathway under simulation. Since a pathway contains
# a list of nodes. The distribution of node score for the pathway in each simulation is measures by maximum value, 
# the 75th quantile, median value and minimum value. The distribution of node score for 
# the pathway in the real data is highlighted.
#
# - B) Histogram of simulated pathway scores. 
#
# The function is always called through `plot.cepa.all` and `plot.cepa`.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `cepa`, `plot.cepa`
#
# == example
# \dontrun{
# data(PID.db)
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
# ora = get.cepa(res.ora, id = 5, cen = 3)
# plotNull(ora)
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# gsa = get.cepa(res.gsa, id = 5, cen = 3)
# plotNull(gsa)
# }
plotNull = function(x) {
    s = x$score
    s.random = x$score.random
    ds = x$score.distribution
    ds.random = x$score.distribution.random
    weight = x$weight
    centrality = x$centrality
    p.value = x$p.value
    pathway = x$pathway              # igraph object
    
    n.node = length(weight)

    
    if(length(weight) == 0) {
        plot(1,1, axes=FALSE, ann=FALSE, type="n")
        text(1,1,"empty pathway!", cex=3)
        return(NULL)
    }

    opar = par(no.readonly = TRUE)

    
    # color settigs
    color = character(4)
    names(color) = colnames(ds.random)
    color["max"] = rgb(227, 26, 28, maxColorValue = 255)
    color["q75"] = rgb(255, 191, 111, maxColorValue = 255)
    color["median"] = rgb(127, 201, 127, maxColorValue = 255)
    color["min"] = rgb(31, 120, 180, maxColorValue = 255)
    
    par(mfrow = c(1, 2))
    
    # figure A
    xrange = range(c(s.random, s))
    yrange = range(c(as.vector(ds.random), ds))
    if(yrange[1] == yrange[[2]]) {
        yrange = sort(c(0, 2*yrange[1]))
    }
    matplot(jitt(s.random, xrange), jitt(ds.random, yrange), xlim = xrange, ylim = yrange, col=color, pch=20, cex=0.2, xlab = "Simulated score", ylab=ifelse(x$framework == "ora.univariate", "Centrality", "Node-level scores"), main = ifelse(x$framework == "ora", paste("(A) Distribution of", centrality, "of differential nodes\nin pathway under simulation"), paste("(A) Distribution of node-level scores\nweighted by", centrality, "centrality")))
    points(rep(s, 4), jitt(ds, yrange), col=color, pch=20, cex=5)
    legend(min(xrange), max(yrange), colnames(ds.random), col=color, pch=20, pt.cex=2)
    abline(v = s, col="red", lwd = 2)

    # figure C
    hist(s.random, freq = FALSE, xlim = range(c(s.random, s)), breaks = 50, xlab = "Simulated score", main = paste("(B) Histogram of simulated scores in the pathway\nusing", centrality, "centrality as weight"))
    box()
    abline(v = s, col = "red", lwd = 2)
    
    par(opar)
}

# == title
# Plot graph for the pathway network
#
# == param
# -x                   a `cepa` object
# -node.name           node.name for each node
# -node.type           node.type for each node
# -draw                Whether to draw the graph
# -tool                Use which tool to visualize the graph. Choices are 'igraph' and 'Rgraphviz'
# -graph.node.max.size max size of the node in the graph
# -graph.node.min.size min size of the node in the graph
# -graph.layout.method function of the layout method. For the list
#                      of available methods, see `igraph::layout`
# -...                 other arguments
#
# == details
# Graph view of the pathway where the size of node is proportional to centrality 
# value of the node. 
#
# By default, the layout for the pathway tree-like. If the number of pathway nodes
# is large, the layout would be a random layout.
#
# Two packages can be selected to visualize the graph: ``igraph`` and ``Rgraphviz``.
# Default package is ``igraph`` (in fact, this package just uses the data generated from
# the layout function in `igraph::igraph` package, which are the coordinate of nodes and edges.
# And the I re-wrote the plotting function to generate the graph). From my personal view,
# ``Rgraphviz`` package generated more beautiful graphs.
# 
# If the ``tool`` is set as ``igraph``, the function returns a `igraph::igraph` object. And
# if the ``tool`` is set as ``Rgraphviz``, the function returns a ``graphAM`` class object.
# So if users don't satisfy, they can draw graphs of the network with their
# own settings.
#
# The function is always called through `plot.cepa.all` and `plot.cepa`.
#
# == return
# A `igraph::igraph` object of the pathway
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
# res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
# ora = get.cepa(res.ora, id = 5, cen = 3)
# plotGraph(ora)
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# gsa = get.cepa(res.gsa, id = 5, cen = 3)
# plotGraph(gsa)
# }
plotGraph = function(x, node.name = NULL, node.type = NULL, draw = TRUE, 
    tool = c("igraph", "Rgraphviz"), graph.node.max.size = 20, 
    graph.node.min.size = 3, graph.layout.method = NULL) {
    
    tool = tool[1]
    if(! tool %in% c("igraph", "Rgraphviz")) {
        stop("Must choose visuallization method in 'igraph' and 'Rgraphviz'.")
    }
    
    s = x$score
    s.random = x$score.random
    weight = x$weight
    centrality = x$centrality
    p.value = x$p.value
    pathway = x$pathway              # igraph object
    node.level = x$node.level.t.value
    
    n.node = length(node.level)
    
    node.id = pathway.nodes(pathway)
    # check node.name and node.type
    if(! is.null(node.name)) {
        node.name = node.name[node.id]
        node.name[is.na(node.name)] = ""
    } else {
        node.name = x$node.name   
    }
    
    if(! is.null(node.type)) {
        node.type = node.type[node.id]
        node.type[is.na(node.type)] = "asist.node"
    }
    
    if(length(weight) == 0) {
        plot(1,1, axes=FALSE, ann=FALSE, type="n")
        text(1,1,"empty pathway!", cex=3)
        return(NULL)
    }
    
    
    # node should only be in protein, complex, family, compound, rna, subunit
    if(is.null(node.type)) {
        if(x$framework == "gsa.univariate") {
            v.color = get_color(node.level, colors = c("#91CF60", "#EEEEEE", "#FC8D59"), fc = c(-3, 0, 3), gradient = function(x) x^2)
        } else {
            v.color = get_color(node.level, colors = c("#4DAF4A", "#4DAF4A", "#E41A1C"), fc = c(-1, 0, 1), gradient = function(x) x^2)
        }
        v.color[node.name == ""] = "#377EB8" # non-protein nodes
        
        v.frame.color = "white"
    } else {    
        if(x$framework == "gsa.univariate") {
            v.color = get_color(node.level, colors = c("#91CF60", "#EEEEEE", "#FC8D59"), fc = c(-3, 0, 3), gradient = function(x) x^2)
        } else {
            v.color = get_color(node.level, colors = c("#4DAF4A", "#4DAF4A", "#E41A1C"), fc = c(-1, 0, 1), gradient = function(x) x^2)
        }
        
        all.node.types = unique(node.type)
        other.node.types = all.node.types[!grepl("^(protein|complex|asist.node)$", all.node.types, ignore.case=TRUE)]
        other.colors = c("#377EB8", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
        for(i in seq_along(other.node.types)) {
            v.color[node.type == other.node.types[i]] = other.colors[i]   
        }
    }
    
    # size
    max.size = graph.node.max.size
    min.size = graph.node.min.size
    v.size = rep(min.size, n.node)
    if(max(weight) != min(weight)) {
        v.size = min.size + (max.size - min.size)/(max(weight) - min(weight))*(weight - min(weight))
    }
    
    # label
    if(is.null(node.name)) {
        v.label = x$node.name
    }    else {
        v.label = node.name
        v.label = gsub("/", "\n", v.label)
    }
    
    if(! is.null(node.type)) {
        v.size[node.type == "asist.node"] = min.size
        v.color[node.type == "asist.node"] = "#CCCCCC"
    }

    V(pathway)$shape = "circle"
    V(pathway)$color = v.color
    V(pathway)$size = v.size
    V(pathway)$label = v.label
    V(pathway)$label.cex = 0.8
    V(pathway)$label.color = "black"
    E(pathway)$arrow.size = 0.5
    #E(pathway)$color = "#CCCCCC"
    
    if(! draw) {
        if(tool == "igraph") {
            return(pathway)
        } else if(tool == "Rgraphviz") {
            return(plotByRGraphviz(pathway, draw = FALSE))
        }
    }
    
    if(is.null(graph.layout.method)) {
        layout.method = function(x) layout.reingold.tilford(x, mode = "all")
    }    else {
        layout.method = graph.layout.method
    }
    
    opar = par(c("mar", "xpd"))
    
    layout(cbind(1:2, 3:4), heights = c(9,cm(1)), widths = c(9, cm(0)))
    if(tool == "igraph") {
        plot.igraph2(pathway, layout.method=layout.method)
    } else if(tool == "Rgraphviz") {
        ra = plotByRGraphviz(pathway)
    }
    
    par(mar = c(1, 4, 0, 4), xpd=FALSE)
    plot(0,0, xlim = c(0, 1), type = "n", ylim = c(0, 1), ann = FALSE, axes = FALSE, xaxs = "i", yaxs = "i")
    if(is.null(node.type)) {
        if(x$framework == "ora") {
            legend(0, 0, c("diff nodes", "non-diff nodes", "non-protein nodes", "asist nodes"),
                   col=c("#E41A1C", "#4DAF4A", "#377EB8", "#CCCCCC"), pch=16, pt.cex=2, yjust=0, cex=0.8)
        } else {
            e = seq(-4, 4, by=0.5)
            color = get_color(e, colors = c("#91CF60", "#EEEEEE", "#FC8D59"), fc = c(-3, 0, 3), gradient = function(x) x^2)
            startx = 0.8
            endx = 1
            for(i in 1:length(e)) {
                x_left = startx + (i-1)*(endx - startx)/length(e)
                x_right = startx + (i)*(endx - startx)/length(e)
                rect(x_left, 0.1, x_right, 0.4, col = color[i], border = NA)
            }
            text(startx, 0.1, min(e), cex=0.8, adj=c(0.5, 1))
            text((startx+endx)/2, 0.1, "0", cex=0.8, adj=c(0.5, 1))
            text(endx, 0.1, max(e), cex=0.8, adj=c(0.5, 1))
            text((startx+endx)/2, 0.4, "nodes t-value", cex=0.8, adj=c(0.5, -1))
            #rect(startx-0.04, -0.05, endx+0.04, 0.15, border=TRUE)
            
            legend(0, 0, c("non-protein nodes", "asist nodes"),
                   col=c("#377EB8", "#CCCCCC"), pch=16, pt.cex=2, yjust=0, cex=0.8)
        }
    } else {
        if(x$framework == "ora") {
            all.node.types = unique(node.type)
            other.node.types = all.node.types[!grepl("^(protein|complex|asist.node)$", all.node.types, ignore.case=TRUE)]
            other.colors = c("#377EB8", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
            
            legend(0, 0, c("diff nodes", "non-diff nodes", other.node.types, "asist nodes"),
                   col=c("#E41A1C", "#4DAF4A", other.colors[seq_along(other.node.types)], "#CCCCCC"),
                    pch=16, pt.cex=2, cex=0.8, yjust=0)
        } else {
            e = seq(-4, 4, by=0.5)
            color = get_color(e, colors = c("#91CF60", "#EEEEEE", "#FC8D59"), fc = c(-3, 0, 3), gradient = function(x) x^2)
            startx = 0.8
            endx = 1
            for(i in 1:length(e)) {
                x_left = startx + (i-1)*(endx - startx)/length(e)
                x_right = startx + (i)*(endx - startx)/length(e)
                rect(x_left, 0.1, x_right, 0.4, col = color[i], border = NA)
            }
            text(startx, 0.1, min(e), cex=0.8, adj=c(0.5, 1))
            text((startx+endx)/2, 0.1, "0", cex=0.8, adj=c(0.5, 1))
            text(endx, 0.1, max(e), cex=0.8, adj=c(0.5, 1))
            text((startx+endx)/2, 0.4, "nodes t-value", cex=0.8, adj=c(0.5, -1))
            #rect(startx-0.04, -0.05, endx+0.04, 0.15, border=TRUE)
                
            all.node.types = unique(node.type)
            other.node.types = all.node.types[!grepl("^(protein|complex|asist.node)$", all.node.types, ignore.case=TRUE)]
            other.colors = c("#377EB8", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
            
            legend(0, 0, c(other.node.types, "asist nodes"),
                col=c(other.colors[seq_along(other.node.types)], "#CCCCCC"),
                pch=16, pt.cex=2, cex=0.8, yjust=0)
        }
    }
    layout(cbind(1))
    par(opar)
    
    if(tool == "igraph") {
        return(invisible(pathway))
    } else if(tool == "Rgraphviz") {
        return(invisible(ra))
    }
}

plot.igraph2 = function(g, layout.method = layout.random) {
    
    v.color = V(g)$color
    v.shape = V(g)$shape
    v.shape2 = rep(16, length(v.shape))
    v.shape2[v.shape == "square"] = 15
    v.shape = v.shape2
    v.size = V(g)$size
    names(v.size) = pathway.nodes(g)
    v.label = V(g)$label
    v.frame.color = V(g)$frame.color
    v.label.cex = V(g)$label.cex
    v.label.font = V(g)$label.font
    v.label.color = V(g)$label.color
    e.arrow.size = E(g)$arrow.size
    e.color = E(g)$color
    
    ly = layout.method(g)
    r1 = max(ly[, 1]) - min(ly[, 1])
    if(r1 == 0) {
        ly[, 1] = 0.5
    } else { 
        ly[, 1] = (ly[, 1] - min(ly[, 1]))/r1
    }
    r2 = max(ly[, 2]) - min(ly[, 2])
    if(r2 == 0) {
        ly[, 2] = 0.5
    } else {
        ly[, 2] = (ly[, 2] - min(ly[, 2]))/r2
    }
    rownames(ly) = pathway.nodes(g)
    par(mar = c(0, 4.1, 2, 2.1))
    # points
    plot(ly, cex = v.size,
             col = v.color,
             pch = v.shape,
             ann = FALSE,
             axes = FALSE,
             xlim = c(-0.1, 1.1),
             ylim = c(-0.1, 1.1),
             xaxs = "i",
             yaxs = "i"
             )
    frame.pch = v.shape
    frame.pch[frame.pch == 16] = 21
    frame.pch[frame.pch == 15] = 0
    points(ly, pch = frame.pch, cex = v.size,
                         col = v.frame.color)
    text(ly, v.label, cex = v.label.cex,
                       col = v.label.color)
                                              
    # edges
    el = get.edgelist(g)
    from = ly[el[, 1], ,drop=FALSE]
    from.node = rownames(from)
    to = ly[el[, 2], ,drop=FALSE]
    to.node = rownames(to)
    new.from = matrix(0, nrow=dim(from)[1], ncol=2)
    new.to = matrix(0, nrow=dim(to)[1], ncol=2)
    for(i in 1:length(from.node)) {
        foo = reedge(from[i, 1], from[i, 2], to[i, 1], to[i, 2], v.size[from.node[i]]/21/par("pin")[1], v.size[to.node[i]]/21/par("pin")[1])
        new.from[i, ] = foo[1:2]
        new.to[i, ] = foo[3:4]
    }
    
    segments(new.from[, 1], new.from[, 2], new.to[, 1], new.to[, 2])
    
    for(i in 1:ecount(g)) {
        a.coor = arrow.coor(new.from[i, 1], new.from[i, 2],
                            new.to[i, 1], new.to[i, 2])
        polygon(a.coor, col="black")
    }
    
}



arrow.coor = function(x1, y1, x2, y2, length = 0.02) {
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    st = -(y2 - y1)/d
    ct = (x2 - x1)/d
    A = matrix(c(ct, -st, st, ct), 2)

    M = c(x2, y2)
    p1 = c(-cos(pi/12)*length, sin(pi/12)*length)
    a1 = A %*% p1 + M
    a1 = as.vector(a1)
    p2 = c(-cos(pi/12)*length, -sin(pi/12)*length)
    a2 = A %*% p2 + M
    a2 = as.vector(a2)
    
    return(rbind(a1, a2, c(x2, y2)))
}

reedge = function(x1, y1, x2, y2, r1, r2) {
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    st = -(y2 - y1)/d
    ct = (x2 - x1)/d
    A = matrix(c(ct, -st, st, ct), 2)

    M = c(x1, y1)
    p1 = c(r1, 0)
    a1 = A %*% p1 + M
    a1 = as.vector(a1)
    p2 = c(d - r2, 0)
    a2 = A %*% p2 + M
    a2 = as.vector(a2)
    
    return(c(a1, a2))
}

# == title
#   print the cepa object
#
# == param
# -x   a `cepa` object
# -... other arguments
#
# == details
# The function print procedure of the analysis, the centrality and the p-value for the pathway.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == seealso
# `cepa`
#
# == example
# \dontrun{
#
# data(PID.db)
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, id = 2)
# res.ora
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa(mat = eset, label = label, pc = PID.db$NCI, id = 2)
# res.gsa
# }
print.cepa = function(x, ...) {
    cat("\n")
    cat("  procedure:", x$framework, "\n")
    cat("  weight:", x$centrality, "\n")
    cat("  p-value:", ifelse(x$p.value < 0.001, sprintf("%.3e", x$p.value), sprintf("%.3f", x$p.value)), "\n")
    cat("\n")
}


# add random noise to a point in a xy coordinate system
# similar to function jitter
jitt = function(x, r = range(x[x != Inf], na.rm=TRUE)) {
    return(x + runif(length(x), -0.5, 0.5)*(r[2] - r[1])/50)
}

plotByRGraphviz = function(g, draw = TRUE) {
    m = get.adjacency(g)
    m = as.matrix(m)
    g2 = new("graphAM", m, edgemode = "directed")
    node = V(g)$name
    node.size = V(g)$size

    if(max(node.size) == min(node.size)) {
        node.size = rep(0.5, length(node.size)) 
    } else {
        node.size = (node.size - min(node.size))/(max(node.size) - min(node.size))+0.2
    }
    names(node.size) = node
    node.color = V(g)$color
    names(node.color) = node
    node.label.size = rep(0.7, length(node))
    names(node.label.size) = node
    node.label = V(g)$label
    names(node.label) = node
    node.shape = rep("circle", length(node))
    node.shape[node.color == "#CCCCCC"] = "circle"
    
    names(node.shape) = node
    
    fix = rep(FALSE, length(node))
    names(fix) = fix
    
    nAttr = list(width = node.size, height = node.size, label = node.label, shape = node.shape)
    x = layoutGraph(g2, nodeAttrs = nAttr)
    nodeRenderInfo(x) = list(fill = node.color, cex = node.label.size)
    if(draw) {
        renderGraph(x)
    }
    return(invisible(x))
}
