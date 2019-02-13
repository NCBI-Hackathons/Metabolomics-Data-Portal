# == title
# Table of p-values of pathways
#
# == param
# -x a `cepa.all` object
# -adj.method methods in `stats::p.adjust`, available methods are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
# -cutoff cutoff for significance
#
# == details
# Since the p-values for each pathway are calculated for several centralities, the
# whole p-values are represented as a table.
#
# Also it can extract significant pathways only.
#
# == value
# A data matrix where rows are pathways and columns are centralities.
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
# res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
# p.table(res.ora)
# p.table(res.ora, adj.method = "BH")
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# p.table(res.gsa)
# }
p.table = function(x, adj.method = NA, cutoff = ifelse(adj.method == "none", 0.01, 0.05)) {

    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    n.pathway = length(x)
    p.value = matrix(0, nrow=length(x), ncol= length(x[[1]]))
    for(i in 1:length(x)) {
        p.value[i, ] = sapply(x[[i]], function(x) x$p.value)
    }

    rownames(p.value) = names(x)
    colnames(p.value) = names(x[[1]])
    
    if(!is.na(adj.method)) {
        p.value = apply(p.value, 2, p.adjust, adj.method)
        l = apply(p.value, 1, function(x) { sum(x <= cutoff) > 0})
        p.value = p.value[l, , drop=FALSE]
    }
    
    return(p.value)
}

# == title
# get single cepa object from cepa.all object
#
# == param
# -x a `cepa.all` object
# -id index or the name of the pathway
# -cen index or the name of the centrality
#
# == details
# The `cepa.all object contains the result for pathways under several centrality
# measurements. In `cepa.all` object, each pathway under a specific centrality
# is a single `cepa` object. The `get.cepa` function is used to get the `cepa`
# object from the `cepa.all` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `cepa`, `cepa.all`
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
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# gsa = get.cepa(res.gsa, id = 5, cen = 3)
# }
get.cepa = function(x, id = NULL, cen = 1) {
    
    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    if(is.null(id)) {
        stop("id cannot be null")
    }
    
    if(length(cen) > 1) {
        stop("Length of cen must be equal to 1.\n")
    }
    
    if(is.function(cen)) {
        cen = deparse(substitute(cen))
    }
    else if(mode(cen) == "name") {
        cen = deparse(cen)
    }
    
    return(x[[id]][[cen]])
}
    
