
# == title
# Apply CePa algorithm on a list of pathways under multiple centralities
#
# == param
# -dif differential gene list
# -bk background gene list. If background gene list are not specified, use whole human genes
# -mat expression matrix in which rows are genes and columns are samples
# -label a `sampleLabel` object identify the design of the microarray experiment
# -pc a ``pathway.catalogue`` object storing information of pathways
# -cen centrality measuments, it can ce a string, or a function
# -cen.name centrality measurement names. By default it is parsed from ``cen`` argument
# -nlevel node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs".
#                Also self-defined functions are allowed, see `cepa.univariate.all` for detail.
# -plevel pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank".
#                Also, self-defined functions are allowed, see `cepa.univariate.all` for detail.
# -iter number of simulations
#
# == details
# All the calculation can be achieved by this function. The function is wrapper of 
# both ORA extension and GSA extension. It chooses corresponding procedure according 
# to the arguments specified. If the arguments contain gene lists, then the calculation 
# is sent to functions doing ORA extension. While if the arguments contain an expression 
# matrix and a phenotype label, the GSA extension is evoked. 
#
# The function is a wrapper of `cepa.ora.all` and `cepa.univariate.all`.
#
# This is the core function of the package. User can refer to the vignette to find
# how to use it (``vignette("CePa")``).
#
# If ``dif``, ``bk``, ``pc``, ``cen``, ``cen.name`` and ``iter``
# are specified, the arguments are passed to ``cepa.ora.all``. The centrality-extension 
# of over-representation analysis (ORA) will be applied on the list of differential genes.
#
# If ``mat``, ``label``, ``pc``, ``cen``, ``cen.name``, ``nlevel``,
# ``plevel`` and ``iter`` are specified, the arguments are passed to ``cepa.univariate.all``.
# The centrality-extension of gene-set analysis (GSA) will be applied on the whole gene expressions.
#
# There is a parallel version of the function: `cepa.all.parallel`.
#
# == value
# A `cepa.all` class object
#
# == reference
# Gu Z, Liu J, Cao K, Zhang J, Wang J. Centrality-based pathway enrichment: a systematic 
# approach for finding significant pathways dominated by key genes. BMC Syst Biol. 2012 Jun 6;6(1):56.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `cepa`, `cepa.ora.all`, `cepa.univariate.all`, `cepa.all.parallel`
#
# == examples
# \dontrun{
#
# data(PID.db)
#
# # ORA extension
# data(gene.list)
# # will spend about 20 min
# res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53_symbol.gct")
# label = read.cls("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", treatment="MUT", control="WT")
# # will spend about 45 min
# res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
# }
cepa.all = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000 ) {
    
    # for those who are lazy to specify argument names
    # if the first argument is a vector, then it is ora method
    if(is.vector(dif)) {
        # do nothing
        if(is.null(bk)) {
            dir = system.file(package = "CePa")
            bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        }
    } else if(is.matrix(dif)) {
        mat = dif
        dif = NULL
    } else if(is.data.frame(dif)) {
        mat = as.matrix(dif)
        dif = NULL
    }
    
    if(! is.null(dif)) {     # if dif is specified
        res = cepa.ora.all(dif = dif, bk = bk, pc = pc, cen = cen, cen.name = cen.name, iter = iter)
    } else {
        res = cepa.univariate.all(mat = mat, label = label, pc = pc, cen = cen, cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter)
    }
    
    return(res)
}


# == title
# Apply CePa algorithm on a single pathway
#
# == param
# -dif differential gene list
# -bk background gene list. If background gene list are not specified, use whole human genes
# -mat expression matrix in which rows are genes and columns are samples
# -label a `sampleLabel` object identify the design of the microarray experiment
# -pc a ``pathway.catalogue`` object storing information of pathways
# -pathway an `igraph::igraphtest` object or edge list
# -id identify which pathway should be analysis in the pathway catalogue
# -cen centrality measuments, it can ce a string, or function has been quote
# -cen.name centrality measurement names. This argument should be set if the ``cen`` is a function.
# -nlevel node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs".
#                 Also self-defined functions are allowed, see `cepa.univariate` for detail.
# -plevel pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank".
#                 Also, self-defined functions are allowed, see `cepa.univariate` for detail.
# -iter number of simulations
#
# == details
# The function is a wrapper of `cepa.ora` and `cepa.univariate`.
# Selection of which function depends on the arguments specified.
#
# If ``dif``, ``bk``, ``pc``, ``pathway``, ``id``, ``cen``, ``cen.name`` and ``iter``
# are specified, the arguments are passed to `cepa.ora`. The centrality-extension 
# of over-representation analysis (ORA) will be applied on the list of differential genes.
#
# If ``mat``, ``label``, ``pc``, ``pathway``, ``id``, ``cen``, ``cen.name``, ``nlevel``,
# ``plevel`` and ``iter`` are specified, the arguments are passed to `cepa.univariate`.
# The centrality-extension of gene-set analysis (GSA) will be applied on the whole gene expressions.
#
# This function is always called by `cepa.all`. But you can still use it
# if you want to analysis a single pathway under a specific centrality.
#
# == value
# A `cepa` class object
#
# == seealso
# `cepa.all`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
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
#
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("P53_symbol.gct")
# label = read.cls("P53.cls", treatment="MUT", control="WT")
# # will take about 45 min
# res.gsa = cepa(mat = eset, label = label, pc = PID.db$NCI, id = 2)
# }
cepa = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, pathway = NULL, 
    id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen)) 
               else if(mode(cen) == "name") deparse(cen) else cen,
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000) {
    
    # if the first argument is a vector, then it is ora method
    if(is.vector(dif)) {
        if(is.null(bk)) {
            dir = system.file(package = "CePa")
            bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        }
    } else if(is.matrix(dif)) {
        mat = dif
        dif = NULL
    } else if(is.data.frame(dif)) {
        mat = as.matrix(dif)
        dif = NULL
    }
    
    if(! is.null(dif)) {     # if dif is specified
        cat("  Applying Over-representative analysis.\n")
        res = cepa.ora(dif = dif, bk = bk, pc = pc, pathway = pathway, id = id, cen = cen, cen.name = cen.name, iter = iter)
    } else {
        cat("  Applying Gene-set analysis (univariate procedure).\n")
        res = cepa.univariate(mat = mat, label = label, pathway = pathway, pc = pc, id = id, cen = cen,
                   cen.name = cen.name, iter = iter, nlevel = nlevel, plevel = plevel)
    }
    
    return(res)
    
}

is.ora = function(x) {
    return(class(x) == "cepa.all" && x[[1]][[1]]$framework == "ora")
}

is.gsa = function(x) {
    return(class(x) == "cepa.all" && x[[1]][[1]]$framework == "gsa.univariate")
}



divide = function(x, k) {
    if(length(x) ==1 && is.numeric(x)) {
        x = 1:x
    }
    if(length(x) < k) {
        k = length(x)
    }
    w = floor(length(x)/k)
    q = length(x) - k*w
    d = matrix(0, nrow=k, ncol=2)
    n = 1
    for(i in 1:k) {
        d[i, 1] = n
        d[i, 2] = n+w-1+ifelse(q>0, 1, 0)
        n = d[i,2]+1
        q = ifelse(q > 0, q-1, 0)
    }
    d[k,2] = length(x)
    return(d)
}

combine.cepa.all = function(res) {        
    obj = list()
    for(i in 1:length(res)) {
         obj = c(obj, res[[i]])
    }
    class(obj) = "cepa.all"
    return(obj)
}

# == title
# use CePa package through parallel computing
#
# == param
# -dif differential gene list
# -bk background gene list. If background gene list are not specified, use whole human genes
# -mat expression matrix in which rows are genes and columns are samples
# -label a `sampleLabel` object identify the design of the microarray experiment
# -pc a ``pathway.catalogue`` object storing information of pathways
# -cen centrality measuments, it can ce a string, or a function
# -cen.name centrality measurement names. By default it is parsed from ``cen`` argument
# -nlevel node level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs".
#             Also self-defined functions are allowed, see `cepa.univariate.all` for detail.
# -plevel pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank".
#             Also, self-defined functions are allowed, see `cepa.univariate.all` for detail.
# -iter number of simulations
# -ncores number of cores for parallel computing
#
# == details
# The function divides the pathway list into several parts and each part is sent to a core for 
# parallel computing.
#
# The package for parallel computing is \code{snow}.
#
# Note: there may be warnings saying connections not closed. In fact I have closed
# connections after the parallel computing is done. I don't know why this
# happens. Maybe you breaked the computing ahead manually. However it does not matter 
# unless you have obsessive compulsive disorder.
#
# == value
# A `cepa.all` class object
#
# == reference
# Gu Z, Liu J, Cao K, Zhang J, Wang J. Centrality-based pathway enrichment: a systematic 
# approach for finding significant pathways dominated by key genes. BMC Syst Biol. 2012 Jun 6;6(1):56.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# cepa.all
#
# == example
# \dontrun{
# data(PID.db)
# # ORA extension
# data(gene.list)
# res.ora = cepa.all.parallel(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, ncores = 4)
# # GSA extension
# # P53_symbol.gct and P53_cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53_symbol.gct")
# label = read.cls("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", 
#     treatment="MUT", control="WT")
# res.gsa = cepa.all.parallel(mat = eset, label = label, pc = PID.db$NCI, ncores = 4)
# }
cepa.all.parallel = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, 
    pc, cen = default.centralities, 
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2) {
    
    if(length(pc$pathList) < ncores) {
        stop("Number of cores should not be larger than the number of pathways.")
    }
    
    cat("Use snow package to do parallel computing (", ncores, "cores ) ...\n")
    cepa.all.parallel.by.snow(dif = dif, bk = bk, mat = mat, label = label, pc = pc, cen = cen,
        cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter, ncores = ncores)

}

cepa.all.parallel.by.snow = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2) {
    
    
    # to avoid lazy evaluation
    mode(dif)
    mode(bk)
    mode(mat)
    mode(label)
    mode(pc)
    mode(cen)
    mode(cen.name)
    mode(nlevel)
    mode(plevel)
    mode(iter)

    d = divide(1:length(pc$pathList), ncores)
    if(dim(d)[1] < ncores) {
        ncores = dim(d)[1]
        cat("Since your task can only be divided into", ncores, "parts, modify the number of cores to", ncores, "\n")
    }
    
    cl = makeCluster(ncores, type="SOCK")
    ignore = clusterCall(cl, function() {library(CePa); NULL})
    
    res = clusterApply(cl, 1:ncores, function(i) {
                pc = set.pathway.catalogue(pathList = pc$pathList[d[i, 1]:d[i, 2]],
                                           interactionList = pc$interactionList,
                                           mapping = pc$mapping)
                cepa.all(dif = dif, bk = bk, mat = mat, label = label, pc = pc, cen = cen,
                         cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter)
            })
    stopCluster(cl)
    res.p = combine.cepa.all(res)
    return(res.p)
}
