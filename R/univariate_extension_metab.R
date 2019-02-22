#' Align z-scores
#'
#' @param ig igraph object containing pathway
#' @param z_vec Named vector of metabolites z-scores 
#'
#' @return Subset of named z-scores, aligned to i-graph indices
#' 
align_z_ig <- function(ig, z_vec){
  #align
  z_tm <- data.table(z = z_vec, lab = names(z_vec), key = "lab")[, ind_z := .I]
  iglab <- unique(data.table(lab = V(ig)$label, key = "lab"))[, ind_g := .I]
  
  dt_align <- iglab[z_tm][!is.na(ind_g)]
  
  z_out <- dt_align[, z]
  names(z_out) <- dt_align[, lab]
  z_out
}

#' Return z-scores weighted to graph topology
#'
#' @param ig igraph object containing pathway
#' @param z_vec Named vector of metabolites z-scores 
#' @param thresh Value used for thresholding
#' @param cen centrality paremeter used
#'
#' @return Named vector of weighted z-scores
#' 
weight_z_central <- function(ig, z_vec, thresh = 1, cen = "equal.weight"){
  #get weights
  weight = CePa:::centrality(pathway, cen)
  #fix non-zero weights ? huh
  add = 0
  if(sum(weight == 0) != length(weight)) {
    add = min(weight[weight > 0])/100
  }
  weight = weight + ifelse(sum(weight == weight[1]) == length(weight), 0, add)
  
  #align
  z_tm <- data.table(z = z_vec, lab = names(z_vec), key = "lab")[, ind_z := .I]
  iglab <- unique(data.table(lab = V(ig)$label, weight = weight, key = "lab"))[, ind_g := .I]
  iglab[, maxw := max(weight), lab]
  iglab <- iglab[maxw == weight][, maxw := NULL]
  dt_align <- iglab[z_tm][!is.na(ind_g)]
  
  zw_out <- dt_align[, weight * z]
  names(zw_out) <- dt_align[, lab]
  zw_out
}

#' Implementation of univariate CEPA for clinical metabolites data
#'
#' @param ig The pathway igraph
#' @param z_vec The named vector of metabolites z-scores 
#' @param thresh The z-score threshold
#' @param plevel The function used to summarize weighted z-scores (mean by default)
#' @param cen The centrality parameter (equal.weights by default)
#' @param cen.name
#' @param iter The number of iterations
#'
#' @import data.table
#' @return List describing effect across pathway
#' @export
#'
#' @examples
#' cepa.univariate_fix(ig, z_vec)
cepa.univariate_metab <- function(ig, z_vec, thresh = 1.96,
  plevel = "mean",
  cen = "equal.weight",
  cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen),
  iter = 1000){
  
  #Storage for monte-carlo step
  s.random <- numeric(iter)
  ds.random <- matrix(0, iter, 4)
  
  
  #Compute centrality-based weights 
  z_al <- align_z_ig(ig, z_vec)
  #Restrict to values exceeding threshold
  z_th <- z_al[abs(z_al) > thresh]
  
  #if no value remain, return NA's
  if(length(z_th)){
    wz <- weight_z_central(ig, z_th, 1)
    
    ds = quantile(wz, c(1, 0.75, 0.5, 0))
    names(ds) = c("max", "q75", "median", "min")
    s = plevelFun(wz)
    
    for(i in 1:iter) {
      new_z <- sample(z_al, length(z_th)) #resample z-scores
      names(new_z) <- sample(names(z_al), length(z_th), replace = FALSE) #scramble labels
      
      #calculate weights
      new_wz <- weight_z_central(ig, new_z, 1)
      s.random[i] = plevelFun(new_wz)
      ds.random[i, ] = quantile(new_wz, c(1, 0.75, 0.5, 0))
    }
    p.value = (sum(s.random >= s) + 1) / (iter + 1)  
  }else{
    p.value = NA

  }
  
  res = list("score" = s,                                  # pathway score
    "score.distribution" = ds,                    # distribution of node value in the pathway
    "score.random" = s.random,                    # simulated pathway scores
    "score.distribution.random" = ds.random,      # distribution of node value in the pathway in each simulation
    "p.value" = p.value,                          # p value
    "centrality" = cen,                      # centrality name
    "weight" = wz/z_th,                            # weight for each node
    "node.level.t.value" = z_al ,    # value for each node, exclude the centrality part
    "node.level" = z_al,                    # value for each node, exclude the centrality part
    "node.name" = names(z_al),                      # node names
    "pathway" = pathway,                          # pathway in igraph format
    "framework" = "gsa.univariate")                          
  
  class(res) = "cepa"
  
  return(invisible(res))
}



init_cepa_allresults <- function(pathway.name, cen.name){
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
  pathway.result
}

cepa.univariate.metab.all <- function(z_vec, pmap.path, pathway.name, cen = "betweenness",
  cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
  iter = 1000) {
  
  if(length(cen) < 1) stop("cen argument must be specified.\n")
  for(ce in cen) {
    if(is.function(ce)) stop("Functions cannot be used directly because we need the function name, use quote or substitute.\n")
  }
  
  cat("Calculating pathway scores...\n")
  pathway.result <- init_cepa_allresults(pathway.name, cen.name)
  n.pathway <- length(pathway.name)
  
  for(i in seq_along(pathway.name)) {
    cat("    ", i, "/", n.pathway, ", ", pathway.name[i], "...\n", sep="")
    pathway = getPathwayIgraph(Pathway.Name=pathway.name[i], pmap.path=pmap.path)
    
    for(j in seq_along(cen)) {
      pathway.result[[i]][[j]] = cepa.univariate_metab(z_vec, pathway, 
        thresh, plevel, cen = cen[j], cen.name = cen.name[j], iter = iter)
      
      cat("      - ", cen[j], ": ", round(pathway.result[[i]][[j]]$p.value, 3), "\n", sep = "")
    }
    
  }
}
