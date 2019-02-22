#' getPathwayIgraph
#'
#' @param input - A list object of parameters (esp. from R shiny app). Required parameters are ptIDs, diagClass and pathwayMapId.
#' @param Pathway.Name - The name of the pathway map for which you want the topological information.
#' @return template.ig - Igraph object of selected pathway map.
#' @export getPathwayIgraph
#' @import igraph
#' @examples
#'  data(Miller2015_Heparin)
#' # Input is supplied by R shiny app, but you can hard code parameters as a list object, too, to test functionality.
#' input = list()
#' input$ptIDs = colnames(Miller2015_Heparin)[4]
#' input$diagClass = "paa"
#' input$pathwayMapId = "All"
#' ig = getPathwayIgraph(input, Miller2015_Heparin)
#' # Returns a blank template for selected pathway.
#' plot.igraph(ig, edge.arrow.size = 0.01)
getPathwayIgraph = function(input, Pathway.Name) {
  Pathway.Name = gsub(" ", "-", input$pathwayMapId)
  pmap.path = "./inst/extdata"
  if (Pathway.Name=="All") {
    load(sprintf("%s/RData/allPathways2.RData", pmap.path))
    V(ig)$label[which(V(ig)$label %in% c("DSGEGDFXAEGGGVR", "Dsgegdfxaegggvr"))] = ""
    Pathway.Name = "allPathways"
  } else {
    load(sprintf("%s/RData/%s.RData", pmap.path, Pathway.Name))
  }
  template.ig = ig
  
  # Load id to display label mappings
  nodeDisplayNames= read.table(sprintf("%s/%s/DisplayName-%s.txt", pmap.path, Pathway.Name, Pathway.Name),
                               header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = suppressWarnings(as.numeric(tmp.nms))
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeDisplayNames = as.character(tmp)
  names(nodeDisplayNames) = tmp.nms
  nodeDisplayNames = gsub("\\+", " ", nodeDisplayNames)
  # Load id to node types mappings
  nodeType = read.table(sprintf("%s/%s/CompoundType-%s.txt", pmap.path, Pathway.Name, Pathway.Name),
                        header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = suppressWarnings(as.numeric(tmp.nms))
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeType = as.character(tmp)
  names(nodeType) = tmp.nms
  nodeType = nodeType[which(names(nodeType) %in% names(nodeDisplayNames))]
  
  node.labels = vector("character", length = length(V(template.ig)$name))
  node.types = vector("character", length = length(V(template.ig)$name))
  for (n in 1:length(V(template.ig)$name)) {
    node.labels[n] = URLdecode(as.character(nodeDisplayNames[V(template.ig)$name[n]]))
    node.types[n] = as.character(nodeType[V(template.ig)$name[n]])
  }
  
  V(template.ig)$label = node.labels
  V(template.ig)$shape = node.types
  V(template.ig)$shape[grep("Enzyme", V(template.ig)$shape)] = "rectangle"
  V(template.ig)$shape[grep("Metabolite", V(template.ig)$shape)] = "circle"
  template.ig = delete.vertices(template.ig, v=V(template.ig)$name[-which(V(template.ig)$shape %in% c("rectangle", "circle"))])
  
  return(template.ig)
}
