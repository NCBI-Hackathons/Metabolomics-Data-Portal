# Loop through all Metabolon Pathways and get a list of unique metabolites.


#' Get All Metabolites In Metabolon's Pathway Knowledgebase
#' @return pathway.names - a character vector of all pathways defined in Metabolon's pathway knowledgebase.
#' @export list_all_pathways
#' @import igraph
#' @examples
#' metabolon_knowledgebase_path = sprintf("%s/extdata", find.package("MetabolomicsDataPortal"))
#' pwys = list_all_pathways(metabolon_knowledgebase_path)
#' print(pwys)
list_all_pathways = function(Pathway.Knowledgebase.Path) {
  ig_files = list.files(sprintf("%s/extdata/RData/", getwd()), pattern = ".RData")
  ig_files = ig_files
  pwys = unlist(sapply(ig_files, function(i) unlist(strsplit(i, split=".RData"))[1]))
  pwys = as.character(pwys)
  pwys = gsub("-", " ", pwys)
  pwys[which(pwys=="allPathways")] = "All Pathways"
  
  return (pwys)
}
