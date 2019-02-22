# Loop through all Metabolon Pathways and get a list of unique metabolites.


#' Get All Metabolites In Metabolon's Pathway Knowledgebase
#' @return mets - a character vector of unique metabolites and enzymes found in at least 1 pathway in Metabolon's pathway knowledgebase.
#' @export get_all_metabolites_metabolon
#' @import igraph
#' @examples
#' mets = get_all_metabolites_metabolon()
#' print(mets)
get_all_metabolites_metabolon = function() {
  ig_files = list.files(sprintf("%s/extdata/RData/", getwd()), pattern = ".RData")
  ig_files = ig_files

  mets = c()
  for (pwy in 1:length(ig_files)) {
    print(pwy)
    load(sprintf("%s/extdata/RData/%s", getwd(), ig_files[pwy]))
    length(V(ig)$label)
    mets = c(mets, V(ig)$label)
  }
  mets = unique(mets)
  
  return (mets)
}
