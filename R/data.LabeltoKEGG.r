#' Convert compound names to KEGG compound IDs.
#'
#' A function that converts HMDB IDs to KEGG compounds.
#' @param compound.names - A character vector of metabolite (compound) names.
#' @export data.LabeltoKEGG
#' @examples
#' data.kegg = data.HMDBtoKEGG(data)
data.LabeltoKEGG = function(compound.names) {
  KEGG.ids = data.frame(Name=character(), KEGG=character(), stringsAsFactors=FALSE)
  for (met in 1:length(compound.names)) {
    id = compound.names[met]
    # Query by metabolite label
    # Remove any part of string that has ()*:
    if (!is.na(id)) {
      print(sprintf("Processing metabolite (%d/%d) by compound name for %s", met, length(compound.names), id))
      response = getURL(sprintf("http://webservice.bridgedb.org/Human/attributeSearch/%s", id))
      tmp = as.matrix(unlist(lapply(strsplit(response, split="\n"), rbind)))
      tmp = t(apply(tmp, 1, function(i) cbind(unlist(strsplit(i, split="\t")))))
      if (response!="" && ("HMDB" %in% tmp[,2])) {
        res = tmp[which(tmp[,2]=="HMDB"),1]
        map = data.HMDBtoKEGG(res)
        KEGG.ids[met, "Name"] = id
        KEGG.ids[met, "KEGG"] = paste(map[which(map[,"KEGG"]!=""),"KEGG"], collapse=";")
      } else {
        KEGG.ids[met, "Name"] = id
        KEGG.ids[met, "KEGG"] = ""
      }
    }
  }
  return(KEGG.ids)
}
