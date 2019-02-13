#' Convert HMDB IDs to KEGG compound IDs.
#'
#' A function that converts HMDB IDs to KEGG compounds.
#' @param hmdb.ids - A character vector of HMDB IDs.
#' @export data.HMDBtoKEGG
#' @examples
#' data.kegg = data.HMDBtoKEGG(data)
data.HMDBtoKEGG = function(hmdb.ids) {
  KEGG.ids = data.frame(HMDB=character(), KEGG=character(), stringsAsFactors=FALSE)
  print("Mapping HMDB IDs to KEGG IDs using BridgeDB.org...")
  for (met in 1:length(hmdb.ids)) {
    id = hmdb.ids[met]
    if (grepl("HMDB", id, perl=TRUE)) {
      print(sprintf("Processing metabolite (%d/%d) by HMDB ID for %s", met, length(hmdb.ids), id))
      response = getURL(sprintf("http://webservice.bridgedb.org/Human/xrefs/Ch/%s", id))
      tmp = as.matrix(unlist(lapply(strsplit(response, split="\n"), rbind)))
      tmp = t(apply(tmp, 1, function(i) cbind(unlist(strsplit(i, split="\t")))))
      if (response!="" && ("KEGG Compound" %in% tmp[,2])) {
        KEGG.ids[met, "HMDB"] = id
        KEGG.ids[met, "KEGG"] = paste(tmp[which(tmp[,2]=="KEGG Compound"),1], collapse=";")
      } else {
        KEGG.ids[met, "HMDB"] = id
        KEGG.ids[met, "KEGG"] = ""
      }
    }
  }

  return(KEGG.ids)
}
