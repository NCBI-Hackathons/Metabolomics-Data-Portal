#' Metabolite set enrichment analysis (MSEA) (using a hypergeometric test) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' @param met.profile - A character vector of a patient's metabolomic profile, including KEGG IDs
#' and the associated z-score or p-value describing the level of the metabolite compared to controls.
#' @param threshold - A cutoff to select metabolites with a zscore > threshold or < -1*threshold.
#' @param type - Either "p-value" or "z-score".
#' @param gene.profile - Default set to NULL, meaning the default enrichment analysis only considers metabolites. However, if you have gene data, too, set
#' this parameter to a character vector of the gene names with found variants in the patient's record. Gene IDs must be converted to Entrez Identifiers.
#' @export data.getMSEA_Metabolon
#' @examples
#'
#' pathway.data = data.getMSEA_Metabolon(met.profile, threhold=3, "z-score", NULL)
data.getMSEA_Metabolon = function(met.profile, threshold=3, type="zscore", gene.profile=NULL) {
  met.profile = met.profile[which(!(is.na(met.profile)))]
  if (type=="zscore") {
    perturbed.mets = met.profile[which(abs(met.profile) > threshold)]
  } else {
    perturbed.mets = met.profile[which(met.profile < 0.05)]
  }
  nms.perturbed.mets = unique(as.character(unlist(sapply(names(perturbed.mets), function(i) strsplit(i, split=";")))))
  
  # The size of the population of total possible metabolites to draw from
  population = names(met.profile)
  paths.hsa = list.dirs(path="../extdata", full.names = FALSE)
  paths.hsa = paths.hsa[-which(paths.hsa %in% c("", "RData", "allPathways"))]
  row = 1
  pathway.data = data.frame(Pathway=character(), FDR=numeric(), Pvalue=numeric(), Hits=integer(), Size=integer(), stringsAsFactors = FALSE)
  for (pathway in 1:length(paths.hsa)) {
    load(sprintf("../extdata/RData/%s.RData", paths.hsa[pathway]))
    
    pathway.compounds = V(ig)$label[which(V(ig)$shape=="circle")]
    pathCompIDs = unique(tolower(pathway.compounds[which(pathway.compounds %in% population)]))
    # q (sample successes), m (population successes), n (population failures), k (sample size)
    sampleSuccesses = length(which(nms.perturbed.mets %in% pathCompIDs))
    populationSuccesses = length(intersect(pathCompIDs, population))
    N = length(population)
    populationFailures=N-populationSuccesses
    numDraws=length(perturbed.mets)
    if (populationSuccesses>0) {
      pathway.data[row, "Pathway"] = paths.hsa[pathway]
      pathway.data[row, "Pvalue"] = phyper(q=sampleSuccesses-1, m=populationSuccesses, n=populationFailures, k=numDraws, lower.tail=FALSE)
      pathway.data[row, "Hits"] = sampleSuccesses
      pathway.data[row, "Size"] = populationSuccesses
      row = row + 1
    }
  }
  pathway.data[,"FDR"] = p.adjust(pathway.data[,"Pvalue"], method="fdr")
  # Sort by FDR, then Pvalue
  pathway.data = pathway.data[order(pathway.data[,"FDR"], pathway.data[,"Pvalue"]),]
  
  # Then, only return rows that have Pvalue < 0.25
  pathway.data = pathway.data[which(pathway.data[,"Pvalue"]<0.25),]
  
  return(pathway.data)
}

