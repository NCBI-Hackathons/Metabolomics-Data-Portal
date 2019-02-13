#' Generate pathway map with patient perturbation data superimposed.
#'
#' @param Pathway.Name - The name of the pathway map you want to plot patient data on.
#' @param PatientID - An identifier string associated with the patient.
#' @param patient.zscore - A named vector of metabolites with corresponding z-scores.
#' @param scalingFactor - Integer associated with increase in node size.
#' @param outputFilePath - The directory in which you want to store image files.
#' @export plot.pathwayMap
#' @examples
#' plot.pathwayMap(simMat, path)
plot.pathwayMap = function(Pathway.Name, PatientID, patient.zscore, scalingFactor, outputFilePath) {
  gmlPath = "../inst/extdata"
  load(sprintf("%s/complexNodes.RData", gmlPath))
  if (Pathway.Name=="allPathways") {
    load(sprintf("%s/RData/%s2.RData", gmlPath, Pathway.Name))
  } else {
    load(sprintf("%s/RData/%s.RData", gmlPath, Pathway.Name))
  }
  if (Pathway.Name=="allPathways") {
    V(template.g)$label[which(V(template.g)$label=="Dsgegdfxaegggvr")] = ""
    scalingFactor=1
  } else {
    template.g = ig
  }

  patient.zscore[which(is.na(patient.zscore))] = 0
  patient.zscore = patient.zscore[-which(abs(patient.zscore)<2)]
  if (length(which(names(patient.zscore)=="3-ureidopropionate"))>0) {
    names(patient.zscore)[which(names(patient.zscore)=="3-ureidopropionate")] = "ureidopropionate"
  }

  nodeDisplayNames= read.table(sprintf("%s/%s/DisplayName-%s.txt", gmlPath, Pathway.Name, Pathway.Name), header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = as.numeric(tmp.nms)
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeDisplayNames = as.character(tmp)
  names(nodeDisplayNames) = tmp.nms
  nodeDisplayNames = gsub("\\+", " ", nodeDisplayNames)
  # Load id to node types mappings
  nodeType = read.table(sprintf("%s/%s/CompoundType-%s.txt", gmlPath, Pathway.Name, Pathway.Name), header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = as.numeric(tmp.nms)
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeType = as.character(tmp)
  names(nodeType) = tmp.nms
  nodeType = nodeType[which(names(nodeType) %in% names(nodeDisplayNames))]

  node.labels = vector("character", length = length(V(template.g)$name))
  node.types = vector("character", length = length(V(template.g)$name))
  for (n in 1:length(V(template.g)$name)) {
    node.labels[n] = tolower(as.character(nodeDisplayNames[V(template.g)$name[n]]))
    node.types[n] = as.character(nodeType[V(template.g)$name[n]])
  }
  node.labels = as.character(sapply(node.labels, URLdecode))
  rm(ind, ind2, n, tmp, tmp.nms, nodeDisplayNames, nodeType)

  metabolon_to_data = read.csv(sprintf("%s/metabolon_to_data.txt", gmlPath), sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
  metabolon_to_data = metabolon_to_data[,c(1,2,3,4)]
  metabolon_to_data = apply(metabolon_to_data, 2, tolower)
  # Relabel nodes that have different names in dataset
  for (i in 1:nrow(metabolon_to_data)) {
    if (metabolon_to_data[i, "Data_Label"]!="") {
      lbl = metabolon_to_data[,"PathwayMap_Label"][i]
      ind = which(node.labels==lbl)
      if (length(ind)>0) {
        node.labels[ind] = metabolon_to_data[i, "Data_Label"]
      }
    }
  }

  nms = node.labels[which(node.labels %in% names(patient.zscore))]
  minZscore = ceiling(min(patient.zscore[which(names(patient.zscore) %in% nms)]))-1
  maxZscore = ceiling(max(patient.zscore[which(names(patient.zscore) %in% nms)]))
  blues = colorRampPalette(c("blue", "white"))(abs(minZscore)+1)
  reds = colorRampPalette(c("white", "red"))(abs(maxZscore)+1)
  redblue = c(blues[1:(length(blues)-1)], reds[2:length(reds)])
  mapped = 1
  for (i in 1:length(node.labels)) {
    if (node.labels[i] %in% nms) {
      mapped = mapped + 1
      V(template.g)$size[i] = 1+abs(patient.zscore[which(names(patient.zscore)==node.labels[i])])
      V(template.g)$color[i] = redblue[abs(minZscore)+ceiling(patient.zscore[which(names(patient.zscore)==node.labels[i])])]
    } else {
      V(template.g)$size[i] = 1
      V(template.g)$color[i] = "#000000"
    }
  }

  names(complexNodes) = tolower(names(complexNodes))
  mapped=0
  complexNodes = complexNodes[which(names(complexNodes) %in% node.labels)]
  # Next do the complex nodes
  if (length(which(names(complexNodes) %in% node.labels))>0) {
    for (n in 1:length(complexNodes)) {
      metsInComplex = as.character(sapply(complexNodes[[n]], tolower))
      metsInComplex = gsub("\\*", "", metsInComplex)
      mapped.mets = metsInComplex[which(metsInComplex %in% names(patient.zscore))]
      if (length(mapped.mets)>0) {
        print(sprintf("%s: %f", names(complexNodes)[n], length(mapped.mets)/length(metsInComplex)))
        mapped = mapped + length(mapped.mets)
        nodeSize = max(abs(na.omit(patient.zscore[which(names(patient.zscore) %in% mapped.mets)])))
        if (is.na(nodeSize)) {
          V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 1
          V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 1
          V(template.g)$color[which(node.labels==names(complexNodes[n]))] = "#000000"
        } else {
          V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 1 + nodeSize
          V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 1 + nodeSize
          V(template.g)$color[which(node.labels==names(complexNodes[n]))] = redblue[abs(minZscore)+ceiling(nodeSize)]
        }
      } else {
        V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 1
        V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 1
        V(template.g)$color[which(node.labels==names(complexNodes[n]))] = "#000000"
      }
    }
  }

  V(template.g)$size[which(node.types=="Class")]
  V(template.g)$label = capitalize(tolower(V(template.g)$label))
  wrap_strings = function(vector_of_strings,width){
    as.character(sapply(vector_of_strings, FUN=function(x){
      paste(strwrap(x, width=width), collapse="\n")
    }))
  }
  V(template.g)$label = wrap_strings(V(template.g)$label, 15)
  V(template.g)$label.cex = 0.75
  template.g = delete.vertices(template.g, v=grep(unlist(strsplit(Pathway.Name, split="-"))[1], V(template.g)$label))
  V(template.g)$color[which(V(template.g)$shape=="rectangle")] = rep("#32CD32", length(which(V(template.g)$shape=="rectangle")))

  #png(sprintf("%s/%s-%s.png", outputFilePath, Pathway.Name, PatientID), 1000, 1000)
  #plot.igraph(template.g, layout=cbind(V(template.g)$x, V(template.g)$y), edge.arrow.size = 0.01, edge.width = 1,
  #            vertex.frame.color=V(template.g)$color, main = gsub("-", " ", Pathway.Name))
  #legend('bottom',legend=1:max(ceiling(V(template.g)$size/scalingFactor)),
  #       pt.cex=seq(1, ceiling(max(V(template.g)$size)), scalingFactor),
  #       col='black',pch=21, pt.bg='white', cex=2, horiz=TRUE)
  #dev.off()
  V(template.g)$label.cex = 0.25
  V(template.g)$label = rep("", length(V(template.g)$name))
  svg(sprintf("%s/%s-%s.svg", outputFilePath, Pathway.Name, PatientID), width=15, height=15)
  plot.igraph(template.g, layout=cbind(V(template.g)$x, V(template.g)$y), edge.arrow.size = 0.01, edge.width = 1,
              vertex.frame.color=V(template.g)$color, main = gsub("-", " ", Pathway.Name))
  legend('bottom',legend=1:max(ceiling(V(template.g)$size/scalingFactor)),
         pt.cex=seq(1, ceiling(max(V(template.g)$size)), scalingFactor),
         col='black',pch=21, pt.bg='white', cex=2, horiz=TRUE)
  dev.off()
}


