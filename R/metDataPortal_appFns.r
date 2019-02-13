
getData = function(input) {
  print("called getData()...")
  data = .GlobalEnv$all_data
  pts = as.character(unlist(sapply(input$showThese, function(i) cohorts[[i]])))
  ind = which(colnames(data) %in% pts)
  data = data[,ind]
  data = apply(data, c(1,2), function(i) round(i, 2))
  res = cbind(rownames(data), data)
  colnames(res) = c("Metabolite", colnames(data))
  res = as.matrix(res)
  return(res)
}

getPatientReport = function(input, all_data) {
  print(input$diagClass)
  print(input$ptIDs)
  
  # MetaboliteName Zscore
  all_data = data.matrix(all_data)
  tmp.zscore = rownames(all_data)
  if (length(input$ptIDs)>1) {
    zscore.data = apply(all_data[ , which(colnames(all_data) %in% input$ptIDs)], 1, function(i) mean(na.omit(i)))
  } else {
    zscore.data = all_data[ , which(colnames(all_data)==input$ptIDs)]
  }
  names(zscore.data) = tmp.zscore
  print(head(zscore.data))
  
  ind = which(is.na(zscore.data))
  if (length(ind)>0) {
    zscore.data = zscore.data[-ind]
  }

  data = data.frame(Metabolite=character(), Zscore=numeric(), stringsAsFactors = FALSE)
  for (row in 1:length(zscore.data)) {
    data[row, "Metabolite"] = names(zscore.data)[row]
    data[row, "Zscore"] = round(zscore.data[names(zscore.data)[row]], 2)
  }

  # Remove mets that were NA in zscore 
  ind0 = which(is.na(data[,"Zscore"]))
  if (length(ind0)>0) {
    data = data[-ind0,]
  }
  print(dim(data))
  
  # Order by abs(Zscore)
  class(data[,"Zscore"]) = "numeric"
  data = data[order(abs(data[,"Zscore"]), decreasing = TRUE), ]
  names(data) = c("Metabolite", "Z-score")

  return(list(patientReport=data))
}

getPathwayMap = function(input, zscore.data) {
  #' Generate pathway map with patient data superimposed.
  #' @param Pathway.Name - The name of the pathway map you want to plot patient data on.
  #' @param PatientID - An identifier string associated with the patient.
  #' @param patient.zscore - A named vector of metabolites with corresponding z-scores.
  #' @param scalingFactor - Integer associated with increase in node size.
  #' @param outputFilePath - The directory in which you want to store image files.

  if (length(input$ptIDs)==0) {
    return(list(pmap = list(src="", contentType = 'image/svg+xml'), colorbar = NULL))
  } else {
    Pathway.Name = gsub(" ", "-", input$pathwayMapId)
    PatientID = input$ptIDs
    scalingFactor = input$scalingFactor
    tmp = rownames(zscore.data)
    if (length(input$ptIDs)>1) {
      patient.zscore = zscore.data[,which(colnames(zscore.data) %in% input$ptIDs)]
      patient.zscore = apply(patient.zscore, 1, function(i) mean(na.omit(i)))
    } else {
      patient.zscore = zscore.data[,which(colnames(zscore.data)==input$ptIDs)]
    }
    names(patient.zscore) = tmp
    #print(patient.zscore)

    pmap.path = "../inst/extdata"
    load(sprintf("%s/complexNodes.RData", pmap.path))
    if (Pathway.Name=="All") {
      load(sprintf("%s/RData/allPathways2.RData", pmap.path))
      V(ig)$label[which(V(ig)$label %in% c("DSGEGDFXAEGGGVR", "Dsgegdfxaegggvr"))] = ""
      scalingFactor=1
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
    nms = node.labels[which(node.labels %in% names(patient.zscore))]
    patient.zscore = patient.zscore[which(names(patient.zscore) %in% nms)]
    granularity = 2
    blues = colorRampPalette(c("blue", "white"))(granularity*ceiling(abs(min(na.omit(patient.zscore))))+1)
    reds = colorRampPalette(c("white", "red"))(granularity*ceiling(max(abs(na.omit(patient.zscore)))))
    redblue = c(blues, reds[2:length(reds)])
    for (i in 1:length(node.labels)) {
      if ((node.labels[i] %in% nms) && (!is.na(patient.zscore[node.labels[i]]))) {
        if (!is.na(1 + abs(patient.zscore[node.labels[i]]))) {
          V(template.ig)$size[i] = scalingFactor*ceiling(abs(patient.zscore[node.labels[i]]))
          V(template.ig)$color[i] = redblue[1+granularity*(ceiling(patient.zscore[node.labels[i]])-ceiling(min(na.omit(patient.zscore))))]
        } else {
          V(template.ig)$size[i] = 1
          V(template.ig)$color[i] = "#D3D3D3"
        }
      } else {
        V(template.ig)$size[i] = 1
        V(template.ig)$color[i] = "#D3D3D3"
      }
    }

    mapped=0
    complexNodes = complexNodes[which(names(complexNodes) %in% node.labels)]
    # Next do the complex nodes
    if (length(which(names(complexNodes) %in% node.labels))>0) {
      for (n in 1:length(complexNodes)) {
        metsInComplex = complexNodes[[n]]
        mapped.mets = metsInComplex[which(metsInComplex %in% names(patient.zscore))]
        if (length(na.omit(patient.zscore[mapped.mets]))>0) {
          print(sprintf("%s: %f", names(complexNodes)[n], length(mapped.mets)/length(metsInComplex)))
          mapped = mapped + length(mapped.mets)
          nodeSize = ceiling(max(abs(na.omit(patient.zscore[mapped.mets]))))
          if (is.na(nodeSize)) {
            V(template.ig)$size[which(node.labels==names(complexNodes[n]))] = 1
            V(template.ig)$size2[which(node.labels==names(complexNodes[n]))] = 1
            V(template.ig)$color[which(node.labels==names(complexNodes[n]))] = "#FFFFFF"
          } else {
            V(template.ig)$size[which(node.labels==names(complexNodes[n]))] = scalingFactor*nodeSize
            V(template.ig)$size2[which(node.labels==names(complexNodes[n]))] = scalingFactor*nodeSize
            V(template.ig)$color[which(node.labels==names(complexNodes[n]))] = redblue[1+granularity*(nodeSize-ceiling(min(na.omit(patient.zscore))))]
          }
        } else {
          V(template.ig)$size[which(node.labels==names(complexNodes[n]))] = 1
          V(template.ig)$size2[which(node.labels==names(complexNodes[n]))] = 1
          V(template.ig)$color[which(node.labels==names(complexNodes[n]))] = "#FFFFFF"
        }
      }
    }
    V(template.ig)$size[which(node.types=="Class")]
    V(template.ig)$label = capitalize(tolower(V(template.ig)$label))
    wrap_strings = function(vector_of_strings,width){
      as.character(sapply(vector_of_strings, FUN=function(x){
        paste(strwrap(x, width=width), collapse="\n")
      }))
    }
    V(template.ig)$label = wrap_strings(V(template.ig)$label, 15)
    V(template.ig)$label.cex = 0.75
    template.ig = delete.vertices(template.ig, v=grep(unlist(strsplit(Pathway.Name, split="-"))[1], V(template.ig)$label))

    svg_filename = sprintf("%s/pmap-%s_%s.svg", getwd(), Pathway.Name, input$diagClass)
    svg(filename = svg_filename, width=10, height=10)
    par(mar=c(1,0.2,1,1))
    plot.igraph(template.ig, layout=cbind(V(template.ig)$x, V(template.ig)$y), edge.arrow.size = 0.01, edge.width = 1,
                vertex.frame.color=V(template.ig)$color) #main = gsub("-", " ", Pathway.Name)
    legend('bottom',legend=1:max(ceiling(V(template.ig)$size/scalingFactor)),
           pt.cex=seq(1, ceiling(max(V(template.ig)$size)), scalingFactor),
           col='black',pch=21, pt.bg='white', cex=1, horiz=TRUE)
    dev.off()

    # Get colorbar
    z = seq(floor(min(na.omit(patient.zscore))), ceiling(max(na.omit(patient.zscore))), 1/granularity)
    df = data.frame(Zscores = z[1:length(redblue)],
                    Colors = redblue)
    if (length(which(apply(df, 1, function(i) any(is.na(i)))))>0) {
      df = df[-which(apply(df, 1, function(i) any(is.na(i)))),]
    }
    cb = ggplot(df, aes(x=1:nrow(df), y=Zscores, colour=Zscores)) + geom_point() + #ggtitle(input$diagClass) +
      scale_colour_gradient2(guide = "colourbar", low = "blue", mid="white", high="red") +
      guides(colour = guide_colourbar(draw.llim = min(df$Zscores), draw.ulim = max(df$Zscores),
                                      direction="horizontal", title.position = "top", barwidth = 10, barheight = 2, reverse = FALSE))
    g_legend=function(a.gplot){
      tmp = ggplot_gtable(ggplot_build(a.gplot))
      leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend = tmp$grobs[[leg]]
      return(legend)}
    leg = g_legend(cb);
    return(list(pmap = list(src=svg_filename, contentType = 'image/svg+xml'), colorbar = leg))
  }
}





