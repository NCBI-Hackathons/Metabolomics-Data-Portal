
getMDS = function(input) {
  patientSim = graphs[[input$diagnosis]][["patientSim"]]
  d = graphs[[input$diagnosis]][["diagnoses"]]
  if (input$diagnosis=="zsd") {
    d = c(rep("PEX1", 19), rep("PEX7", 21), rep("negCntl", 40))
  }
  if (input$dim==2) {
    fitSim = cmdscale(patientSim, eig=FALSE, k=2)
    x = round(fitSim[,1], 2)
    y = round(fitSim[,2], 2)
    df = data.frame(x=x, y=y, color=d)
    mds_plot = plot_ly(df, x=~x, y=~y, color=~color, marker = list(size = 20))
  } else {
    fitSim = cmdscale(patientSim, eig=FALSE, k=3)
    x = round(fitSim[,1], 2)
    y = round(fitSim[,2], 2)
    z = round(fitSim[,3], 2)
    df = data.frame(x=x, y=y, z=z, color=d, label=colnames(patientSim))
    mds_plot = plot_ly(df, x=~x, y=~y, z=~z, color=~color, text=~label)
  }
  return(mds_plot)
}

getData = function(input) {
  print("called getData()...")
  if (input$raworZscore == "Raw") {
    data = .GlobalEnv$all_raw_data
  } else if (input$raworZscore == "Normalized") {
    data = .GlobalEnv$all_norm_data
  } else if (input$raworZscore == "Zscored") {
    data = .GlobalEnv$all_data
  }
  pts = as.character(unlist(sapply(input$showThese, function(i) cohorts[[i]])))
  ind = which(colnames(data) %in% pts)
  data = data[,ind]
  data = apply(data, c(1,2), function(i) round(i, 2))
  res = cbind(rownames(data), data)
  colnames(res) = c("Metabolite", colnames(data))
  res = as.matrix(res)
  return(res)
}

getMetList = function(input) {
  # First, get rid of metabolites that have below fil rate
  ref = data.matrix(.GlobalEnv$all_norm_data)
  r.ids = unlist(sapply(ref.ids, function(i) sprintf("X%s", i)))
  ref = ref[,which(colnames(ref) %in% r.ids)]
  ref.fil = apply(ref, 1, function(i) 1-(sum(is.na(i))/length(i)))
  ref = ref[which(ref.fil>0.66),]
  metClass = .GlobalEnv$metClass[which(ref.fil>0.66)]

  if (input$metClass=="Lipid") {
    return(rownames(ref)[which(metClass=="Lipid")])
  } else if (input$metClass=="Unknown") {
    return(rownames(ref)[which(metClass=="Unknown")])
  } else if (input$metClass=="Nucleotide") {
    return(rownames(ref)[which(metClass=="Nucleotide")])
  } else if (input$metClass=="Amino Acid") {
    return(rownames(ref)[which(metClass=="Amino Acid")])
  } else if (input$metClass=="Cofactors and Vitamins") {
    return(rownames(ref)[which(metClass=="Cofactors and Vitamins")])
  } else if (input$metClass=="Xenobiotics") {
    return(rownames(ref)[which(metClass=="Xenobiotics")])
  } else if (input$metClass=="Carbohydrate") {
    return(rownames(ref)[which(metClass=="Carbohydrate")])
  } else if (input$metClass=="Energy") {
    return(rownames(ref)[which(metClass=="Energy")])
  } else if (input$metClass=="Peptide") {
    return(rownames(ref)[which(metClass=="Peptide")])
  }
}

neg = function(x) { return(-x) }

# Get reference population statistics & plots
getRefPop = function(input, norm.data) {
  print("getRefPop() called.")
  r.ids = unlist(sapply(ref.ids, function(i) sprintf("X%s", i)))
  ref = norm.data[, which(colnames(norm.data) %in% r.ids)]
  ref.fil = apply(ref, 1, function(i) 1-(sum(is.na(i))/length(i)))
  ref = ref[which(ref.fil>0.66),]
  print(dim(ref))
  print(input$metSelect)
  print(input$metSelect %in% rownames(ref))

  # Histogram plot of metabolite, with outliers that were removed during z-score calculation highlighted in red
  outlierSamples = which(ref[input$metSelect,] %in% boxplot.stats(ref[input$metSelect,])$out)
  print(sprintf("Length outlier samples = %d", length(outlierSamples)))
  if (length(outlierSamples)>0) {
    df = data.frame(x=as.numeric(ref[input$metSelect, -outlierSamples]))
  } else {
    df = data.frame(x=as.numeric(ref[input$metSelect,]))
  }
  hst = ggplot(data=df, aes(x=x)) + geom_histogram() +
    ggtitle(sprintf("%s (%s)", input$metSelect, input$metClass)) + labs(x="log(NormScaledImputed)", y="Count")
  y = quantile(df$x[!is.na(df$x)], c(0.05, 0.95))
  x = qnorm(c(0.05, 0.95))
  slope = diff(y)/diff(x)
  int = y[1L] - slope * x[1L]
  qq = ggplot(data=df, aes(sample=x)) + stat_qq() + ggtitle(sprintf("Normal QQ-Plot for %s (%s)", input$metSelect, input$metClass)) +
    geom_abline(slope = slope, intercept = int)

  per = list(up=length(which(.GlobalEnv$all_data[input$metSelect,]>2))/nrow(.GlobalEnv$all_data),
             down=length(which(neg(.GlobalEnv$all_data[input$metSelect,]) > 2))/nrow(.GlobalEnv$all_data))
  df = data.frame(Sample=1:ncol(.GlobalEnv$all_data), Zscore=.GlobalEnv$all_data[input$metSelect,])
  rare = ggplot(data=df, aes(x=Sample, y=Zscore)) + geom_point(size=1) +
    ggtitle(sprintf("Percentage with zscore >2 = %.2f.\nPercentage with zscore <-2 = %.2f.", per$up, per$down)) +
    geom_hline(yintercept=2, color="red") + geom_hline(yintercept=-2, color="red")

  x = ref[input$metSelect,-outlierSamples]
  d = qqnorm(as.numeric(x), plot.it = FALSE)
  xx = d$y
  zz = d$x
  t = lm(xx~zz, data=as.data.frame(x=xx, z=zz))
  mn.est = as.numeric(t$coefficients[1])
  sd.est = as.numeric(t$coefficients[2])

  samples = colnames(ref)[which(ref[input$metSelect,] %in% boxplot.stats(ref[input$metSelect,])$out)]
  values = ref[input$metSelect, which(ref[input$metSelect,] %in% boxplot.stats(ref[input$metSelect,])$out)]
  outlierSamples = cbind(samples, round(as.numeric(values),2))
  colnames(outlierSamples) = c("Samples Outliers", "Sample Value")

  return (list(hst=hst, outliers=outlierSamples, qq=qq, ests=list(mean=mn.est, std=sd.est), rare=rare, per=per))
}

getPatientReport = function(input, raw.data, norm.data, zscore.data) {
  raw.data = data.matrix(raw.data)
  norm.data = data.matrix(norm.data)
  zscore.data = data.matrix(zscore.data)

  r.ids = unlist(sapply(ref.ids, function(i) sprintf("X%s", i)))
  ref = norm.data[,which(colnames(norm.data) %in% r.ids)]
  ref.fil = apply(ref, 1, function(i) 1-(sum(is.na(i))/length(i)))

  # MetaboliteName  RawIonIntensity Anchor(CMTRX.5 median value)  Zscore
  tmp = rownames(raw.data)
  tmp.zscore = rownames(zscore.data)
  ind = which(rownames(zscore.data) %in% rownames(norm.data))
  if (length(input$ptIDs)>1) {
    if (length(which(colnames(raw.data) %in% input$ptIDs))>1) {
      raw.data = apply(raw.data[,which(colnames(raw.data) %in% input$ptIDs)], 1, function(i) mean(na.omit(i)))
      norm.data = apply(norm.data[,which(colnames(norm.data) %in% input$ptIDs)], 1, function(i) mean(na.omit(i)))
    } else {
      raw.data = rep(NA, nrow(raw.data))
      norm.data = rep(NA, nrow(norm.data))
    }
    zscore.data = apply(zscore.data[ind, which(colnames(zscore.data) %in% input$ptIDs)], 1, function(i) mean(na.omit(i)))
  } else {
    raw.data = rep(NA, nrow(raw.data))
    norm.data = rep(NA, nrow(norm.data))
    zscore.data = zscore.data[ind, which(colnames(zscore.data)==input$ptIDs)]
  }
  names(raw.data) = tmp
  names(norm.data) = tmp
  names(zscore.data) = tmp.zscore[ind]

  data = data.frame(Metabolite=character(), Raw=numeric(), Anchor=numeric(), Zscore=numeric(), stringsAsFactors = FALSE)
  for (row in 1:length(raw.data)) {
    data[row, "Metabolite"] = names(raw.data)[row]
    data[row, "Raw"] = round(raw.data[row], 2)
    if (length(which(names(norm.data)==names(raw.data)[row]))>0) {
      data[row, "Anchor"] = round(norm.data[which(names(norm.data)==names(raw.data)[row])], 2)
    } else {
      data[row, "Anchor"] = NA
    }
    if (length(which(names(zscore.data)==names(raw.data)[row]))>0) {
      data[row, "Zscore"] = round(zscore.data[which(names(zscore.data)==names(raw.data)[row])], 2)
    } else {
      data[row, "Zscore"] = NA
    }
  }

  # Remove mets that were NA in raw, norm and zscore AND
  # Next, Remove mets that were NA in raw, but not in Anchor. These will be displayed in separate table.
  # Note, these values were imputed and therefore should not be included in patient report, but should
  # be noted that these metabolites were normally found.
  # Last, remove mets that were NA in raw, but not in Zscore.
  ind0 = intersect(intersect(which(is.na(data[,"Raw"])), which(is.na(data[,"Anchor"]))), which(is.na(data[,"Zscore"])))
  ind1 = intersect(which(is.na(data[,"Raw"])), which(!is.na(data[,"Anchor"])))
  #ind2 = intersect(which(is.na(data[,"Raw"])), which(!is.na(data[,"Zscore"])))
  ind_all = data[,"Metabolite"][unique(c(ind0, ind1))] #ind2
  tmp = ref.fil[ind_all]
  report_these = tmp[which(tmp>0.80)]
  # Report these metabolites
  missingMets = data.frame(Metabolite=character(), Reference.FillRate=numeric(), stringsAsFactors = FALSE)
  if (length(report_these)>0) {
    for (i in 1:length(report_these)) {
      met = names(report_these)[i]
      missingMets[i, "Metabolite"] = met
      missingMets[i, "Reference.FillRate"] = ref.fil[which(names(ref.fil)==met)]
    }
    colnames(missingMets) = c("Compound", "Reference Fill Rate")
  } else {
    missingMets = NULL
  }

  if (length(ind_all)>0) {
    data = data[-unique(c(ind0, ind1)),]
  }
  print(dim(data))

  # Order by Fill Rate
  missingMets = missingMets[order(missingMets[,"Reference Fill Rate"], decreasing = TRUE),]

  # Order by abs(Zscore)
  class(data[,"Zscore"]) = "numeric"
  data = data[order(abs(data[,"Zscore"]), decreasing = TRUE), ]
  names(data) = c("Metabolite", "Raw Ion Intensity", "Anchor", "Z-score")

  return(list(patientReport=data, missingMets=missingMets))
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
    patient.zscore = zscore.data[,which(colnames(zscore.data) %in% input$ptIDs)]
    patient.zscore = apply(patient.zscore, 1, function(i) mean(na.omit(i)))
    names(patient.zscore) = tmp
    #print(patient.zscore)

    pmap.path = "../../inst/extdata"
    load(sprintf("%s/complexNodes.RData", pmap.path))
    if (Pathway.Name=="All") {
      load(sprintf("%s/RData/allPathways.RData", pmap.path))
      V(template.ig)$label[which(V(template.ig)$label %in% c("DSGEGDFXAEGGGVR", "Dsgegdfxaegggvr"))] = ""
      scalingFactor=1
      Pathway.Name = "allPathways"
    } else {
      load(sprintf("%s/RData/%s.RData", pmap.path, Pathway.Name))
      template.ig = ig
    }

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





