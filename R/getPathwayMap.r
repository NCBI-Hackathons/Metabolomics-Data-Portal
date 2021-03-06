#' Get Pathway Map: Generate pathway map with patient data superimposed.
#'
#' @param input - List of parameters passed from R shiny from UI. Important parameters will be input$ptIDs, input$pathwayMapId, and
#'                input$scalingFactor.
#' @param zscore.data - A named vector of metabolites with corresponding z-scores.
#'
#' @return An SVG image of the selected pathway map in input$pathwayMapId, with the mean profile of input$ptIDs
#'         superimposed.
#' @import igraph
#' @import grid
#' @export getPathwayMap
#' @examples
#' data(Miller2015_Heparin)
#' # Input is supplied by R shiny app, but you can hard code parameters as a list object, too, to test functionality.
#' input = list()
#' input$ptIDs = colnames(Miller2015_Heparin)[4]
#' input$diagClass = "paa"
#' input$pathwayMapId = "allPathways"
#' input$scalingFactor = 1
#' svg_img = getPathwayMap(input, Miller2015_Heparin)
#' # SVG file will be saved at following path
#' svg_img$src
#' # Colorbar can be visualized this way
#' plot(svg_img$colorbar)
getPathwayMap = function(input, zscore.data) {
  Pathway.Name = input$pathwayMapId
  if (length(input$ptIDs)==0) {
    return(list(pmap = list(src="", contentType = 'image/svg+xml'), colorbar = NULL))
  } else {
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
    template.ig = getPathwayIgraph(input, Pathway.Name)
    
    node.labels = V(template.ig)$label
    node.types = V(template.ig)$shape
    pmap.path = "./extdata"
    load(sprintf("%s/complexNodes.RData", pmap.path))
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
    V(template.ig)$label = tolower(V(template.ig)$label)
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




