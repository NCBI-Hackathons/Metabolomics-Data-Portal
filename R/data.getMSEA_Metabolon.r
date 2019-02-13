#' Metabolite set enrichment analysis (MSEA) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' @param met.profile - A character vector of a patient's metabolomic profile, including KEGG IDs
#'                      and the associated z-score or p-value describing the level of the metabolite compared to controls.
#' @param p.corr - A vector of correlations representing the correlation of the metabolite level with a disease phenotype.
#'                     The names of the vector are the metabolite names.
#' @export data.getMSEA_Metabolon
#' @examples
#' data(Miller2015_Heparin)
#' diagnoses = gsub("[[:digit:]]", "", colnames(Miller2015_Heparin))
#' diag.ind = diagnoses
#' diag.ind[which(diag.ind!="Argininemia")] = 0
#' diag.ind[which(diag.ind=="Argininemia")] = 1
#' diag.ind = as.numeric(diag.ind)
#' profile.ind = which(diagnoses=="Argininemia")[1]
#' met.profile = Miller2015_Heparin[,profile.ind]
#' names(met.profile) = rownames(Miller2015_Heparin)
#' 
#' population = names(met.profile)
#' paths.hsa = list.dirs(path="../inst/extdata", full.names = FALSE)
#' paths.hsa = paths.hsa[-which(paths.hsa %in% c("", "RData", "allPathways"))]
#' sink("MetaboliteSetDatabases/Miller2015.gmt")
#' for (p in 1:length(paths.hsa)) {
#'   load(sprintf("../inst/extdata/RData/%s.RData", paths.hsa[p]))
#'   pathway.compounds = V(ig)$label[which(V(ig)$shape=="circle")]
#'   pathCompIDs = unique(tolower(pathway.compounds[which(pathway.compounds %in% population)]))
#'   print(sprintf("%s         %s", paths.hsa[p], paste(pathCompIDs, collapse="    ")), quote=FALSE)
#' }
#' sink()
#' print("test")
#' abs_filename_dataset = "Datasets/Miller2015.gct"
#' abs_filename_classes = "Datasets/Miller2015_arg.cls"
#' pathway.data = data.getMSEA_Metabolon(met.profile, diag.ind, Miller2015_Heparin)
data.getMSEA_Metabolon = function(abs_filename_dataset, abs_filename_classes, pathway_knowledgebase = "Metabolon", output_dir = getwd(), expt_name="msea_results") {
  if (pathway_knowledgebase=="Metabolon") {
    met.db = "Pathway_GMTs/Metabolon.gmt"
  } else if (pathway_knowledgebase=="KEGG") {
    met.db = "Pathway_GMTs/KEGG.gmt"
  } else if (pathway_knowledgebase=="SMPDB") {
    met.db = "Pathway_GMTs/SMPDB.gmt"
  } else if (pathway_knowledgebase=="Reactome") {
    met.db = "Pathway_GMTs/Reactome.gmt"
  } else {
    # WikiPathways
    met.db = "Pathway_GMTs/WikiPathways.gmt"
  }
  res = MSEA(input.ds =  abs_filename_dataset, input.cls = abs_filename_classes, met.db = met.db,
             output.directory = output_dir, doc.string=expt_name, 
             reshuffling.type="sample.labels", nperm=1000, weighted.score.type=1, 
             nom.p.val.threshold=-1, fwer.p.val.threshold=-1, fdr.q.val.threshold=0.25, topmet = 15, 
             adjust.FDR.q.val = F, met.size.threshold.min = 5, met.size.threshold.max = 50, preproc.type = 0, 
             random.seed = 760435, perm.type = 0, fraction = 1.0, replace = F)
  return(res)
}

# M S E A -- Metabolite Set Enrichment Analysis
MSEA.MetaboliteRanking_SingleProfile = function(A, class.labels, metabolite.labels, nperm, permutation.type = 0, sigma.correction = "MetaboliteCluster", fraction=1.0, replace=F) { 
  A = A + 0.00000001
  N = length(A[,1])
  Ns = length(A[1,])
  
  subset.mask = matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels1 = matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels2 = matrix(0, nrow=Ns, ncol=nperm)
  class.labels1 = matrix(0, nrow=Ns, ncol=nperm)
  class.labels2 = matrix(0, nrow=Ns, ncol=nperm)
  
  order.matrix = matrix(0, nrow = N, ncol = nperm)
  obs.order.matrix = matrix(0, nrow = N, ncol = nperm)
  s2n.matrix = matrix(0, nrow = N, ncol = nperm)
  obs.s2n.matrix = matrix(0, nrow = N, ncol = nperm)
  
  obs.metabolite.labels = vector(length = N, mode="character")
  obs.metabolite.descs = vector(length = N, mode="character")
  obs.metabolite.symbols = vector(length = N, mode="character")
  
  M1 = matrix(0, nrow = N, ncol = nperm)
  M2 = matrix(0, nrow = N, ncol = nperm)
  S1 = matrix(0, nrow = N, ncol = nperm)
  S2 = matrix(0, nrow = N, ncol = nperm)
  
  C = split(class.labels, class.labels)
  class1.size = length(C[[1]])
  class2.size = length(C[[2]])
  class1.index = seq(1, class1.size, 1)
  class2.index = seq(class1.size + 1, class1.size + class2.size, 1)
  for (r in 1:nperm) {
    class1.subset = sample(class1.index, size = ceiling(class1.size*fraction), replace = replace)
    class2.subset = sample(class2.index, size = ceiling(class2.size*fraction), replace = replace)
    class1.subset.size = length(class1.subset)
    class2.subset.size = length(class2.subset)
    subset.class1 = rep(0, class1.size)
    for (i in 1:class1.size) {
      if (is.element(class1.index[i], class1.subset)) {
        subset.class1[i] = 1
      }
    }
    subset.class2 = rep(0, class2.size)
    for (i in 1:class2.size) {
      if (is.element(class2.index[i], class2.subset)) {
        subset.class2[i] = 1
      }
    }
    subset.mask[, r] = as.numeric(c(subset.class1, subset.class2))
    fraction.class1 = class1.size/Ns
    fraction.class2 = class2.size/Ns
    if (permutation.type == 0) { # random (unbalanced) permutation
      full.subset = c(class1.subset, class2.subset)
      label1.subset = sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] = rep(0, Ns)
      reshuffled.class.labels2[, r] = rep(0, Ns)
      class.labels1[, r] = rep(0, Ns)
      class.labels2[, r] = rep(0, Ns)
      for (i in 1:Ns) {
        m1 = sum(!is.na(match(label1.subset, i)))
        m2 = sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] = m1
        reshuffled.class.labels2[i, r] = m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] = m2
          class.labels2[i, r] = 0
        } else {
          class.labels1[i, r] = 0
          class.labels2[i, r] = m2
        }
      }
    } else if (permutation.type == 1) { # proportional (balanced) permutation
      class1.label1.subset = sample(class1.subset, size = ceiling(class1.subset.size*fraction.class1))
      class2.label1.subset = sample(class2.subset, size = floor(class2.subset.size*fraction.class1))
      reshuffled.class.labels1[, r] = rep(0, Ns)
      reshuffled.class.labels2[, r] = rep(0, Ns)
      class.labels1[, r] = rep(0, Ns)
      class.labels2[, r] = rep(0, Ns)
      for (i in 1:Ns) {
        if (i <= class1.size) {
          m1 = sum(!is.na(match(class1.label1.subset, i)))
          m2 = sum(!is.na(match(class1.subset, i)))
          reshuffled.class.labels1[i, r] = m1
          reshuffled.class.labels2[i, r] = m2 - m1
          class.labels1[i, r] = m2
          class.labels2[i, r] = 0
        } else {
          m1 = sum(!is.na(match(class2.label1.subset, i)))
          m2 = sum(!is.na(match(class2.subset, i)))
          reshuffled.class.labels1[i, r] = m1
          reshuffled.class.labels2[i, r] = m2 - m1
          class.labels1[i, r] = 0
          class.labels2[i, r] = m2
        }
      }
    }
  }
  # compute S2N for the random permutation matrix
  P = reshuffled.class.labels1 * subset.mask
  n1 = sum(P[,1])         
  M1 = A %*% P
  M1 = M1/n1      
  A2 = A*A        
  S1 = A2 %*% P   
  S1 = S1/n1 - M1*M1    
  if (n1>1) {
    S1 = sqrt(abs((n1/(n1-1)) * S1))
  }
  P = reshuffled.class.labels2 * subset.mask
  n2 = sum(P[,1])           
  M2 = A %*% P           
  M2 = M2/n2          
  A2 = A*A           
  S2 = A2 %*% P      
  S2 = S2/n2 - M2*M2 
  if (n2>1) {
    S2 = sqrt(abs((n2/(n2-1)) * S2))
  }
  if (sigma.correction == "MetaboliteCluster") {  # small sigma "fix" as used in MetaboliteCluster
    S2 = ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 = ifelse(S2 == 0, 0.2, S2)
    S1 = ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 = ifelse(S1 == 0, 0.2, S1)
  }
  M1 = M1 - M2
  S1 = S1 + S2
  s2n.matrix = M1/S1
  for (r in 1:nperm) {order.matrix[, r] = order(s2n.matrix[, r], decreasing=T)}
  # compute S2N for the "observed" permutation matrix
  P = class.labels1 * subset.mask
  n1 = sum(P[,1])  
  if (n1>1) {
    M1 = A %*% P
    M1 = M1/n1      
    A2 = A*A        
    S1 = A2 %*% P   
    S1 = S1/n1 - M1*M1 
    S1 = sqrt(abs((n1/(n1-1)) * S1))   
  } else {
    M1 = A %*% P
    A2 = A*A        
    S1 = A2 %*% P   
    S1 = S1 - M1*M1
    S1 = sqrt(abs(S1))   
  }
  P = class.labels2 * subset.mask
  n2 = sum(P[,1]) 
  if (n2>1) {
    M2 = A %*% P           
    M2 = M2/n2          
    A2 = A*A           
    S2 = A2 %*% P      
    S2 = S2/n2 - M2*M2 
    S2 = sqrt(abs((n2/(n2-1)) * S2))
  } else {
    M2 = A %*% P           
    A2 = A*A           
    S2 = A2 %*% P      
    S2 = S2 - M2*M2 
    S2 = sqrt(abs(S2))
  }
  if (sigma.correction == "MetaboliteCluster") {  # small sigma "fix" as used in MetaboliteCluster
    S2 = ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 = ifelse(S2 == 0, 0.2, S2)
    S1 = ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 = ifelse(S1 == 0, 0.2, S1)
  } 
  M1 = M1 - M2
  S1 = S1 + S2
  obs.s2n.matrix = M1/S1
  for (r in 1:nperm) {obs.order.matrix[,r] = order(obs.s2n.matrix[,r], decreasing=T)            }
  
  return(list(s2n.matrix = s2n.matrix, obs.s2n.matrix = obs.s2n.matrix, order.matrix = order.matrix, obs.order.matrix = obs.order.matrix))
}

MSEA.EnrichmentScore = function(metabolite.list, metabolite.set, weighted.score.type = 1, correl.vector = NULL) {  
  tag.indicator = sign(match(metabolite.list, metabolite.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator = 1 - tag.indicator 
  N = length(metabolite.list) 
  Nh = length(metabolite.set) 
  Nm =  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector = rep(1, N)
  }
  alpha = weighted.score.type
  correl.vector[which(is.na(correl.vector))] = 0
  correl.vector = abs(correl.vector**alpha)
  sum.correl.tag    = sum(na.omit(correl.vector[tag.indicator == 1]))
  norm.tag    = 1.0/sum.correl.tag
  norm.no.tag = 1.0/Nm
  RES = cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES = max(RES)
  min.ES = min(RES)
  if (max.ES > - min.ES) {
    ES = signif(max.ES, digits = 5)
    arg.ES = which.max(RES)
  } else {
    ES = signif(min.ES, digits=5)
    arg.ES = which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

MSEA.EnrichmentScore2 = function(metabolite.list, metabolite.set, weighted.score.type = 1, correl.vector = NULL) {  
  # Computes the weighted MSEA score of metabolite.set in metabolite.list. It is the same calculation as in 
  # MSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  N = length(metabolite.list) 
  Nh = length(metabolite.set) 
  Nm =  N - Nh 
  
  loc.vector = vector(length=N, mode="numeric")
  peak.res.vector = vector(length=Nh, mode="numeric")
  valley.res.vector = vector(length=Nh, mode="numeric")
  tag.correl.vector = vector(length=Nh, mode="numeric")
  tag.loc.vector = vector(length=Nh, mode="numeric")
  tag.diff.vector = vector(length=Nh, mode="numeric")
  
  loc.vector[metabolite.list] = seq(1, N)
  tag.loc.vector = loc.vector[metabolite.set]
  tag.loc.vector = sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector = rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector = correl.vector[tag.loc.vector]
    tag.correl.vector = abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector = correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector = abs(tag.correl.vector)
  } else {
    tag.correl.vector = correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector = abs(tag.correl.vector)
  }
  tag.correl.vector[is.na(tag.correl.vector)] = 1
  
  norm.tag = 1.0/sum(tag.correl.vector)
  tag.correl.vector = tag.correl.vector * norm.tag
  norm.no.tag = 1.0/Nm
  
  tag.diff.vector[1] = (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] = tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector = tag.diff.vector * norm.no.tag
  peak.res.vector = cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector = peak.res.vector - tag.correl.vector
  max.ES = max(peak.res.vector)
  min.ES = min(valley.res.vector)
  ES = signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
}

MSEA.HeatMapPlot = function(V, row.names = F, col.labels, col.classes, col.names = F, main = " ", xlab=" ", ylab=" ") {
  # Plots a heatmap "pinkogram" of a metabolite expression matrix including phenotype vector and metabolite, sample and phenotype labels
  n.rows = length(V[,1])
  n.cols = length(V[1,])
  row.mean = apply(V, MARGIN=1, function(i) mean(na.omit(i)))
  row.sd = apply(V, MARGIN=1, function(i) sd(na.omit(i)))
  row.n = length(V[,1])
  for (i in 1:n.rows) {
    if (row.sd[i] == 0) {V[i,] = 0} else {V[i,] = (V[i,] - row.mean[i])/(0.5 * row.sd[i])}
    V[i,] = ifelse(V[i,] < -6, -6, V[i,])
    V[i,] = ifelse(V[i,] > 6, 6, V[i,])
  }
  # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype).
  mycol = c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")
  mid.range.V = mean(range(V)) - 0.1
  heatm = matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows,] = V[seq(n.rows, 1, -1),]
  heatm[n.rows + 1,] = ifelse(col.labels == 0, 7, -7)
  image(1:n.cols, 1:(n.rows + 1), t(heatm), col=mycol, axes=FALSE, main=main, xlab= xlab, ylab=ylab)
  
  if (length(row.names) > 1) {
    numC = nchar(row.names)
    size.row.char = 35/(n.rows + 5)
    size.col.char = 25/(n.cols + 5)
    maxl = floor(n.rows/1.6)
    for (i in 1:n.rows) {row.names[i] = substr(row.names[i], 1, maxl)}
    row.names = c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at=1:(n.rows + 1), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
  }
  if (length(col.names) > 1) {
    axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
  }
  C = split(col.labels, col.labels)
  class1.size = length(C[[1]])
  class2.size = length(C[[2]])
  axis(3, at=c(floor(class1.size/2),class1.size + floor(class2.size/2)), labels=col.classes, tick=FALSE, las = 1, cex.axis=1.25, font.axis=2, line=-1)
  return()
}

MSEA.ReadClsFile = function(file = "NULL") { 
  # Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a metabolite expression file (RES or GCT format)
  cls.cont = readLines(file)
  num.lines = length(cls.cont)
  class.list = unlist(strsplit(cls.cont[[3]], " "))
  s = length(class.list)
  t = table(class.list)
  l = length(t)
  phen = vector(length=l, mode="character")
  phen.label = vector(length=l, mode="numeric")
  class.v = vector(length=s, mode="numeric")
  for (i in 1:l) {
    phen[i] = noquote(names(t)[i])
    phen.label[i] = i - 1
  }
  for (i in 1:s) {
    for (j in 1:l) {
      if (class.list[i] == phen[j]) {class.v[i] = phen.label[j]}
    }
  }
  return(list(phen = phen, class.v = class.v))
}

# ----------------------------------------------------------------------------------------
# Main MSEA Analysis Function that implements the entire methodology
# This is a methodology for the analysis of global molecular profiles called Metabolite Set Enrichment Analysis (MSEA). It determines 
# whether an a priori defined set of metabolites shows statistically significant, concordant differences between two biological 
# states (e.g. phenotypes). MSEA operates on all metabolites from an experiment, rank ordered by the signal to noise ratio and 
# determines whether members of an a priori defined metabolite set are nonrandomly distributed towards the top or bottom of the 
# list and thus may correspond to an important biological process. To assess significance the program uses an empirical 
# permutation procedure to test deviation from random that preserves correlations between metabolites. 
#
# For details see Subramanian et al 2005
MSEA = function(input.ds, input.cls, met.db, output.directory = "", doc.string = "MSEA.analysis", 
                reshuffling.type = "sample.labels", nperm = 1000, weighted.score.type = 1, 
                nom.p.val.threshold = -1, fwer.p.val.threshold = -1, fdr.q.val.threshold = -1, 
                topmet = 10, adjust.FDR.q.val = F, met.size.threshold.min = 5, met.size.threshold.max = 100, 
                preproc.type = 0, random.seed = 123456, perm.type = 0, 
                fraction = 1.0, replace = F) {
  # Inputs:
  #   input.ds: Input metabolite expression dataset file in GCT format 
  #   input.cls:  Input class vector (phenotype) file in CLS format 
  #   met.file: Metabolite set database in GMT format 
  #   output.directory: Directory where to store output and results (default: .) 
  #   doc.string:  Documentation string used as a prefix to name result files (default: "MSEA.analysis") 
  #   reshuffling.type: Type of permutation reshuffling: "sample.labels" or "metabolite.labels" (default: "sample.labels") 
  #   nperm: Number of random permutations (default: 1000) 
  #   weighted.score.type: Enrichment correlation-based weighting: 0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1) 
  #   nom.p.val.threshold: Significance threshold for nominal p-vals for metabolite sets (default: -1, no thres) 
  #   fwer.p.val.threshold: Significance threshold for FWER p-vals for metabolite sets (default: -1, no thres) 
  #   fdr.q.val.threshold: Significance threshold for FDR q-vals for metabolite sets (default: 0.25) 
  #   topmet: Besides those passing test, number of top scoring metabolite sets used for detailed reports (default: 10) 
  #   adjust.FDR.q.val: Adjust the FDR q-vals (default: F) 
  #   met.size.threshold.min: Minimum size (in metabolites) for database metabolite sets to be considered (default: 25) 
  #   met.size.threshold.max: Maximum size (in metabolites) for database metabolite sets to be considered (default: 500) 
  #   preproc.type: Preprocessing normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (default: 0) 
  #   random.seed: Random number metaboliterator seed. (default: 123456) 
  #   perm.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
  #   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
  #   replace: Resampling mode (replacement or not replacement). For experts only (default: F) 
  #   use.fast.enrichment.routine: if true it uses a faster version to compute random perm. enrichment "MSEA.EnrichmentScore2"  
  #
  #   Output:
  #    The results of the method are stored in the "output.directory" specified by the user as part of the input parameters. 
  #      The results files are:
  #    - Two tab-separated global result text files (one for each phenotype). These files are labeled according to the doc 
  #      string prefix and the phenotype name from the CLS file: <doc.string>.SUMMARY.RESULTS.REPORT.<phenotype>.txt
  #    - One set of global plots. They include a.- metabolite list correlation profile, b.- global observed and null densities, c.- heat map 
  #      for the entire sorted dataset, and d.- p-values vs. NES plot. These plots are in a single JPEG file named 
  #      <doc.string>.global.plots.<phenotype>.jpg. When the program is run interactively these plots appear on a window in the R GUI.
  #    - A variable number of tab-separated metabolite result text files according to how many sets pass any of the significance thresholds 
  #      ("nom.p.val.threshold," "fwer.p.val.threshold," and "fdr.q.val.threshold") and how many are specified in the "topmet" 
  #      parameter. These files are named: <doc.string>.<metabolite set name>.report.txt. 
  #   - A variable number of metabolite set plots (one for each metabolite set report file). These plots include a.- Metabolite set running enrichment
  #      "mountain" plot, b.- metabolite set null distribution and c.- heat map for metabolites in the metabolite set. These plots are stored in a 
  #      single JPEG file named <doc.string>.<metabolite set name>.jpg.
  # The format (columns) for the global result files is as follows.
  # MS : Metabolite set name.
  # SIZE : Size of the set in metabolites.
  # SOURCE : Set definition or source.
  # ES : Enrichment score.
  # NES : Normalized (multiplicative rescaling) normalized enrichment score.
  # NOM p-val : Nominal p-value (from the null distribution of the metabolite set).
  # FDR q-val: False discovery rate q-values
  # FWER p-val: Family wise error rate p-values.
  # Tag %: Percent of metabolite set before running enrichment peak.
  # Metabolite %: Percent of metabolite list before running enrichment peak.
  # Signal : enrichment signal strength.
  # FDR (median): FDR q-values from the median of the null distributions.
  # glob.p.val: P-value using a global statistic (number of sets above the set's NES).
  # 
  # The rows are sorted by the NES values (from maximum positive or negative NES to minimum)
  # The format (columns) for the metabolite set result files is as follows.
  # #: Metabolite number in the (sorted) metabolite set
  # METABOLITE : metabolite name. For example the probe accession number, metabolite symbol or the metabolite identifier gin the dataset.
  # SYMBOL : metabolite symbol from the metabolite annotation file.
  # DESC : metabolite description (title) from the metabolite annotation file.
  # LIST LOC : location of the metabolite in the sorted metabolite list.
  # S2N : signal to noise ratio (correlation) of the metabolite in the metabolite list.
  # RES : value of the running enrichment score at the metabolite location.
  # CORE_ENRICHMENT: is this metabolite is the "core enrichment" section of the list? Yes or No variable specifying in the metabolite location is before (positive ES) or after (negative ES) the running enrichment peak.
  # 
  # The rows are sorted by the metabolite location in the metabolite list.
  # The function call to MSEA returns a  two element list containing the two global result reports as data frames ($report1, $report2).
  # 
  # results1: Global output report for first phenotype 
  # result2:  Global putput report for second phenotype
  print(" *** Running MSEA Analysis...")
  # Copy input parameters to log file
  filename = paste(output.directory, doc.string, "_params.txt", sep="", collapse="")  
  time.string = as.character(as.POSIXlt(Sys.time(),"GMT"))
  write(paste("Run of MSEA on ", time.string), file=filename)
  write(paste("input.ds=", input.ds, sep=" "), file=filename, append=T)
  write(paste("input.cls=", input.cls, sep=" "), file=filename, append=T) 
  write(paste("met.db=", met.db, sep=" "), file=filename, append=T) 
  write(paste("output.directory =", output.directory, sep=" "), file=filename, append=T) 
  write(paste("doc.string = ", doc.string, sep=" "), file=filename, append=T) 
  write(paste("reshuffling.type =", reshuffling.type, sep=" "), file=filename, append=T) 
  write(paste("nperm =", nperm, sep=" "), file=filename, append=T) 
  write(paste("weighted.score.type =", weighted.score.type, sep=" "), file=filename, append=T) 
  write(paste("nom.p.val.threshold =", nom.p.val.threshold, sep=" "), file=filename, append=T) 
  write(paste("fwer.p.val.threshold =", fwer.p.val.threshold, sep=" "), file=filename, append=T) 
  write(paste("fdr.q.val.threshold =", fdr.q.val.threshold, sep=" "), file=filename, append=T) 
  write(paste("topmet =", topmet, sep=" "), file=filename, append=T)
  write(paste("adjust.FDR.q.val =", adjust.FDR.q.val, sep=" "), file=filename, append=T) 
  write(paste("met.size.threshold.min =", met.size.threshold.min, sep=" "), file=filename, append=T) 
  write(paste("met.size.threshold.max =", met.size.threshold.max, sep=" "), file=filename, append=T) 
  write(paste("preproc.type =", preproc.type, sep=" "), file=filename, append=T) 
  write(paste("random.seed =", random.seed, sep=" "), file=filename, append=T) 
  write(paste("perm.type =", perm.type, sep=" "), file=filename, append=T) 
  write(paste("fraction =", fraction, sep=" "), file=filename, append=T) 
  write(paste("replace =", replace, sep=" "), file=filename, append=T)
  
  # Start of MSEA methodology 
  # Read input data matrix
  set.seed(seed=random.seed, kind = NULL)
  adjust.param = 0.5
  time1 = proc.time()
  dataset = read.table(input.ds, sep="\t", check.names=FALSE)
  metabolite.labels = row.names(dataset)
  sample.names = colnames(dataset)
  A = data.matrix(dataset)
  dim(A) 
  cols = length(A[1,])
  rows = length(A[,1])
  # preproc.type control the type of pre-processing: threshold, variation filter, normalization
  if (preproc.type == 1) {  # Column normalize (Z-score)
    A = MSEA.NormalizeCols(A)
  } else if (preproc.type == 2) { # Column (rank) and row (Z-score) normalize 
    for (j in 1:cols) {A[,j] = rank(A[,j])}
    A = MSEA.NormalizeRows(A)
  } else if (preproc.type == 3) { # Column (rank) norm.
    for (j in 1:cols) {A[,j] = rank(A[,j])}
  }
  # Read input class vector
  CLS = MSEA.ReadClsFile(file=input.cls)
  class.labels = CLS$class.v
  class.phen = CLS$phen
  phen1 = class.phen[1]
  phen2 = class.phen[2]
  # sort samples according to phenotype
  col.index = order(class.labels, decreasing=F)
  class.labels = class.labels[col.index]
  sample.names = sample.names[col.index]
  for (j in 1:rows) {A[j, ] = A[j, col.index]}
  names(A) = sample.names
  
  # Read input metabolite set database
  temp = readLines(met.db)
  max.Ng = length(temp)
  temp.size.G = vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] = length(unlist(strsplit(temp[[i]], "    "))) - 2
  }
  
  max.size.G = max(temp.size.G)      
  gs = matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names = vector(length = max.Ng, mode = "character")
  temp.desc = vector(length = max.Ng, mode = "character")
  met.count = 1
  for (i in 1:max.Ng) {
    metabolite.set.size = length(unlist(strsplit(temp[[i]], "    "))) - 2
    met.line = noquote(unlist(strsplit(temp[[i]], "    ")))
    metabolite.set.name = met.line[1] 
    metabolite.set.desc = met.line[2] 
    metabolite.set.tags = vector(length = metabolite.set.size, mode = "character")
    for (j in 1:metabolite.set.size) {
      metabolite.set.tags[j] = trimws(met.line[j + 2])
    } 
    existing.set = is.element(metabolite.set.tags, metabolite.labels)
    set.size = length(existing.set[existing.set == T])
    if ((set.size < met.size.threshold.min) || (set.size > met.size.threshold.max)) next
    temp.size.G[met.count] = set.size
    gs[met.count,] = c(metabolite.set.tags[existing.set], rep("null", max.size.G - temp.size.G[met.count]))
    temp.names[met.count] = metabolite.set.name
    temp.desc[met.count] = metabolite.set.desc
    met.count = met.count + 1
  } 
  Ng = met.count - 1
  met.names = vector(length = Ng, mode = "character")
  met.desc = vector(length = Ng, mode = "character")
  size.G = vector(length = Ng, mode = "numeric") 
  met.names = temp.names[1:Ng]
  met.desc = temp.desc[1:Ng] 
  size.G = temp.size.G[1:Ng]
  
  N = length(A[,1])
  Ns = length(A[1,])
  print(c("Number of metabolites in dataset:", N))
  print(sprintf("Number of Metabolite Pathway Sets between size %d-%d: %d", met.size.threshold.min, met.size.threshold.max, Ng))
  print(c("Number of samples:", Ns))
  print(c("Original number of Metabolite Pathway Sets:", max.Ng))
  print(c("Maximum metabolite pathway set size:", max.size.G))
  
  # Read metabolite and metabolite set annotations if metabolite annotation file was provided
  all.metabolite.descs = vector(length = N, mode ="character") 
  all.metabolite.symbols = vector(length = N, mode ="character") 
  all.met.descs = vector(length = Ng, mode ="character") 
  for (i in 1:N) {
    all.metabolite.descs[i] = metabolite.labels[i]
    all.metabolite.symbols[i] = metabolite.labels[i]
  }
  for (i in 1:Ng) {
    all.met.descs[i] = met.desc[i]
  }
  Obs.indicator = matrix(nrow= Ng, ncol=N)
  Obs.RES = matrix(nrow= Ng, ncol=N)
  Obs.ES = vector(length = Ng, mode = "numeric")
  Obs.arg.ES = vector(length = Ng, mode = "numeric")
  Obs.ES.norm = vector(length = Ng, mode = "numeric")
  time2 = proc.time()
  
  # MSEA methodology
  # Compute observed and random permutation metabolite rankings
  obs.s2n = vector(length=N, mode="numeric")
  signal.strength = vector(length=Ng, mode="numeric")
  tag.frac = vector(length=Ng, mode="numeric")
  metabolite.frac = vector(length=Ng, mode="numeric")
  coherence.ratio = vector(length=Ng, mode="numeric")
  obs.phi.norm = matrix(nrow = Ng, ncol = nperm)
  correl.matrix = matrix(nrow = N, ncol = nperm)
  obs.correl.matrix = matrix(nrow = N, ncol = nperm)
  order.matrix = matrix(nrow = N, ncol = nperm)
  obs.order.matrix = matrix(nrow = N, ncol = nperm)
  
  nperm.per.call = 100
  n.groups = nperm %/% nperm.per.call
  n.rem = nperm %% nperm.per.call
  n.perms = c(rep(nperm.per.call, n.groups), n.rem)
  n.ends = cumsum(n.perms)
  n.starts = n.ends - n.perms + 1
  if (n.rem == 0) {n.tot = n.groups} else {n.tot = n.groups + 1}
  
  for (nk in 1:n.tot) {
    call.nperm = n.perms[nk]
    print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", n.starts[nk], "--", n.ends[nk], sep=" "))
    if (sum(CLS$class.v)==1) {
      O = MSEA.MetaboliteRanking_SingleProfile(A, class.labels, metabolite.labels, call.nperm, permutation.type = perm.type, sigma.correction = "MetaboliteCluster", fraction=fraction, replace=replace)
    } else {
      O = MSEA.MetaboliteRanking(A, class.labels, metabolite.labels, call.nperm, permutation.type = perm.type, sigma.correction = "MetaboliteCluster", fraction=fraction, replace=replace)
    }
    order.matrix[,n.starts[nk]:n.ends[nk]] = O$order.matrix
    obs.order.matrix[,n.starts[nk]:n.ends[nk]] = O$obs.order.matrix
    correl.matrix[,n.starts[nk]:n.ends[nk]] = O$s2n.matrix
    obs.correl.matrix[,n.starts[nk]:n.ends[nk]] = O$obs.s2n.matrix
    rm(O)
  }
  obs.s2n = apply(obs.correl.matrix, 1, function(i) median(na.omit(i)))  # using median to assign enrichment scores
  obs.index = order(obs.s2n, decreasing=T)            
  obs.s2n   = sort(obs.s2n, decreasing=T, na.last = TRUE)            
  
  obs.metabolite.labels = metabolite.labels[obs.index]       
  obs.metabolite.descs = all.metabolite.descs[obs.index]       
  obs.metabolite.symbols = all.metabolite.symbols[obs.index]       
  for (r in 1:nperm) {correl.matrix[, r] = correl.matrix[order.matrix[,r], r]}
  for (r in 1:nperm) {obs.correl.matrix[, r] = obs.correl.matrix[obs.order.matrix[,r], r]}
  metabolite.list2 = obs.index
  for (i in 1:Ng) {
    print(paste("Computing observed enrichment for metabolite set:", i, met.names[i], sep=" ")) 
    metabolite.set = gs[i,gs[i,] != "null"]
    metabolite.set2 = vector(length=length(metabolite.set), mode = "numeric")
    metabolite.set2 = match(metabolite.set, metabolite.labels)
    MSEA.results = MSEA.EnrichmentScore(metabolite.list=metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
    Obs.ES[i] = MSEA.results$ES
    Obs.arg.ES[i] = MSEA.results$arg.ES
    Obs.RES[i,] = MSEA.results$RES
    Obs.indicator[i,] = MSEA.results$indicator
    if (Obs.ES[i] >= 0) {  # compute signal strength
      tag.frac[i] = sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.G[i]
      metabolite.frac[i] = Obs.arg.ES[i]/N
    } else {
      tag.frac[i] = sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
      metabolite.frac[i] = (N - Obs.arg.ES[i] + 1)/N
    }
    signal.strength[i] = tag.frac[i] * (1 - metabolite.frac[i]) * (N / (N - size.G[i]))
  }
  
  # Compute enrichment for random permutations 
  phi = matrix(nrow = Ng, ncol = nperm)
  phi.norm = matrix(nrow = Ng, ncol = nperm)
  obs.phi = matrix(nrow = Ng, ncol = nperm)
  if (reshuffling.type == "sample.labels") { # reshuffling phenotype labels
    for (i in 1:Ng) {
      print(paste("Computing random permutations' enrichment for metabolite set:", i, met.names[i], sep=" ")) 
      metabolite.set = gs[i,gs[i,] != "null"]
      metabolite.set2 = vector(length=length(metabolite.set), mode = "numeric")
      metabolite.set2 = match(metabolite.set, metabolite.labels)
      for (r in 1:nperm) {
        metabolite.list2 = order.matrix[,r]
        MSEA.results = MSEA.EnrichmentScore2(metabolite.list=metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
        phi[i, r] = MSEA.results$ES
      }
      if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
        for (r in 1:nperm) {
          obs.metabolite.list2 = obs.order.matrix[,r]
          MSEA.results = MSEA.EnrichmentScore2(metabolite.list=obs.metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
          obs.phi[i, r] = MSEA.results$ES
        }
      } else { # if no resampling then compute only one column (and fill the others with the same value)
        obs.metabolite.list2 = obs.order.matrix[,1]
        MSEA.results = MSEA.EnrichmentScore2(metabolite.list=obs.metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
        obs.phi[i, 1] = MSEA.results$ES
        for (r in 2:nperm) {obs.phi[i, r] = obs.phi[i, 1]}
      }
      gc()
    } # if (reshuffling.type == "sample.labels")
  } else if (reshuffling.type == "metabolite.labels") { # reshuffling metabolite labels
    for (i in 1:Ng) {
      metabolite.set = gs[i,gs[i,] != "null"]
      metabolite.set2 = vector(length=length(metabolite.set), mode = "numeric")
      metabolite.set2 = match(metabolite.set, metabolite.labels)
      for (r in 1:nperm) {
        reshuffled.metabolite.labels = sample(1:rows)
        MSEA.results = MSEA.EnrichmentScore2(metabolite.list=reshuffled.metabolite.labels, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=obs.s2n)   
        phi[i, r] = MSEA.results$ES
      }
      if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
        for (r in 1:nperm) {
          obs.metabolite.list2 = obs.order.matrix[,r]
          MSEA.results = MSEA.EnrichmentScore2(metabolite.list=obs.metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
          obs.phi[i, r] = MSEA.results$ES
        }
      } else { # if no resampling then compute only one column (and fill the others with the same value)
        obs.metabolite.list2 = obs.order.matrix[,1]
        MSEA.results = MSEA.EnrichmentScore2(metabolite.list=obs.metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
        obs.phi[i, 1] = MSEA.results$ES
        for (r in 2:nperm) {obs.phi[i, r] = obs.phi[i, 1]}
      }
      gc()
    }
  }
  
  # Compute 3 types of p-values
  # Find nominal p-values       
  print("Computing nominal p-values...")
  p.vals = matrix(1, nrow = Ng, ncol = 2)
  for (i in 1:Ng) {
    pos.phi = NULL
    neg.phi = NULL
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        pos.phi = c(pos.phi, phi[i, j]) 
      } else {
        neg.phi = c(neg.phi, phi[i, j]) 
      }
    }
    ES.value = Obs.ES[i]
    if (ES.value >= 0) {
      p.vals[i, 1] = signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=5)
    } else {
      p.vals[i, 1] = signif(sum(neg.phi <= ES.value)/length(neg.phi), digits=5)
    }
  }
  
  # Find effective size 
  erf = function (x) {2 * pnorm(sqrt(2) * x)}
  KS.mean = function(N) { # KS mean as a function of set size N
    S = 0
    for (k in -100:100) {
      if (k == 0) next
      S = S + 4 * (-1)**(k + 1) * (0.25 * exp(-2 * k * k * N) - sqrt(2 * pi) *  erf(sqrt(2 * N) * k)/(16 * k * sqrt(N)))
    }
    return(abs(S))
  }
  # Rescaling normalization for each metabolite set null
  print("Computing rescaling normalization for each metabolite set null...")
  for (i in 1:Ng) {
    pos.phi = NULL
    neg.phi = NULL
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        pos.phi = c(pos.phi, phi[i, j]) 
      } else {
        neg.phi = c(neg.phi, phi[i, j]) 
      }
    }
    pos.m = mean(pos.phi)
    neg.m = mean(abs(as.numeric(neg.phi)))
    
    pos.phi = pos.phi/pos.m
    neg.phi = neg.phi/neg.m
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        phi.norm[i, j] = phi[i, j]/pos.m
      } else {
        phi.norm[i, j] = phi[i, j]/neg.m
      }
    }
    for (j in 1:nperm) {
      if (obs.phi[i, j] >= 0) {
        obs.phi.norm[i, j] = obs.phi[i, j]/pos.m
      } else {
        obs.phi.norm[i, j] = obs.phi[i, j]/neg.m
      }
    }
    if (Obs.ES[i] >= 0) {
      Obs.ES.norm[i] = Obs.ES[i]/pos.m
    } else {
      Obs.ES.norm[i] = Obs.ES[i]/neg.m
    }
  }
  
  # Compute FWER p-vals
  print("Computing FWER p-values...")
  max.ES.vals.p = NULL
  max.ES.vals.n = NULL
  for (j in 1:nperm) {
    pos.phi = NULL
    neg.phi = NULL
    for (i in 1:Ng) {
      if (phi.norm[i, j] >= 0) {
        pos.phi = c(pos.phi, phi.norm[i, j]) 
      } else {
        neg.phi = c(neg.phi, phi.norm[i, j]) 
      }
    }
    if (length(pos.phi) > 0) {
      max.ES.vals.p = c(max.ES.vals.p, max(pos.phi))
    }
    if (length(neg.phi) > 0) {
      max.ES.vals.n = c(max.ES.vals.n, min(neg.phi))
    }
  }
  for (i in 1:Ng) {
    ES.value = Obs.ES.norm[i]
    if (Obs.ES.norm[i] >= 0) {
      p.vals[i, 2] = signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), digits=5)
    } else {
      p.vals[i, 2] = signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), digits=5)
    }
  }
  
  # Compute FDRs 
  print("Computing FDR q-values...")
  NES = vector(length=Ng, mode="numeric")
  phi.norm.mean  = vector(length=Ng, mode="numeric")
  obs.phi.norm.mean  = vector(length=Ng, mode="numeric")
  phi.norm.median  = vector(length=Ng, mode="numeric")
  obs.phi.norm.median  = vector(length=Ng, mode="numeric")
  phi.norm.mean  = vector(length=Ng, mode="numeric")
  obs.phi.mean  = vector(length=Ng, mode="numeric")
  FDR.mean = vector(length=Ng, mode="numeric")
  FDR.median = vector(length=Ng, mode="numeric")
  phi.norm.median.d = vector(length=Ng, mode="numeric")
  obs.phi.norm.median.d = vector(length=Ng, mode="numeric")
  
  Obs.ES.index = order(Obs.ES.norm, decreasing=T)
  Orig.index = seq(1, Ng)
  Orig.index = Orig.index[Obs.ES.index]
  Orig.index = order(Orig.index, decreasing=F)
  Obs.ES.norm.sorted = Obs.ES.norm[Obs.ES.index]
  met.names.sorted = met.names[Obs.ES.index]
  
  for (k in 1:Ng) {
    NES[k] = Obs.ES.norm.sorted[k]
    ES.value = NES[k]
    count.col = vector(length=nperm, mode="numeric")
    obs.count.col = vector(length=nperm, mode="numeric")
    for (i in 1:nperm) {
      phi.vec = phi.norm[,i]
      obs.phi.vec = obs.phi.norm[,i]
      if (ES.value >= 0) {
        count.col.norm = sum(phi.vec >= 0)
        obs.count.col.norm = sum(obs.phi.vec >= 0)
        count.col[i] = ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
        obs.count.col[i] = ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
      } else {
        count.col.norm = sum(phi.vec < 0)
        obs.count.col.norm = sum(obs.phi.vec < 0)
        count.col[i] = ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
        obs.count.col[i] = ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
      }
    }
    phi.norm.mean[k] = mean(count.col)
    obs.phi.norm.mean[k] = mean(obs.count.col)
    phi.norm.median[k] = median(count.col)
    obs.phi.norm.median[k] = median(obs.count.col)
    FDR.mean[k] = ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
    FDR.median[k] = ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
  }
  FDR.mean[which(is.na(FDR.mean))] = 1
  FDR.median[which(is.na(FDR.median))] = 1
  
  # adjust q-values
  if (adjust.FDR.q.val == T) {
    pos.nes = length(NES[NES >= 0])
    min.FDR.mean = FDR.mean[pos.nes]
    min.FDR.median = FDR.median[pos.nes]
    for (k in seq(pos.nes - 1, 1, -1)) {
      if (FDR.mean[k] < min.FDR.mean) {min.FDR.mean = FDR.mean[k]}
      if (min.FDR.mean < FDR.mean[k]) {FDR.mean[k] = min.FDR.mean}
    }
    neg.nes = pos.nes + 1
    min.FDR.mean = FDR.mean[neg.nes]
    min.FDR.median = FDR.median[neg.nes]
    for (k in seq(neg.nes + 1, Ng)) {
      if (FDR.mean[k] < min.FDR.mean) {min.FDR.mean = FDR.mean[k]}
      if (min.FDR.mean < FDR.mean[k]) {FDR.mean[k] = min.FDR.mean}
    }
  }
  obs.phi.norm.mean.sorted = obs.phi.norm.mean[Orig.index]
  phi.norm.mean.sorted = phi.norm.mean[Orig.index]
  FDR.mean.sorted = FDR.mean[Orig.index]
  FDR.median.sorted = FDR.median[Orig.index]
  
  # Compute global statistic
  glob.p.vals = vector(length=Ng, mode="numeric")
  NULL.pass = vector(length=nperm, mode="numeric")
  OBS.pass = vector(length=nperm, mode="numeric")
  
  for (k in 1:Ng) {
    NES[k] = Obs.ES.norm.sorted[k]
    if (NES[k] >= 0) {
      for (i in 1:nperm) {
        NULL.pos = sum(phi.norm[,i] >= 0)            
        NULL.pass[i] = ifelse(NULL.pos > 0, sum(phi.norm[,i] >= NES[k])/NULL.pos, 0)
        OBS.pos = sum(obs.phi.norm[,i] >= 0)
        OBS.pass[i] = ifelse(OBS.pos > 0, sum(obs.phi.norm[,i] >= NES[k])/OBS.pos, 0)
      }
    } else {
      for (i in 1:nperm) {
        NULL.neg = sum(phi.norm[,i] < 0)
        NULL.pass[i] = ifelse(NULL.neg > 0, sum(phi.norm[,i] <= NES[k])/NULL.neg, 0)
        OBS.neg = sum(obs.phi.norm[,i] < 0)
        OBS.pass[i] = ifelse(OBS.neg > 0, sum(obs.phi.norm[,i] <= NES[k])/OBS.neg, 0)
      }
    }
    glob.p.vals[k] = sum(NULL.pass >= mean(OBS.pass))/nperm
  }
  glob.p.vals.sorted = glob.p.vals[Orig.index]
  
  # Produce results report
  print("Producing result tables and plots...")
  Obs.ES = signif(Obs.ES, digits=5)
  Obs.ES.norm = signif(Obs.ES.norm, digits=5)
  p.vals = signif(p.vals, digits=4)
  signal.strength = signif(signal.strength, digits=3)
  tag.frac = signif(tag.frac, digits=3)
  metabolite.frac = signif(metabolite.frac, digits=3)
  FDR.mean.sorted = signif(FDR.mean.sorted, digits=5)
  FDR.median.sorted =  signif(FDR.median.sorted, digits=5)
  glob.p.vals.sorted = signif(glob.p.vals.sorted, digits=5)
  
  report = data.frame(cbind(met.names, size.G, all.met.descs, Obs.ES, Obs.ES.norm, p.vals[,1], FDR.mean.sorted, p.vals[,2], tag.frac, metabolite.frac, signal.strength, FDR.median.sorted, glob.p.vals.sorted))
  names(report) = c("MS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Metabolite %", "Signal", "FDR (median)", "glob.p.val")
  report2 = report
  report.index2 = order(Obs.ES.norm, decreasing=T)
  for (i in 1:Ng) {report2[i,] = report[report.index2[i],]}   
  report3 = report
  report.index3 = order(Obs.ES.norm, decreasing=F)
  for (i in 1:Ng) {report3[i,] = report[report.index3[i],]}   
  report.phen1 = report2
  report.phen2 = report3
  
  if (output.directory != "")  {
    filename = paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", phen1,".txt", sep="", collapse="")
    write.table(report.phen1, file = filename, quote=F, row.names=F, sep = "\t")
  }
  # Global plots
  if (output.directory != "")  {
    glob.filename = paste(output.directory, doc.string, ".global.plots.pdf", sep="", collapse="")
    pdf(file=glob.filename, height = 10, width = 10)
  }
  nf = layout(matrix(c(1,2,3,4), 2, 2, byrow=T), c(1,1), c(1,1), TRUE)
  # plot S2N correlation profile
  location = 1:N
  max.corr = max(obs.s2n)
  min.corr = min(obs.s2n)
  
  x = plot(location, obs.s2n, ylab = "Signal to Noise Ratio (S2N)", xlab = "Metabolite List Location", main = "Metabolite List Correlation (S2N) Profile", type = "l", lwd = 2, cex = 0.9, col = 1)            
  for (i in seq(1, N, 20)) {
    lines(c(i, i), c(0, obs.s2n[i]), lwd = 3, cex = 0.9, col = colors()[12]) # shading of correlation plot
  }
  x = points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, col = 1)            
  lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1) # zero correlation horizontal line
  temp = order(abs(obs.s2n), decreasing=T)
  arg.correl = temp[N]
  lines(c(arg.correl, arg.correl), c(min.corr, 0.7*max.corr), lwd = 2, lty = 3, cex = 0.9, col = 1) # zero correlation vertical line
  
  area.bias = signif(100*(sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N]))/sum(abs(obs.s2n[1:N])), digits=3)
  area.phen = ifelse(area.bias >= 0, phen1, phen2)
  delta.string = paste("Corr. Area Bias to \"", area.phen, "\" =", abs(area.bias), "%", sep="", collapse="")
  zero.crossing.string = paste("Zero Crossing at location ", arg.correl, " (",  signif(100*arg.correl/N, digits=3), " %)")
  leg.txt = c(delta.string, zero.crossing.string)
  legend(x=N/10, y=max.corr, bty="n", bg = "white", legend=leg.txt, cex = 0.9)
  
  leg.txt = paste("\"", phen1, "\" ", sep="", collapse="")
  text(x=1, y=-0.05*max.corr, adj = c(0, 1), labels=leg.txt, cex = 0.9)
  leg.txt = paste("\"", phen2, "\" ", sep="", collapse="")
  text(x=N, y=0.05*max.corr, adj = c(1, 0), labels=leg.txt, cex = 0.9)
  
  if (Ng > 1) { # make these plots only if there are multiple metabolite sets.
    # compute plots of actual (weighted) null and observed
    phi.densities.pos = matrix(0, nrow=512, ncol=nperm)
    phi.densities.neg = matrix(0, nrow=512, ncol=nperm)
    obs.phi.densities.pos = matrix(0, nrow=512, ncol=nperm)
    obs.phi.densities.neg = matrix(0, nrow=512, ncol=nperm)
    phi.density.mean.pos = vector(length=512, mode = "numeric")
    phi.density.mean.neg = vector(length=512, mode = "numeric")
    obs.phi.density.mean.pos = vector(length=512, mode = "numeric")
    obs.phi.density.mean.neg = vector(length=512, mode = "numeric")
    phi.density.median.pos = vector(length=512, mode = "numeric")
    phi.density.median.neg = vector(length=512, mode = "numeric")
    obs.phi.density.median.pos = vector(length=512, mode = "numeric")
    obs.phi.density.median.neg = vector(length=512, mode = "numeric")
    x.coor.pos =  vector(length=512, mode = "numeric")
    x.coor.neg =  vector(length=512, mode = "numeric")
    for (i in 1:nperm) {
      pos.phi = phi.norm[phi.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp = density(pos.phi, adjust=adjust.param, n = 512, from=0, to=3.5)
      } else {
        temp = list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      phi.densities.pos[, i] = temp$y
      norm.factor = sum(phi.densities.pos[, i])
      phi.densities.pos[, i] = phi.densities.pos[, i]/norm.factor
      if (i == 1) {
        x.coor.pos = temp$x
      }
      neg.phi = phi.norm[phi.norm[, i] < 0, i]
      if (length(neg.phi) > 2) {
        temp = density(neg.phi, adjust=adjust.param, n = 512, from=-3.5, to=0)
      } else {
        temp = list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      phi.densities.neg[, i] = temp$y
      norm.factor = sum(phi.densities.neg[, i])
      phi.densities.neg[, i] = phi.densities.neg[, i]/norm.factor
      if (i == 1) {
        x.coor.neg = temp$x
      }
      pos.phi = obs.phi.norm[obs.phi.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp = density(pos.phi, adjust=adjust.param, n = 512, from=0, to=3.5)
      } else {
        temp = list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      obs.phi.densities.pos[, i] = temp$y
      norm.factor = sum(obs.phi.densities.pos[, i])
      obs.phi.densities.pos[, i] = obs.phi.densities.pos[, i]/norm.factor
      neg.phi = obs.phi.norm[obs.phi.norm[, i] < 0, i]
      if (length(neg.phi)> 2) {  
        temp = density(neg.phi, adjust=adjust.param, n = 512, from=-3.5, to=0)
      } else {
        temp = list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      obs.phi.densities.neg[, i] = temp$y
      norm.factor = sum(obs.phi.densities.neg[, i])
      obs.phi.densities.neg[, i] = obs.phi.densities.neg[, i]/norm.factor
    }
    phi.density.mean.pos = apply(phi.densities.pos, 1, mean)
    phi.density.mean.neg = apply(phi.densities.neg, 1, mean)
    obs.phi.density.mean.pos = apply(obs.phi.densities.pos, 1, mean)
    obs.phi.density.mean.neg = apply(obs.phi.densities.neg, 1, mean)
    phi.density.median.pos = apply(phi.densities.pos, 1, median)
    phi.density.median.neg = apply(phi.densities.neg, 1, median)
    obs.phi.density.median.pos = apply(obs.phi.densities.pos, 1, median)
    obs.phi.density.median.neg = apply(obs.phi.densities.neg, 1, median)
    
    x = c(x.coor.neg, x.coor.pos)
    x.plot.range = range(x)
    y1 = c(phi.density.mean.neg, phi.density.mean.pos)
    y2 = c(obs.phi.density.mean.neg, obs.phi.density.mean.pos)
    y.plot.range = c(-0.3*max(c(y1, y2)),  max(c(y1, y2)))
    print(c(y.plot.range, max(c(y1, y2)), max(y1), max(y2)))
    plot(x, y1, xlim = x.plot.range, ylim = 1.5*y.plot.range, type = "l", lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", main = "Global Observed and Null Densities (Area Normalized)")
    
    y1.point = y1[seq(1, length(x), 2)]
    y2.point = y2[seq(2, length(x), 2)]
    x1.point = x[seq(1, length(x), 2)]
    x2.point = x[seq(2, length(x), 2)]
    
    points(x, y1, type = "l", lwd = 2, col = colors()[555])
    points(x, y2, type = "l", lwd = 2, col = colors()[29])
    for (i in 1:Ng) {
      col = ifelse(Obs.ES.norm[i] > 0, 2, 3) 
      lines(c(Obs.ES.norm[i], Obs.ES.norm[i]), c(-0.2*max(c(y1, y2)), 0), lwd = 1, lty = 1, col = 1)
    }
    leg.txt = paste("Neg. ES: \"", phen2, " \" ", sep="", collapse="")
    text(x=x.plot.range[1], y=-0.25*max(c(y1, y2)), adj = c(0, 1), labels=leg.txt, cex = 0.9)
    leg.txt = paste(" Pos. ES: \"", phen1, "\" ", sep="", collapse="")
    text(x=x.plot.range[2], y=-0.25*max(c(y1, y2)), adj = c(1, 1), labels=leg.txt, cex = 0.9)
    
    leg.txt = c("Null Density", "Observed Density", "Observed NES values")
    c.vec = c(colors()[555], colors()[29], 1)
    lty.vec = c(1, 1, 1)
    lwd.vec = c(2, 2, 2)
    legend(x=0, y=1.5*y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 0.9)
    B = A[obs.index,]
    if (N > 300) {
      C = rbind(B[1:100,], rep(0, Ns), rep(0, Ns), B[(floor(N/2) - 50 + 1):(floor(N/2) + 50),], rep(0, Ns), rep(0, Ns), B[(N - 100 + 1):N,])
    } 
    rm(B)
    MSEA.HeatMapPlot(V = C, col.labels = class.labels, col.classes = class.phen, main = "Heat Map for Metabolites in Dataset")
    
    # p-vals plot
    nom.p.vals = p.vals[Obs.ES.index,1]
    FWER.p.vals = p.vals[Obs.ES.index,2]
    plot.range = 1.25*range(NES)
    plot(NES, FDR.mean, ylim = c(0, 1), xlim = plot.range, col = 1, bg = 1, type="p", pch = 22, cex = 0.75, xlab = "NES", main = "p-values vs. NES", ylab ="p-val/q-val")
    points(NES, nom.p.vals, type = "p", col = 2, bg = 2, pch = 22, cex = 0.75)
    points(NES, FWER.p.vals, type = "p", col = colors()[577], bg = colors()[577], pch = 22, cex = 0.75)
    leg.txt = c("Nominal p-value", "FWER p-value", "FDR q-value")
    c.vec = c(2, colors()[577], 1)
    pch.vec = c(22, 22, 22)
    legend(x=-0.5, y=0.5, bty="n", bg = "white", legend=leg.txt, pch = pch.vec, col = c.vec, pt.bg = c.vec, cex = 0.9)
    lines(c(min(NES), max(NES)), c(nom.p.val.threshold, nom.p.val.threshold), lwd = 1, lty = 2, col = 2) 
    lines(c(min(NES), max(NES)), c(fwer.p.val.threshold, fwer.p.val.threshold), lwd = 1, lty = 2, col = colors()[577]) 
    lines(c(min(NES), max(NES)), c(fdr.q.val.threshold, fdr.q.val.threshold), lwd = 1, lty = 2, col = 1) 
    dev.off()
  } # if Ng > 1
  
  #----------------------------------------------------------------------------
  # Produce report for each metabolite set passing the nominal, FWER or FDR test or the top topmet in each side
  if (topmet > floor(Ng/2)) {topmet = floor(Ng/2)}
  for (i in 1:Ng) {
    if ((p.vals[i, 1] <= nom.p.val.threshold) || (p.vals[i, 2] <= fwer.p.val.threshold) ||
        (FDR.mean.sorted[i] <= fdr.q.val.threshold) || (is.element(i, c(Obs.ES.index[1:topmet], Obs.ES.index[(Ng - topmet + 1): Ng])))) {
      #  produce report per metabolite set
      kk = 1
      metabolite.number = vector(length = size.G[i], mode = "character")
      metabolite.names = vector(length = size.G[i], mode = "character")
      metabolite.symbols = vector(length = size.G[i], mode = "character")
      metabolite.descs = vector(length = size.G[i], mode = "character")
      metabolite.list.loc = vector(length = size.G[i], mode = "numeric")
      core.enrichment = vector(length = size.G[i], mode = "character")
      metabolite.s2n = vector(length = size.G[i], mode = "numeric")
      metabolite.RES = vector(length = size.G[i], mode = "numeric")
      rank.list = seq(1, N)
      if (Obs.ES[i] >= 0) {
        set.k = seq(1, N, 1)
        phen.tag = phen1
        loc = match(i, Obs.ES.index)
      } else {
        set.k = seq(N, 1, -1)
        phen.tag = phen2
        loc = Ng - match(i, Obs.ES.index) + 1
      }
      
      for (k in set.k) {
        if (Obs.indicator[i, k] == 1) {
          metabolite.number[kk] = kk
          metabolite.names[kk] = obs.metabolite.labels[k]
          metabolite.symbols[kk] = substr(obs.metabolite.symbols[k], 1, 15)
          metabolite.descs[kk] = substr(obs.metabolite.descs[k], 1, 40)
          metabolite.list.loc[kk] = k
          metabolite.s2n[kk] = signif(obs.s2n[k], digits=3)
          metabolite.RES[kk] = signif(Obs.RES[i, k], digits = 3)
          if (Obs.ES[i] >= 0) {
            core.enrichment[kk] = ifelse(metabolite.list.loc[kk] <= Obs.arg.ES[i], "YES", "NO")
          } else {
            core.enrichment[kk] = ifelse(metabolite.list.loc[kk] > Obs.arg.ES[i], "YES", "NO")
          }
          kk = kk + 1
        }
      }
      metabolite.report = data.frame(cbind(metabolite.number, metabolite.names, metabolite.symbols, metabolite.descs, metabolite.list.loc, metabolite.s2n, metabolite.RES, core.enrichment))
      names(metabolite.report) = c("#", "GENE", "SYMBOL", "DESC", "LIST LOC", "S2N", "RES", "CORE_ENRICHMENT")
      if (output.directory != "")  {
        filename = paste(output.directory, doc.string, ".", met.names[i], ".report.", phen.tag, ".", loc, ".txt", sep="", collapse="")
        write.table(metabolite.report, file = filename, quote=F, row.names=F, sep = "\t")
        met.filename = paste(output.directory, doc.string, ".", met.names[i], ".plot.", phen.tag, ".", loc, ".pdf", sep="", collapse="")
        pdf(file=met.filename, height = 6, width = 14)
      }
      
      nf = layout(matrix(c(1,2,3), 1, 3, byrow=T), 1, c(1, 1, 1), TRUE)
      ind = 1:N
      min.RES = min(na.omit(Obs.RES[i,]))
      max.RES = max(na.omit(Obs.RES[i,]))
      if (max.RES < 0.3) max.RES = 0.3
      if (min.RES > -0.3) min.RES = -0.3
      delta = (max.RES - min.RES)*0.50
      min.plot = min.RES - 2*delta
      max.plot = max.RES
      max.corr = max(obs.s2n)
      min.corr = min(obs.s2n)
      Obs.correl.vector.norm = (obs.s2n - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
      zero.corr.line = (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
      col = ifelse(Obs.ES[i] > 0, 2, 4)
      
      # Running enrichment plot
      sub.string = paste("Number of metabolites: ", N, " (in list), ", size.G[i], " (in metabolite set)", sep = "", collapse="")
      main.string = paste("Metabolite Set ", i, ":", met.names[i])
      plot(ind, Obs.RES[i,], main = main.string, sub = sub.string, xlab = "Metabolite List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
      for (j in seq(1, N, 20)) {
        lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
      }
      lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
      lines(c(Obs.arg.ES[i], Obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
      for (j in 1:N) {
        if (Obs.indicator[i, j] == 1) {
          lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
        }
      }
      lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
      lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
      temp = order(abs(obs.s2n), decreasing=T)
      arg.correl = temp[N]
      lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line
      
      leg.txt = paste("\"", phen1, "\" ", sep="", collapse="")
      text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
      leg.txt = paste("\"", phen2, "\" ", sep="", collapse="")
      text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
      adjx = ifelse(Obs.ES[i] > 0, 0, 1)
      
      leg.txt = paste("Peak at ", Obs.arg.ES[i], sep="", collapse="")
      text(x=Obs.arg.ES[i], y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
      leg.txt = paste("Zero crossing at ", arg.correl, sep="", collapse="")
      text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
      
      # nominal p-val histogram
      sub.string = paste("ES =", signif(Obs.ES[i], digits = 3), " NES =", signif(Obs.ES.norm[i], digits=3), "Nom. p-val=", signif(p.vals[i, 1], digits = 3),"FWER=", signif(p.vals[i, 2], digits = 3), "FDR=", signif(FDR.mean.sorted[i], digits = 3))
      temp = density(phi[i,], adjust=adjust.param)
      x.plot.range = range(temp$x)
      y.plot.range = c(-0.125*max(temp$y), 1.5*max(temp$y))
      plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, ylim = y.plot.range, lwd = 2, col = 2, main = "Metabolite Set Null Distribution", xlab = "ES", ylab="P(ES)")
      x.loc = which.min(abs(temp$x - Obs.ES[i]))
      lines(c(Obs.ES[i], Obs.ES[i]), c(0, temp$y[x.loc]), lwd = 2, lty = 1, cex = 1, col = 1)
      lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)
      
      leg.txt = c("Metabolite Set Null Density", "Observed Metabolite Set ES value")
      c.vec = c(2, 1)
      lty.vec = c(1, 1)
      lwd.vec = c(2, 2)
      legend(x=-0.2, y=y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 1.0)
      leg.txt = paste("Neg. ES \"", phen2, "\" ", sep="", collapse="")
      text(x=x.plot.range[1], y=-0.1*max(temp$y), adj = c(0, 0), labels=leg.txt, cex = 1.0)
      leg.txt = paste(" Pos. ES: \"", phen1, "\" ", sep="", collapse="")
      text(x=x.plot.range[2], y=-0.1*max(temp$y), adj = c(1, 0), labels=leg.txt, cex = 1.0)
      
      # create pinkogram for each metabolite set
      kk = 1
      pinko = matrix(0, nrow = size.G[i], ncol = cols)
      pinko.metabolite.names = vector(length = size.G[i], mode = "character")
      for (k in 1:rows) {
        if (Obs.indicator[i, k] == 1) {
          pinko[kk,] = A[obs.index[k],]
          pinko.metabolite.names[kk] = obs.metabolite.symbols[k]
          kk = kk + 1
        }
      }
      MSEA.HeatMapPlot(V = pinko, row.names = pinko.metabolite.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main =" Heat Map for Metabolites in Metabolite Set", xlab=" ", ylab=" ")
      dev.off()
    } # if p.vals thres
  } # loop over metabolite sets
  return(list(report1 = report.phen1, report2 = report.phen2))
}  # end of definition of MSEA.analysis




