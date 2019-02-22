source("R/getPathwayIgraph.r")
source("R/CePa_v0.6/R/ora_extension_metab.R")
require(CePa)
require(igraph)
require(data.table)


# for(i in 1:length(pc$pathList)) {
  
  # cat("    ", i, "/", length(pc$pathList), ", ", pathway.name[i], "...\n", sep="")
  # 
  # path = pc$pathList[[i]]
  # inter = pc$interactionList[pc$interactionList[, 1] %in% path, 2:3]
  # 
  # pathway = generate.pathway(as.matrix(inter))
  # 
  load("inst/extdata/RData/Arginine-Metabolism.RData")
  
    
  cat("      Calculate node level value and permutate sample labels...\n")

A <- read.delim("data/Wangler2017_EDTA.txt")

A <- tibble::rownames_to_column(A, var="metabolite")


B <- subset(A, metabolite %in% intersect(metabolite, V(ig)$label))

# ora <- cepa.ora.metab(dif = B$metabolite, pmap.path = tools::file_path_as_absolute("inst/extdata/RData"), bk = A$metabolite, pathway = ig, iter=100 )

pathway.names <- c("Arginine-Metabolism", "Ascorbate-Metabolism", "Asp-Glu-Metabolism", 
                   "BCAA-Metabolism", "Benzoate-Metabolism", "Beta-Oxidation", "Bile-Acid-Metabolism", 
                   "Carnitine-Biosynthesis", "Cholesterol-Synthesis", "Creatine-Metabolism", 
                   "DicarboxylicAcid-Metabolism", "Eicosanoids", "Endocannabinoid-Synthesis", 
                   "FattyAcid-Metabolism", "Fibrinogen-Cleavage-Peptides", "GABA-Shunt", 
                   "Galactose-Metabolism", "Glutathione-Metabolism", "Gly-Ser-Thr-Metabolism", 
                   "Glycogen-Metabolism", "Glycolysis", "Glycosylation", "Hemoglobin-Porphyrin-Metabolism", 
                   "Histidine-Metabolism", "Inositol-Metabolism", "Ketone-Bodies", 
                   "Lysine-Catabolism", "Met-Cys-Metabolism", "Mevalonate-Metabolism", 
                   "Nicotinate-Nicotinamide-Metabolism", "Pantothenate-Metabolism", 
                   "Pentose-Phosphate-Metabolism", "Phe-Tyr-Metabolism", "Phospholipid-Metabolism", 
                   "Polyamine-Metabolism", "Proline-Metabolism", "Protein-Degradation", 
                   "Purine-Metabolism", "Pyridoxal-Metabolism", "Pyrimidine-Metabolism", 
                   "Riboflavin-Metabolism", "Secondary-Bile-Acids", "Sorbitol-Glycerol-Metabolism", 
                   "Sphingolipid-Metabolism", "Steroid-Hormone-Biosynthesis", "TCA-Cycle", 
                   "Thyroid-Hormone-Synthesis", "Tryptophan-Metabolism")


# ora <- cepa.ora.metab(dif = B$metabolite,pathway=ig, bk = A$metabolite, iter=100 )
system.time(ora.all <- cepa.ora.metab.all(dif = B$metabolite, pmap.path = tools::file_path_as_absolute("inst/extdata"), bk = A$metabolite, pathway.name =pathway.names, iter=1000, cen = "betweenness" ))
  

ora.table <- cepa_output(ora.all)