
# == title
# Read CLS file which stores the phenotype data
#
# == param
# -file cls file path
# -treatment string of treatment label in cls file
# -control string of control label in cls file
#
# == details
# The CLS file format defines the phenotype data of microarray experiments.
# The first line is the number of samples, number of classes and the third number always be 1.
# These three numbers are seperated by spaces or tabs.
# The second line begins with #. The next two strings usually are the label of the phenotype. 
# The third line is the label of each samples where same label represents 
# the same class. 
#
# The first and the second line is ignored by this function 
# and class labels are taken from the factor of the vector parsed from the third line.
#
# == value
# A ``sampleLabel`` class object
#
# == source
# http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `read.gct`, `sampleLabel`
#
# == example
# \dontrun{
# # P53.cls can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# label = read.cls("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", 
#     treatment="MUT", control="WT")
# }
read.cls = function (file, treatment, control) {
    label = readLines(file, n = -1)
    label = label[3]
    label = unlist(strsplit(label, " |\t"));
    return(sampleLabel(label = label, treatment = treatment, control = control))
}

# == title
# Read GCT format file which stores the expression values
#
# == param
# -file gct file path
#
# == details
# The GCT format is a tab delimited file format that stores the expression value matrix.
# The first line of the file is the version number which always be #1.2.
# The second line is the number of the size of genes and samples, 
# seperated by space, usually for the initiation of reading the expression matrix.
# The third line contains a list of identifiers for the samples associated 
# with each of the columns in the remainder of the file.
# From the fourth line will be the expression value of each gene. 
#
# GCT file is used together with CLS file.
#
# == value
# A matrix of the expression values, with rows correponding to genes and cols to
# samples.
#
# == source
# http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# # expression data stored in a gct format file
# # P53_symbol.gct can be downloaded from
# # http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
# eset = read.gct("http://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53_symbol.gct")
# head(eset)
# }
read.gct = function (file) {
    expr = read.table(file, skip = 2, header = TRUE, sep = "\t", quote = "")
    rownames(expr) = expr[,1]

    checkName = table(expr[,1])
    if(max(checkName) > 1) {
        stop(paste("Genes in gct file should be unique: ", names(which.max(checkName)), sep = " "))
    }
    expr = expr[,-c(1,2)]
    expr = as.matrix(expr)
    
    return(expr)
}

# == title
# Generate data structure of sample labels
#
# == param
# -label sample label vector
# -treatment treatment label
# -control control label
#
# == details
# Since sample label will not be modified in the analysis, this function is used to
# integrate all the label information in one single data object.
#
# == value
# A ``sampleLabel`` class object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == example
# sampleLabel(c("A", "B", "B", "A", "A", "A", "B", "B"), treatment = "A", control = "B")
#
sampleLabel = function (label, treatment, control) {
    if (sum(label == treatment) == 0 | sum(label == control) == 0) {
        stop("Can not find treatment label or control label.")
    }
    res = list(label = label, treatment = treatment, control = control)
    class(res) = "sampleLabel"
    return(res)
}

.treatment = function(sl) {
    return(which(sl$label == sl$treatment))
}

.control = function(sl) {
    return(which(sl$label == sl$control))
}

.permutate = function(sl) {
    sl$label = sample(sl$label, length(sl$label), replace = FALSE)
    return(sl)
}

.factor = function(sl) {
    return(factor(sl$label))
}