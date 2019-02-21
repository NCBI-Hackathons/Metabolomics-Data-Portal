# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SPIA", version = "3.8")
# BiocManager::install("hgu133plus2.db")


library(SPIA)
data(colorectalcancer)
options(digits=3)
head(top)

library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.1,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal=top$ENTREZ

DE_Colorectal[1:10]
ALL_Colorectal[1:10]

