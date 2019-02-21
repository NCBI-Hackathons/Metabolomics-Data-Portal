### R code from vignette source 'SPIA.Rnw'

###################################################
### code chunk number 1: SPIA.Rnw:75-79
###################################################
library(SPIA)
data(colorectalcancer)
options(digits=3)
head(top)


###################################################
### code chunk number 2: SPIA.Rnw:99-108
###################################################
library(hgu133plus2.db)
x <- hgu133plus2ENTREZID 
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.1,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal=top$ENTREZ


###################################################
### code chunk number 3: SPIA.Rnw:118-120
###################################################
DE_Colorectal[1:10]
ALL_Colorectal[1:10]


###################################################
### code chunk number 4: SPIA.Rnw:127-133
###################################################
# pathway analysis based on combined evidence; # use nB=2000 or more for more accurate results
res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
#make the output fit this screen
res$Name=substr(res$Name,1,10)
#show first 15 pathways, omit KEGG links
res[1:20,-12]


###################################################
### code chunk number 5: SPIA.Rnw:156-159
###################################################
plotP(res,threshold=0.05)
points(I(-log(pPERT))~I(-log(pNDE)),data=res[res$ID=="05210",],col="green",pch=19,cex=1.5)



###################################################
### code chunk number 6: SPIA.Rnw:186-193
###################################################
res$pG=combfunc(res$pNDE,res$pPERT,combine="norminv")
res$pGFdr=p.adjust(res$pG,"fdr")
res$pGFWER=p.adjust(res$pG,"bonferroni")

plotP(res,threshold=0.05)
points(I(-log(pPERT))~I(-log(pNDE)),data=res[res$ID=="05210",],col="green",pch=19,cex=1.5)



###################################################
### code chunk number 7: SPIA.Rnw:209-216
###################################################
data(Vessels)
# pathway analysis based on combined evidence; # use nB=2000 or more for more accurate results
res<-spia(de=DE_Vessels,all=ALL_Vessels,organism="hsa",nB=500,plots=FALSE,beta=NULL,verbose=FALSE)
#make the output fit this screen
res$Name=substr(res$Name,1,10)
#show first 15 pathways, omit KEGG links
res[1:15,-12]


###################################################
### code chunk number 8: SPIA.Rnw:222-223
###################################################
res[,"KEGGLINK"][20]


###################################################
### code chunk number 9: SPIA.Rnw:235-246
###################################################
  rel<-c("activation","compound","binding/association","expression","inhibition",
"activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
"inhibition_dephosphorylation","dissociation","dephosphorylation",
"activation_dephosphorylation","state change","activation_indirect effect",
"inhibition_ubiquination","ubiquination", "expression_indirect effect",
"inhibition_indirect effect","repression","dissociation_phosphorylation",
"indirect effect_phosphorylation","activation_binding/association",
"indirect effect","activation_compound","activation_ubiquination")
beta=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
names(beta)<-rel
cbind(beta)


###################################################
### code chunk number 10: SPIA.Rnw:255-258
###################################################
load(file=paste(system.file("extdata/hsaSPIA.RData",package="SPIA")))
names(path.info[["05210"]])
path.info[["05210"]][["activation"]][25:35,30:40]


###################################################
### code chunk number 11: SPIA.Rnw:267-286
###################################################
library(graph)
library(Rgraphviz)

plotG<-function(B){
 nnms<-NULL;colls<-NULL
 mynodes<-colnames(B)
 L<-list();
 n<-dim(B)[1]
 for (i in 1:n){
 L[i]<-list(edges=rownames(B)[abs(B[,i])>0])
 if(sum(B[,i]!=0)>0){
 nnms<-c(nnms,paste(colnames(B)[i],rownames(B)[B[,i]!=0],sep="~"))
 }
 }
 names(L)<-rownames(B)
 g<-new("graphNEL",nodes=mynodes,edgeL=L,edgemode="directed")
 plot(g)
}



###################################################
### code chunk number 12: SPIA.Rnw:295-296
###################################################
plotG(path.info[["04012"]][["activation"]])


###################################################
### code chunk number 13: SPIA.Rnw:312-317
###################################################
mydir=system.file("extdata/keggxml/hsa",package="SPIA")
dir(mydir)
makeSPIAdata(kgml.path=mydir,organism="hsa",out.path="./")
res<-spia(de=DE_Colorectal, all=ALL_Colorectal, organism="hsa",data.dir="./")
res[,-12]


