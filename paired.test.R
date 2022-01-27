library(Seurat)
library(tidyverse)
library(ggplot2)
library(qs)
library(Matrix)
###Load qs data ###
mBC<-qread("DSF/CCA_DSF_20211220.qs")

#Extraction by orig.ident(exclude DSF041,DSF051,DSF081,DSF123)#
mBC@meta.data$CB=rownames(mBC@meta.data)
keep.CB2=mBC@meta.data[!mBC@meta.data$orig.ident %in% c("DSF041","DSF051","DSF081","DSF123"),c("CB","celltype")]

#Extraction by macrophage, neutrophil , B/T cells and DC populations#
keep.celltype=c("Mac_C1Qhi_SELENOPhi","Mac_int","Mac_AREG1hi","Mac_CD163hi","Mac_SPP1hi","Macrophage_proliferated")
keep.celltype="Neutrophils"
keep.celltype=c("Tcell_CD8T_3","Tcell_CD4T","Tcell_CD8T_0","Tcell_CD8T_1","Tcell_CD8T_2","Tcell_CD8T_4","Tcell_Treg")
keep.celltype="Bcell"
keep.celltype=c("DC","cDC3","pDC")
############################################################
keep.CB2=keep.CB2[keep.CB2$celltype %in% keep.celltype,"CB"]
Paired.test=mBC[["RNA"]]@counts[,keep.CB2]
Paired.test <- CreateSeuratObject(counts = Paired.test, project = "Paired.test")
Paired.test <- NormalizeData(Paired.test)
sub.set=mBC@meta.data[,c("seurat_clusters","celltype","treatment","CB")]
Paired.test@meta.data$CB=rownames(Paired.test@meta.data)
Paired.test@meta.data=inner_join(Paired.test@meta.data,sub.set,by="CB")
rownames(Paired.test@meta.data)=Paired.test@meta.data$CB

#Extract normalized gene-expression table from DSF011 to DSF132 and calculate the mean value of each gene in sepecific cell Paired.testset#
Paire_table = Paired.test[["RNA"]]@data
idents = unique(Paired.test@meta.data$orig.ident)

#Calculate means for the expression of each gene in sepecific cell Paired.testsets#
res=data.frame()

for (i in c(1:length(idents))){
  
  tmp1=Paire_table[,colnames(Paire_table)%in% rownames(Paired.test@meta.data[Paired.test@meta.data$orig.ident %in% idents[i],])]
  tmp1 = Matrix(tmp1, sparse=FALSE)
  tmp2 = rowMeans(tmp1)
  tmp2 = as.data.frame(tmp2)
  rownames(tmp2)=rownames(tmp1)
  colnames(tmp2)=idents[i]
  if(i==1) {
    res=tmp2
  }else{
    res = cbind(res, tmp2) 
  }
 
}
write.table(res,file="cell.population_paired.test.txt",quote=F,sep="\t",row.names=T,col.names=T)


