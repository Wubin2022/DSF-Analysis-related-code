library(Seurat)
library(tidyverse)
library(ggplot2)
library(qs)
library(xlsx)
###Load qs data ###
mBC<-qread("DSF/CCA_DSF_20211220.qs")

###Find out DEGs in each cell population with non-paired test###
#Extraction by orig.ident(exclude DSF041,DSF051,DSF081,DSF123)#
mBC@meta.data$CB=rownames(mBC@meta.data)
keep.CB=mBC@meta.data[!mBC@meta.data$orig.ident %in% c("DSF041","DSF051","DSF081","DSF123"),c("CB","celltype")]

#Extraction by macrophage, neutrophil, T/B cells,DC#
keep.celltype=c("Mac_C1Qhi_SELENOPhi","Mac_int","Mac_AREG1hi","Mac_CD163hi","Mac_SPP1hi","Macrophage_proliferated")
keep.celltype="Neutrophils"
keep.celltype=c("Tcell_CD8T_3","Tcell_CD4T","Tcell_CD8T_0","Tcell_CD8T_1","Tcell_CD8T_2","Tcell_CD8T_4","Tcell_Treg")
keep.celltype="Bcell"
keep.celltype=c("DC","cDC3","pDC")
########################################################
keep.CB=keep.CB[keep.CB$celltype %in% keep.celltype,"CB"]
NP.test=mBC[["RNA"]]@counts[,keep.CB]
NP.test <- CreateSeuratObject(counts = NP.test, project = "NP.test")
NP.test <- NormalizeData(NP.test)
sub.set=mBC@meta.data[,c("seurat_clusters","celltype","treatment","CB")]
NP.test@meta.data$CB=rownames(NP.test@meta.data)
NP.test@meta.data=inner_join(NP.test@meta.data,sub.set,by="CB")
rownames(NP.test@meta.data)=NP.test@meta.data$CB

#DEG analysis with Find Markers#
Idents(NP.test)="treatment"
markers <- FindMarkers(NP.test, ident.1="post-treatment",ident.2="pre-treatment",logfc.threshold = 0.25, min.pct = 0.1,
                       only.pos = FALSE, test.use = "wilcox",return.thresh = 0.05)
markers$gene=rownames(markers)
markers$cluster=ifelse(markers$avg_log2FC>0,"post-treatment_up","post-treatment_down")
markers=markers%>%dplyr::filter( p_val_adj < 0.05 )
markers_df = markers
markers_df = markers_df %>% arrange(cluster,desc(avg_log2FC))
write.xlsx(markers_df,file =paste0("DC.","DEGs.treatment_log2fc0.25_padj0.05_minpct0.1.xlsx"),row.names = F,col.names = T)

###Visualization of GO Term###
GO_table.data=read_excel("GO_Term.data.xlsx")
GO_table.data$log10_p_neg=abs(GO_table.data$log10_p_neg)
GO_table.data=GO_table.data%>%arrange(log10_p_neg)
GO_table.data$Description=factor(GO_table.data$Description,levels = GO_table.data$Description)
ggplot(GO_table.data,aes(x=Description,y=log10_p_neg,fill=log10_p_neg))+
  geom_bar(stat="identity",width = 0.8,color="black")+theme(text = element_text(size = 20))+
  coord_flip()
ggsave("Cell.population_GO_ana.png",width = 40, height = 30, units = "cm")

###FeaturePlot and Visualization###
FeaturePlot(object = mBC, features = "MKI67", cols = c("grey", "red"),
            reduction = "umap", min.cutoff = "q10",pt.size = 0.01)

ggsave("MKI67.png",width = 10,height = 10,units = "cm")


