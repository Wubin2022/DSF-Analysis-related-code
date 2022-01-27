library(tidyverse)

pseudobulk.df=read.table("_paired.test.txt",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
pseudobulk.df=pseudobulk.df[rowSums(pseudobulk.df) > 0,]

p.value.vector=c()
log2FC.postVSpre.vector=c()

for (i in rownames(pseudobulk.df)) {
  
  each.gene.df=pseudobulk.df[i,]
  index.pre=which(str_detect(colnames(each.gene.df),"1$"))
  index.post=which(str_detect(colnames(each.gene.df),"2$"))
  pre.v=as.numeric(each.gene.df[,index.pre])
  post.v=as.numeric(each.gene.df[,index.post])
  
  test.res=wilcox.test(pre.v, post.v, paired = TRUE, alternative = "two.sided")
  
  p.value.vector=append(p.value.vector,test.res$p.value)
  tmp=log2(mean(post.v)+0.0001) - log2(mean(pre.v)+0.0001)
  log2FC.postVSpre.vector=append(log2FC.postVSpre.vector,tmp)
  
}
###bonferroni or BH correction
pseudobulk.df$p.value=p.value.vector
pseudobulk.df$log2FC.postVSpre=log2FC.postVSpre.vector
pseudobulk.df$p.value.adj=p.adjust(pseudobulk.df$p.value,method = "bonferroni")
pseudobulk.df$gene=rownames(pseudobulk.df)

write.table(pseudobulk.df[,c("gene","log2FC.postVSpre","p.value","p.value.adj")],file = "_paired.test_bonferroni_adj.txt",quote = F,sep = "\t",row.names = F,col.names = T)
