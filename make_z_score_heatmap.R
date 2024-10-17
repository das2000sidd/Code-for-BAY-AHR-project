setwd("~/Desktop/PhD_Project_related/COMPARING_WT_LT_BAY_AHRKO_DMSO_SAMP")

ahrko=read.table(file="Comparison_WT_vs_AHRKO_both_DMSO.txt",header = T,sep="\t",stringsAsFactors = F)
ltbay=read.table(file="Comparison_WT_vs_LTBAY_both_DMSO.txt",header = T,sep="\t",stringsAsFactors = F)

ahrko=ahrko[complete.cases(ahrko),]
ltbay=ltbay[complete.cases(ltbay),]

ahrko_sig=subset(ahrko,abs(ahrko$log2FoldChange) > 1 & ahrko$padj < 0.01)
ahrko_sig_up=subset(ahrko_sig,ahrko_sig$log2FoldChange > 0)
ahrko_sig_dn=subset(ahrko_sig,ahrko_sig$log2FoldChange < 0)


ltbay_sig=subset(ltbay,abs(ltbay$log2FoldChange) > 1 & ltbay$padj < 0.01)
ltbay_sig_up=subset(ltbay_sig,ltbay_sig$log2FoldChange > 0)
ltbay_sig_dn=subset(ltbay_sig,ltbay_sig$log2FoldChange < 0)


pdf('Overlap_genes_up_AHRKO_and_LTBAY_relative_to_WT.pdf')
p1=venn.diagram(list(`AHRKO vs WT`=ahrko_sig_up$Ensembl,`LT-BAY vs WT`=ltbay_sig_up$Ensembl),filename=NULL,cex=2,fill=c("white","white"),main.cex=c(1),cat.cex=c(1.5,1.5),cat.dist=c(.05,.05),cat.pos=c(-150,150),main="Overlap of genes significantly up in LT-BAY vs WT and AHRKO vs WT for DMSO treatment")
grid.draw(p1)
dev.off()


ahrko_and_ltbay_changed=intersect(ahrko_sig$Ensembl,ltbay_sig$Ensembl)
ahrko_and_ltbay_up=intersect(ahrko_sig_up$Ensembl,ltbay_sig_up$Ensembl)
ahrko_and_ltbay_dn=intersect(ahrko_sig_dn$Ensembl,ltbay_sig_dn$Ensembl)


  
cpm=read.table(file="Normalised_CPM_count.txt",header = T,sep="\t",stringsAsFactors = F)

wt_samp=paste("PW1_",1:4,sep="")
ko_samp=paste("PK",1:4,sep="")
ko_samp=paste(ko_samp,"_1",sep="")
ltbay_samp=paste("LB",1:4,sep="")
ltbay_samp=paste(ltbay_samp,"_1",sep="")



cpm_plot=cpm[c(ahrko_and_ltbay_up,ahrko_and_ltbay_dn),]

colnames(cpm_plot)[1:4]=c(paste("LTBay Rep",1:4,sep=" "))
colnames(cpm_plot)[5:8]=c(paste("AHRKO Rep",1:4,sep=" "))
colnames(cpm_plot)[9:12]=c(paste("WT Rep",1:4,sep=" "))

# labels_col = c("WT Rep 1","WT Rep 2","WT Rep 3","WT Rep 4","AHRKO Rep1","AHRKO Rep2","AHRKO Rep3","AHRKO Rep4","LT Bay Rep1","LT Bay Rep2","LT Bay Rep3","LT Bay Rep4")
library(pheatmap)

tiff(filename="Z score heatmap of DEG in both AHRKO and LT BAY relative to WT 2.tiff",width = 2500,height = 2000,res=500)
pheatmap(cpm_plot,scale = "row",show_rownames = FALSE,main="",cluster_rows = FALSE,fontsize_col=10,labels_col = c("LT Bay","LT Bay","LT Bay","LT Bay","AHRKO","AHRKO","AHRKO","AHRKO","WT","WT","WT","WT"))
dev.off()


