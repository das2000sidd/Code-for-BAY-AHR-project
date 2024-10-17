setwd("~/Desktop/PhD_Project_related/COMPARING_WT_LT_BAY_AHRKO_DMSO_SAMP")

ahrko=read.table(file="Comparison_WT_vs_AHRKO_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")
ltbay=read.table(file="Comparison_WT_vs_LTBAY_both_DMSO.txt",header = T,stringsAsFactors = F,sep="\t")



ahrko_sig=subset(ahrko,ahrko$padj < 0.01 & abs(ahrko$log2FoldChange) > 1)
ltbay_sig=subset(ltbay,ltbay$padj < 0.01 & abs(ltbay$log2FoldChange) > 1)

ahrko_up=subset(ahrko_sig,ahrko_sig$log2FoldChange > 0) ## 701
ahrko_dn=subset(ahrko_sig,ahrko_sig$log2FoldChange < 0) ## 1159

ltbay_up=subset(ltbay_sig,ltbay_sig$log2FoldChange > 0) ## 895
ltbay_dn=subset(ltbay_sig,ltbay_sig$log2FoldChange < 0) ## 1187


nrow(ahrko_up) + nrow(ahrko_dn) ## 1860
nrow(ltbay_up) + nrow(ltbay_dn) ## 2082

length(intersect(ahrko_up$Ensembl,ltbay_up$Ensembl)) ## 295
length(intersect(ahrko_dn$Ensembl,ltbay_dn$Ensembl)) ## 419


ahrko_up=ahrko_up[order(- ahrko_up$log2FoldChange),]
ahrko_dn=ahrko_dn[order(ahrko_dn$log2FoldChange),]

ltbay_up=ltbay_up[order(- ltbay_up$log2FoldChange),]
ltbay_dn=ltbay_dn[order(ltbay_dn$log2FoldChange),]




library(org.Mm.eg.db)
ahrko$Entrez <- mapIds(org.Mm.eg.db, ahrko$Ensembl,keytype="ENSEMBL", column="ENTREZID")
ltbay$Entrez <- mapIds(org.Mm.eg.db, ltbay$Ensembl,keytype="ENSEMBL", column="ENTREZID")

ahrko$Symbol <- mapIds(org.Mm.eg.db, ahrko$Entrez,keytype="ENTREZID", column="SYMBOL")
ltbay$Symbol <- mapIds(org.Mm.eg.db, ltbay$Entrez,keytype="ENTREZID", column="SYMBOL")

ltbay$GeneName <- mapIds(org.Mm.eg.db, ltbay$Entrez,keytype="ENTREZID", column="GENENAME")
ahrko$GeneName <- mapIds(org.Mm.eg.db, ahrko$Entrez,keytype="ENTREZID", column="GENENAME")


ahrko$Symbol=as.character(ahrko$Symbol)
ltbay$Symbol=as.character(ltbay$Symbol)
ltbay$GeneName=as.character(ltbay$GeneName)
ahrko$GeneName=as.character(ahrko$GeneName)


ahrko_no_null=subset(ahrko,ahrko$Symbol!="NULL")
ahrko_no_null=ahrko_no_null[order(-ahrko_no_null$log2FoldChange),]

ahrko_noGm_riken=ahrko_no_null[- grep("RIKEN",ahrko_no_null$GeneName),]
#ahrko_noGm_riken=ahrko[- grep("predicted",ahrko$Gene_Name),]
ahrko_noGm_riken=ahrko_noGm_riken[complete.cases(ahrko_noGm_riken),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("predicted",ahrko_noGm_riken$GeneName),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("Riken",ahrko_noGm_riken$GeneName),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("pseudogene",ahrko_noGm_riken$GeneName),]


up_ahrko_genes=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange > 1 & ahrko_noGm_riken$padj < 0.01)
down_ahrko_genes=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange < -1 & ahrko_noGm_riken$padj < 0.01)


down_ahrko_genes=down_ahrko_genes[order(down_ahrko_genes$log2FoldChange),]


ahrko_top_10_up=up_ahrko_genes[1:10,]
ahrko_top_10_dn=down_ahrko_genes[1:10,]

ahrko_top_10_up_dn=rbind(ahrko_top_10_up,ahrko_top_10_dn)


ltbay_no_null=subset(ltbay,ltbay$Symbol!="NULL")
ltbay_no_null=ltbay_no_null[order(-ltbay_no_null$log2FoldChange),]

ltbay_noGm_riken=ltbay_no_null[- grep("RIKEN",ltbay_no_null$GeneName),]
#ltbay_noGm_riken=ltbay[- grep("predicted",ltbay$Gene_Name),]
ltbay_noGm_riken=ltbay_noGm_riken[complete.cases(ltbay_noGm_riken),]
ltbay_noGm_riken=ltbay_noGm_riken[- grep("predicted",ltbay_noGm_riken$GeneName),]
ltbay_noGm_riken=ltbay_noGm_riken[- grep("Riken",ltbay_noGm_riken$GeneName),]
ltbay_noGm_riken=ltbay_noGm_riken[- grep("pseudogene",ltbay_noGm_riken$GeneName),]

up_ltbay_genes=subset(ltbay_noGm_riken,ltbay_noGm_riken$log2FoldChange > 1 & ltbay_noGm_riken$padj < 0.01)
down_ltbay_genes=subset(ltbay_noGm_riken,ltbay_noGm_riken$log2FoldChange < -1 & ltbay_noGm_riken$padj < 0.01)

library(VennDiagram)

p1=venn.diagram(list(`AHRKO vs WT`=up_ahrko_genes$Ensembl,`LT-BAY vs WT`=up_ltbay_genes$Ensembl),filename=NULL,cex=2,fill=c("white","white"),main.cex=c(1),cat.cex=c(2,2),cat.dist=c(.02,.02),cat.pos=c(-170,170),main="")
grid.draw(p1)

p2=venn.diagram(list(`AHRKO vs WT`=down_ahrko_genes$Ensembl,`LT-BAY vs WT`=down_ltbay_genes$Ensembl),filename=NULL,cex=2,fill=c("white","white"),main.cex=c(1),cat.cex=c(2,2),cat.dist=c(.02,.02),cat.pos=c(-170,170),main="")
grid.draw(p2)

tiff(file="Overlap_genes_up_AHRKO_and_LTBAY_no_Gm_riken_predicted_and_pseudogenes.tiff",res=300,height = 2000,width = 2000)
grid.draw(p1)
dev.off()

tiff(file="Overlap_genes_down_AHRKO_and_LTBAY_no_Gm_riken_predicted_and_pseudogenes.tiff",res=300,height = 2000,width = 2000)
grid.draw(p2)
dev.off()

overlapping_up=intersect(up_ahrko_genes$Ensembl,up_ltbay_genes$Ensembl)
overlapping_dn=intersect(down_ahrko_genes$Ensembl,down_ltbay_genes$Ensembl)

up_ahrko_overlapping=subset(up_ahrko_genes,up_ahrko_genes$Ensembl %in% overlapping_up)
dn_ahrko_overlapping=subset(down_ahrko_genes,down_ahrko_genes$Ensembl %in% overlapping_dn)

`%ni%`=Negate(`%in%`)
up_ahrko_unique=subset(up_ahrko_genes,up_ahrko_genes$Ensembl %ni% overlapping_up)
dn_ahrko_unique=subset(down_ahrko_genes,down_ahrko_genes$Ensembl %ni% overlapping_dn)

up_dn_ahrko_unique=rbind(up_ahrko_unique,dn_ahrko_unique)
up_dn_ahrko_overlapping=rbind(up_ahrko_overlapping,dn_ahrko_overlapping)


up_ltbay_unique=subset(up_ltbay_genes,up_ltbay_genes$Ensembl %ni% overlapping_up)
dn_ltbay_unique=subset(down_ltbay_genes,down_ltbay_genes$Ensembl %ni% overlapping_dn)

up_dn_ltbay_unique=rbind(up_ltbay_unique,dn_ltbay_unique)

up_ltbay_overlapping=subset(up_ltbay_genes,up_ltbay_genes$Ensembl %in% overlapping_up)
dn_ltbay_overlapping=subset(down_ltbay_genes,down_ltbay_genes$Ensembl %in% overlapping_dn)

up_dn_ltbay_overlapping=rbind(up_ltbay_overlapping,dn_ltbay_overlapping)


write.table(up_dn_ahrko_unique,file="AHRKO_unique_no_predicted_Riken_pseudogene.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(up_dn_ltbay_unique,file="LTBAY_unique_no_predicted_Riken_pseudogene.txt",col.names = T,row.names = F,sep="\t",quote = F)

write.table(up_dn_ahrko_overlapping,file="AHRKO_overlapping_genes_no_predicted_Riken_pseudogene.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(up_dn_ltbay_overlapping,file="LTBAY_overlapping_genes_no_predicted_Riken_pseudogene.txt",col.names = T,row.names = F,sep="\t",quote = F)



library(eulerr)

up.overlap <- euler(c("LT-BAY vs WT" = 518, "AHRKO vs WT" = 392,"LT-BAY vs WT&AHRKO vs WT"=285))
plot(up.overlap, counts = TRUE, font=1, cex=4, alpha=0.5,
     fill=c("white", "white","white"),main="",labels=list(fontsize=15,cex=1),quantities=list(cex=3))




down_ltbay_genes=down_ltbay_genes[order(down_ltbay_genes$log2FoldChange),]


ltbay_top_10_up=up_ltbay_genes[1:10,]
ltbay_top_10_dn=down_ltbay_genes[1:10,]


ltbay_top_10_up_dn=rbind(ltbay_top_10_up,ltbay_top_10_dn)



# ltbay_top_10_up_dn, ahrko_top_10_up_dn

library(ggplot2)
library(forcats)

## LTBAY genes to plot




ltbay_top_10_up_dn$`Neg log10(adj. p-val)`=-log10(ltbay_top_10_up_dn$padj)
ltbay_top_10_dn$`Neg log10(adj. p-val)`=-log10(ltbay_top_10_dn$padj)
ltbay_top_10_up$`Neg log10(adj. p-val)`=-log10(ltbay_top_10_up$padj)

ltbay_top_10_up_dn$Direction=ifelse(ltbay_top_10_up_dn$log2FoldChange > 0,"Up","Down")
ltbay_top_10_dn$Direction="Down"
ltbay_top_10_up$Direction="Up"


intersect(ltbay_top_10_up_dn$Ensembl,ltbay_top_10_up$Ensembl)
intersect(ltbay_top_10_up_dn$Ensembl,ltbay_top_10_dn$Ensembl)


plot_to_save=ggplot(ltbay_top_10_up_dn, aes(y = reorder(Symbol, log2FoldChange), x = log2FoldChange,size=`Neg log10(adj. p-val)`,color=Direction)) + 
  geom_point(data = ltbay_top_10_dn, aes(y = reorder(Symbol, log2FoldChange), x = log2FoldChange,size=`Neg log10(adj. p-val)`)) +
  geom_point(data = ltbay_top_10_up, aes(y = reorder(Symbol, log2FoldChange), x = log2FoldChange,size=`Neg log10(adj. p-val)`)) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue"))+
  xlim(-13,25) + geom_vline(xintercept = 0,linetype = "longdash") + ylab("Gene ID")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("LTBAY vs WT") + theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=5),axis.title = element_text(size=7))+guides(colour=FALSE)+theme(plot.title = element_text(size=5))+theme(legend.title = element_text(size=5),legend.text = element_text(size=5))


tiff(file="LTBAY_vsWT_top_10_genes_logFC_plot.tiff",res=300,height = 1500,width = 1500)
grid.draw(plot_to_save)
dev.off()


write.table(ltbay_top_10_up_dn,file="Top_10_up_down_genes_by_logFC_LT-BAY_vsWT.txt",col.names = T,row.names = F,sep="\t",quote = F)




