
library(tximportData)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(biomaRt)
library(tibble)
library(tximport)

res=read.table(file="GRCm39_gene_and_transcript_stable_ID_version.txt",header = T,sep="\t",stringsAsFactors = F)
res=res[,c(1,4)]
colnames(res)=c("GENEID","TXNAME")
res=res[,c(2,1)]
res_tibble=as_tibble(res)


files=list.files(path="/Users/siddhaduio.no/Desktop/PhD_Project_related/COMPARING_WT_LT_BAY_AHRKO_DMSO_SAMP",full.names = T)

files_wt=files[grep("PW",files)] ## WT
files_ahrko=files[grep("PK",files)] ## AHRKO
files_ltbay=files[grep("LB",files)] ## LT-BAY

all_files=c(files_ltbay,files_ahrko,files_wt)

files=as.data.frame(files)

files_quant=paste(files$files,"/quant.sf",sep="")

sample_order=c("LB1_1","LB2_1","LB3_1","LB4_1","PK1_1","PK2_1","PK3_1","PK4_1","PW1_1","PW2_1","PW3_1","PW4_1")


txi = tximport(files_quant,type = "salmon",tx2gene = res_tibble)
names(txi)
head(txi$counts)
colnames(txi$counts)=sample_order
head(txi$counts)

library(DESeq2)
library(edgeR)

genotype_info=c(rep("LB DMSO",4),rep("AHRKO DMSO",4),rep("WT DMSO",4))
genotype_info=as.factor(genotype_info)
genotype_info=as.data.frame(genotype_info)
model_mat_geno_treat=model.matrix( ~ 0 + genotype_info,data=genotype_info)


write.table(txi$counts,file="Salmon_Counts_Object.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(txi$abundance,file="Salmon_Abundance_Object.txt",col.names = T,row.names = T,sep="\t",quote = F)

library(DESeq2)
library(edgeR)
dds <- DESeqDataSetFromTximport(txi, genotype_info, model_mat_geno_treat)

idx <- rowSums( cpm(dds) > 2) >= 2 ## This is not necessary

#head(all_counts_combined[,c(6:8)])

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")


## Exploratory analysis and visualisation
nrow(dds)

# *** variance stabilizing transformation and the rlog***

## Showing effect on some simulated data

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)


vsd <- vsd(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)


library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)


sampleDists <- dist(t(assay(vsd)))
sampleDists


library("pheatmap")
library("RColorBrewer")


sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- vsd$genotype_info
colnames(sampleDistMatrix) <- vsd$genotype_info
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)



library(ggplot2)
###PCA plot***
tiff(file="PCA plot with 3 groups.tiff",res=300,height = 2000,width = 2000)
plotPCA(vsd, intgroup = c("genotype_info")) +ggtitle("")+theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=12))
dev.off()


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype_info")])
top_20_genes=assay(ntd)[select,]
colnames(df)="Condition"
rownames(df)=colnames(top_20_genes)
pheatmap(top_20_genes, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


library(Glimma)

dds$group <- factor(paste0(dds$genotype_info))
design(dds) = ~ group
dds <- DESeq(dds)
plotMDS(dds)
glimmaMDS(dds)


contrast_bay_dmso <- c("group", "WTBay", "WTDMSO")

#res_Bay_DMSO <- results(dds, contrast=contrast_bay_dmso,pAdjustMethod = "BH",format = "DataFrame")
#mcols(res_Bay, use.names = TRUE)
#res_Bay_DMSO_df=as.data.frame(res_Bay_DMSO)
#res_Bay_DMSO_df$Ensembl=rownames(res_Bay_DMSO_df)
#summary(res_Bay_DMSO)
#res_Bay_DMSO_df=res_Bay_DMSO_df[complete.cases(res_Bay_DMSO_df),]


contrast_kyn_dmso <- c("group", "WTKyn", "WTDMSO")

#res_Kyn_DMSO <- results(dds, contrast=contrast_kyn_dmso,pAdjustMethod = "BH",format = "DataFrame")
#res_Kyn_DMSO_df=as.data.frame(res_Kyn_DMSO)
#res_Kyn_DMSO_df$Ensembl=rownames(res_Kyn_DMSO_df)
#summary(res_Bay_DMSO)
#res_Kyn_DMSO_df=res_Kyn_DMSO_df[complete.cases(res_Kyn_DMSO_df),]


contrast_wt_vs_ahrko <- c("group", "AHRKo DMSO", "WT DMSO")

res_wt_vs_ahrko <- results(dds, contrast=contrast_wt_vs_ahrko,pAdjustMethod = "BH",format = "DataFrame")
res_wt_vs_ahrko_df=as.data.frame(res_wt_vs_ahrko)
res_wt_vs_ahrko_df$Ensembl=rownames(res_wt_vs_ahrko_df)
summary(res_wt_vs_ahrko)
res_wt_vs_ahrko_df=res_wt_vs_ahrko_df[complete.cases(res_wt_vs_ahrko_df),]


contrast_wt_vs_ltbay_dmso <- c("group", "LB DMSO", "WT DMSO")

res_wt_vs_ltbay <- results(dds, contrast=contrast_wt_vs_ltbay,pAdjustMethod = "BH",format = "DataFrame")
res_wt_vs_ltbay_df=as.data.frame(res_wt_vs_ltbay)
res_wt_vs_ltbay_df$Ensembl=rownames(res_wt_vs_ltbay_df)
summary(res_wt_vs_ltbay)
res_wt_vs_ltbay_df=res_wt_vs_ltbay_df[complete.cases(res_wt_vs_ltbay_df),]

View(res_wt_vs_ahrko_df)
View(res_wt_vs_ltbay_df)


#write.table(res_wt_vs_ahrko_df,file="Comparison_WT_vs_AHRKO_both_DMSO.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(res_wt_vs_ltbay_df,file="Comparison_WT_vs_LTBAY_both_DMSO.txt",col.names = T,row.names = T,sep="\t",quote = F)




library(edgeR)
cpm.nor.count=cpm(dds,normalized.lib.sizes = TRUE,log=FALSE)

write.table(cpm.nor.count,file="Normalised_CPM_count.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(rownames(dds),file="Background_genes.txt",col.names = F,row.names = F,sep="\t",quote = F)



