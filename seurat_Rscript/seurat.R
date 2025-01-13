library(SingleCellExperiment)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(clustree)
library(stringr)
library(RColorBrewer)
library(slingshot)
library(tradeSeq)


########### load data ###########

mock.data<-Read10X("../1.STARsolo/mock_Solo.out/GeneFull/filtered")
C26.data<-Read10X("../1.STARsolo/C26_Solo.out/GeneFull/filtered")

mock<-CreateSeuratObject(counts=mock.data,project='Mock',min.cells = 1,min.features = 10)
C26<-CreateSeuratObject(counts=C26.data,project='C26',min.cells = 1,min.features = 10)

Neiss<-merge(x=mock,y=C26)


########### QC, not setting any threshold ###########

Neiss$log10GenesPerUMI<-log10(Neiss$nFeature_RNA)/log10(Neiss$nCount_RNA)

# Create metadata dataframe
metadata <- Neiss@meta.data
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$orig.ident<-factor(metadata$orig.ident,levels=c('Mock','C26'))

#update metadata in the Seurat object
Neiss@meta.data<-metadata

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

summary(metadata)
#Mock:3542     
#C26 :6119

VlnPlot(Neiss, features = c("nGene", "nUMI", "percent.mt"),group.by = 'orig.ident', ncol = 3,raster=F)

saveRDS(Neiss, file="Neiss.RDS")

########### normalization and integration ###########

Neiss <- NormalizeData(Neiss)
Neiss <- FindVariableFeatures(Neiss)
Neiss <- ScaleData(Neiss)
Neiss <- RunPCA(Neiss)

DimPlot(Neiss, reduction = "pca") + NoLegend()
DimHeatmap(Neiss, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(Neiss,ndims=30)

Neiss <- FindNeighbors(Neiss, dims = 1:5, reduction = "pca")
Neiss <- FindClusters(Neiss, resolution = 0.4, cluster.name = "unintegrated_clusters")
Neiss <- RunUMAP(Neiss, dims = 1:5, reduction = "pca", reduction.name = "umap.unintegrated")

#cca integration
Neiss <- IntegrateLayers(
  object= Neiss, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = T
)

saveRDS(Neiss, file="Neiss.RDS")

Neiss <- FindNeighbors(Neiss, reduction = "integrated.cca", dims = 1:5)
#Neiss <- FindClusters(Neiss, resolution = 0.55, cluster.name = "cca_clusters")

for(res in seq(0.1,1.5,0.1)){
  Neiss<-FindClusters(Neiss,
                                 resolution = res,
                                 cluster.name = paste0("clusters_res",res))
}

clustree<-clustree(Neiss,prefix = "clusters_res") 
clustree

Neiss <- FindClusters(Neiss, resolution = 0.3, cluster.name = "used_clusters_res0.3")
Neiss <- RunUMAP(Neiss, reduction = "integrated.cca", dims = 1:5, reduction.name = "umap.cca")

DimPlot(Neiss,reduction='umap.cca',split.by = 'orig.ident',pt.size = 1)
DimPlot(Neiss,reduction='umap.cca',group.by = 'orig.ident',pt.size = 1)


#SCT integration
Neiss <- SCTransform(Neiss)
Neiss <- RunPCA(Neiss, verbose = F)
Neiss <- IntegrateLayers(
  object = Neiss,
  method = CCAIntegration,
  normalization.method = "SCT",
  verbose = T,new.reduction='integrated.sct'
)

saveRDS(Neiss, file="Neiss.RDS")

Neiss <- FindNeighbors(Neiss, dims = 1:5, reduction = "integrated.sct")
#Neiss <- FindClusters(Neiss, resolution = 0.4,cluster.name = 'sct_clusters')


for(res in seq(0.1,1.5,0.1)){
  Neiss<-FindClusters(Neiss,
                      resolution = res,
                      cluster.name = paste0("clusters_sct",res))
}

clustree2<-clustree(Neiss,prefix = "clusters_sct") 
clustree2

Neiss <- FindClusters(Neiss, resolution = 0.4, cluster.name = "used_sct.clusters_res0.4") # use 0.4
Neiss <- RunUMAP(Neiss, reduction = "integrated.sct", dims = 1:5, reduction.name = "umap.sct.cca")

DimPlot(Neiss,reduction='umap.sct.cca',split.by = 'orig.ident',pt.size = 1)
DimPlot(Neiss,reduction='umap.sct.cca',group.by = 'orig.ident',pt.size = 1)

########### set color ###########
col<-brewer.pal(5,'Set1')
pie(rep(1,5),col = col)

########### look for signatures res=0.4 ###########
Neiss <- PrepSCTFindMarkers(Neiss)

Neiss.marker<- FindAllMarkers(Neiss, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Neiss.marker$symbols<-ID2gene[Neiss.marker$gene]
Neiss.marker <- Neiss.marker %>%
  mutate(symbols=case_when(!is.na(symbols)~symbols,
                           T~gene))
Neiss.marker$symbols<-str_split_i(Neiss.marker$symbols,'\\.',1)
Neiss.marker$ID<-gene2ID[Neiss.marker$symbols]
write.csv(Neiss.marker,'Neiss.marker.csv')


Neiss.marker.full<- FindAllMarkers(Neiss, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
Neiss.marker.full$symbols<-ID2gene[Neiss.marker.full$gene]
Neiss.marker.full <- Neiss.marker.full %>%
  mutate(symbols=case_when(!is.na(symbols)~symbols,
                           T~gene))
Neiss.marker.full$symbols<-str_split_i(Neiss.marker.full$symbols,'\\.',1)
Neiss.marker.full$ID<-gene2ID[Neiss.marker.full$symbols]
write.csv(Neiss.marker.full,'Neiss.marker.full.csv')

Neiss.marker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> top3

Neiss$final.sct.clus<-Neiss$used_sct.clusters_res0.4
colnames(Neiss@meta.data) # rm unused clusters
Neiss@meta.data<-select(Neiss@meta.data,c(orig.ident,nUMI,nGene,log10GenesPerUMI,final.sct.clus))

UMAP.clus<-DimPlot(Neiss,pt.size = 1,split.by = 'orig.ident')+
  scale_color_manual(values=col)

dev.new()
pdf("UMAP.clus.pdf",width=10)
UMAP.clus
dev.off()


hm.marker<-DoHeatmap(subset(Neiss,downsample=100), features = top3$gene)
#feature<-as.character(hm.marker$data$Feature)
#feature2<-ID2gene[feature]
#feature<-as.data.frame(cbind(feature,feature2))
#feature<-feature %>%
#  mutate(feature2= case_when(is.na(feature2)~feature,
#                             T~feature2))


#feature$feature2<-factor(feature$feature2,levels=rev(c(top3$symbols)[-11]))

#all(hm.marker$data$Feature==feature$feature)
#hm.marker$data$Feature<-feature$feature2

dev.new()
pdf("hm.marker.pdf")
hm.marker
dev.off()

dev.new()
pdf("UMAP.unintegrated.pdf",width=10)
DimPlot(Neiss,pt.size = 1,split.by = 'orig.ident',reduction = 'umap.unintegrated')+
  scale_color_manual(values=col)
dev.off()

#bar plot cell abundance
Neiss_split<-SplitObject(Neiss,split.by = "orig.ident")

mock_count<-as.data.frame(summary(Neiss_split$Mock@active.ident))
colnames(mock_count)<-"mock.count"
mock_count$mock.percent<-c(round(mock_count$mock.count*100/sum(mock_count$mock.count),2))

C26_count<-as.data.frame(summary(Neiss_split$C26@active.ident))
colnames(C26_count)<-"C26.count"
C26_count$C26.percent<-c(round(C26_count$C26.count*100/sum(C26_count$C26.count),2))

comb_count<-as.data.frame(cbind(mock_count,C26_count))
comb_count$cluster<-rownames(comb_count)

comb_count.long<-comb_count %>%
  select(c('mock.percent','C26.percent','cluster')) %>%
  tidyr::gather('group','percent',1:2) %>%
  mutate(group=case_when(group=='mock.percent'~'mock',
                         T~'C26'))

comb_count.long$group<-factor(comb_count.long$group,levels=c('mock','C26'))

bar.abundance<-ggplot(comb_count.long,aes(fill=cluster,y=percent,x=group))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',size = 0.1)+
  #geom_text(aes(label = percent), position = position_fill(vjust = 0.5),size = 8)+
  scale_fill_manual(values=c(col))+
  guides(fill=guide_legend(reverse=T))+
  #theme(text = element_text(size = 30),axis.text.x = element_text(angle=45,hjust=1),
  #      panel.background = element_blank(),
  #      panel.border = element_blank(),
  #      axis.line=element_line(linewidth = 1,color="black"))+
  #xlab("Stages")+
  ylab("Abundance (%)")

dev.new()
pdf("bar.abundance.pdf")
bar.abundance
dev.off()

########### reclassify, merge 0 and 2 ###########
Neiss<-RenameIdents(Neiss,'1'='1',
                    '2'='2',
                    '0'='2',
                    '3'='3',
                    '4'='4',
                    '5'='5')
Neiss$final.sct.reclus<-Idents(Neiss)
saveRDS(Neiss,'Neiss.rds')

Neiss.marker<- FindAllMarkers(Neiss, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Neiss.marker$symbols<-ID2gene[Neiss.marker$gene]
Neiss.marker <- Neiss.marker %>%
  mutate(symbols=case_when(!is.na(symbols)~symbols,
                           T~gene))
Neiss.marker$symbols<-str_split_i(Neiss.marker$symbols,'\\.',1)
Neiss.marker$ID<-gene2ID[Neiss.marker$symbols]
write.csv(Neiss.marker,'Neiss.marker.csv')


Neiss.marker.full<- FindAllMarkers(Neiss, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
Neiss.marker.full$symbols<-ID2gene[Neiss.marker.full$gene]
Neiss.marker.full <- Neiss.marker.full %>%
  mutate(symbols=case_when(!is.na(symbols)~symbols,
                           T~gene))
Neiss.marker.full$symbols<-str_split_i(Neiss.marker.full$symbols,'\\.',1)
Neiss.marker.full$ID<-gene2ID[Neiss.marker.full$symbols]
write.csv(Neiss.marker.full,'Neiss.marker.full.csv')

Neiss.marker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> top3

UMAP.clus<-DimPlot(Neiss,pt.size = 1,split.by = 'orig.ident')+
  scale_color_manual(values=col)

dev.new()
pdf("UMAP.clus.pdf",width=10)
UMAP.clus
dev.off()


hm.marker<-DoHeatmap(subset(Neiss,downsample=100), features = top3$gene)

dev.new()
pdf("hm.marker.pdf")
hm.marker
dev.off()

dev.new()
pdf("UMAP.unintegrated.pdf",width=10)
DimPlot(Neiss,pt.size = 1,split.by = 'orig.ident',reduction = 'umap.unintegrated')+
  scale_color_manual(values=col)
dev.off()

########### rename clusters based on function ###########
Neiss<-RenameIdents(Neiss,'1'='Pre-proliferative N.subflava',
                    '2'='Transitional N.subflava',
                    '3'='Metal binding N.subflava',
                    '4'='Proliferative N.subflava',
                    '5'='Peptidoglycan-forming N.subflava')
Neiss$final.sct.annot<-Idents(Neiss)
Neiss$final.sct.annot<-factor(Neiss$final.sct.annot,
                              levels=c('Proliferative N.subflava','Transitional N.subflava',
                                       'Pre-proliferative N.subflava','Peptidoglycan-forming N.subflava','Metal binding N.subflava'))
Idents(Neiss)<-Neiss$final.sct.annot
saveRDS(Neiss,'Neiss.rds')


UMAP.annotated<-DimPlot(Neiss,pt.size = 1,split.by = 'orig.ident')+
  scale_color_manual(values=col)

dev.new()
pdf("UMAP.annotated.pdf",width=10)
UMAP.annotated
dev.off()


#down sample
mock_cells = Cells(Neiss)[which(Neiss$orig.ident == "Mock")]
C26_cells = Cells(Neiss)[which(Neiss$orig.ident == "C26")]
mock_cells = sample(mock_cells, size = 2000)
C26_cells = sample(C26_cells, size = 2000)
Neiss.sub = subset(Neiss, cells = c(mock_cells, C26_cells))


UMAP<-DimPlot(Neiss.sub,pt.size = 1,split.by = 'orig.ident')+
  scale_color_manual(values=col)

dev.new()
pdf("UMAP.annotated.sub.pdf",width=12)
UMAP
dev.off()

hm.marker<-DoHeatmap(subset(Neiss,downsample=100), features = top3$gene)

dev.new()
pdf("hm.marker.annotated.pdf")
hm.marker
dev.off()

dev.new()
pdf("UMAP.unintegrated.pdf",width=10)
DimPlot(Neiss,pt.size = 1,split.by = 'orig.ident',reduction = 'umap.unintegrated')+
  scale_color_manual(values=col)
dev.off()


#bar plot cell abundance
Neiss_split<-SplitObject(Neiss,split.by = "orig.ident")

mock_count<-as.data.frame(summary(Neiss_split$Mock@active.ident))
colnames(mock_count)<-"mock.count"
mock_count$mock.percent<-c(round(mock_count$mock.count*100/sum(mock_count$mock.count),2))

C26_count<-as.data.frame(summary(Neiss_split$C26@active.ident))
colnames(C26_count)<-"C26.count"
C26_count$C26.percent<-c(round(C26_count$C26.count*100/sum(C26_count$C26.count),2))

comb_count<-as.data.frame(cbind(mock_count,C26_count))
comb_count$cluster<-rownames(comb_count)

comb_count.long<-comb_count %>%
  select(c('mock.percent','C26.percent','cluster')) %>%
  tidyr::gather('group','percent',1:2) %>%
  mutate(group=case_when(group=='mock.percent'~'mock',
                         T~'C26'))

comb_count.long$group<-factor(comb_count.long$group,levels=c('mock','C26'))
comb_count.long$cluster<-factor(comb_count.long$cluster,levels=c('Proliferative N.subflava','Transitional N.subflava',
                                                         'Pre-proliferative N.subflava','Peptidoglycan-forming N.subflava','Metal binding N.subflava'))

bar.abundance<-ggplot(comb_count.long,aes(fill=cluster,y=percent,x=group))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',size = 0.1)+
  #geom_text(aes(label = percent), position = position_fill(vjust = 0.5),size = 8)+
  scale_fill_manual(values=c(col))+
  guides(fill=guide_legend(reverse=T))+
  #theme(text = element_text(size = 30),axis.text.x = element_text(angle=45,hjust=1),
  #      panel.background = element_blank(),
  #      panel.border = element_blank(),
  #      axis.line=element_line(linewidth = 1,color="black"))+
  #xlab("Stages")+
  ylab("Abundance (%)")

dev.new()
pdf("bar.abundance.annot.pdf")
bar.abundance
dev.off()

write.csv(comb_count,'comb_count.csv')

########### pseudotime using slingshot ###########
dimred <- Neiss@reductions$umap.sct.cca@cell.embeddings
clustering <- Neiss$final.sct.reclus
counts <- as.data.frame(Neiss@assays$RNA$counts)
counts$var<-rowVars(as.matrix(counts))
counts<-counts %>%
  arrange(desc(var)) %>%
  slice_head(n=100) %>%
  select(-c(var))
counts<-as.matrix(counts)



pal <- col
set.seed(123)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        #end.clus = c('1',"3","4","5"), #define how many branches/lineages to consider
                        start.clus = "4") #define where to start the trajectories
lineages

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}

dev.new()
pdf('pst.neiss.pdf')
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(curves), lwd = 3, col = "black")
dev.off()

# DEG between lineages # plan to use this for different groups
sce <- fitGAM(counts = as.matrix(counts), sds = curves)
plotGeneCount(curves, counts, clusters = clustering, models = sce)

plot_differential_expression <- function(feature_id) {
  #feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, counts, gene = feature_id[1], clusters = clustering, models = sce) + 
                       ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

#Genes that change with pseudotime
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

#Genes that change between two pseudotime points
pseudotime_start_end_association <- startVsEndTest(sce, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)

feature_id <- pseudotime_start_end_association %>% filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)

plot_differential_expression('QFG95-05765')
plotSmoothers(sce, as.matrix(counts), gene = feature_id[1])+
  geom_point(fill=NA,size=0)


#Genes that are different between lineages
different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)

branch_point_association <- earlyDETest(sce)
branch_point_association$feature_id <- rownames(branch_point_association)

feature_id <- branch_point_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)



########### Annotation ###########
background.gene<-as.data.frame(Neiss@assays$RNA$counts)
background.gene$sum<-rowSums(as.matrix(background.gene))
background.gene<-subset(background.gene,!background.gene$sum==0)
background.gene<-rownames(background.gene)
background.gene<-as.data.frame(background.gene)
background.gene$background.gene<-str_split_i(background.gene$background.gene,'\\.',1)

background.gene<-background.gene %>%
  mutate(ID=case_when(grepl(background.gene,pattern='QFG')~background.gene,
                      T~gene2ID[background.gene]))

write.csv(background.gene,'background.gene.csv')

########### genemap ###########
gff <- read.delim("C:/Users/hongs/OneDrive - Nanyang Technological University/Project/LiLiang_neisseria/Mock.gff", header=FALSE)

ID<-str_split_i(gff$V9,pattern=';',1)
gene<-str_split_i(gff$V9,pattern=';',2)
other<-str_split_i(gff$V9,pattern=';',3)
other2<-str_split_i(gff$V9,pattern=';',4)
other3<-str_split_i(gff$V9,pattern=';',5)
other4<-str_split_i(gff$V9,pattern=';',6)
other5<-str_split_i(gff$V9,pattern=';',7)
other6<-str_split_i(gff$V9,pattern=';',8)

df<-as.data.frame(cbind(ID,gene))
df$ID<-str_replace(df$ID,pattern="ID=",'')

#add symbol
df <- df %>%
  mutate(symbol = case_when(grepl(gene,pattern='gene=')~df$gene,
                            T~NA))
df$symbol<-str_replace(df$symbol,pattern="gene=",'')
df$gene<-NULL

# add symbol in other columns
df$other<-other
df <- df %>%
  mutate(symbol = case_when(
    !is.na(df$symbol)~df$symbol,
    grepl(other,pattern='gene=')~df$other,
    grepl(other,pattern='product=')~df$other,
    T~NA))
df$symbol<-str_replace(df$symbol,pattern="product=",'')


df$other<-other3
df <- df %>%
  mutate(symbol = case_when(
    !is.na(df$symbol)~df$symbol,
    grepl(other,pattern='gene=')~df$other,
    grepl(other,pattern='product=')~df$other,
    T~NA))
df$symbol<-str_replace(df$symbol,pattern="product=",'')

df$other<-other4
df <- df %>%
  mutate(symbol = case_when(
    !is.na(df$symbol)~df$symbol,
    grepl(other,pattern='gene=')~df$other,
    grepl(other,pattern='product=')~df$other,
    T~NA))
df$symbol<-str_replace(df$symbol,pattern="product=",'')

df$other<-other5
df <- df %>%
  mutate(symbol = case_when(
    !is.na(df$symbol)~df$symbol,
    grepl(other,pattern='gene=')~df$other,
    grepl(other,pattern='product=')~df$other,
    T~NA))
df$symbol<-str_replace(df$symbol,pattern="product=",'')


df$other<-other6
df <- df %>%
  mutate(symbol = case_when(
    !is.na(df$symbol)~df$symbol,
    grepl(other,pattern='gene=')~df$other,
    grepl(other,pattern='product=')~df$other,
    T~NA))
df$symbol<-str_replace(df$symbol,pattern="product=",'')
df$other<-NULL

df$genetype<-gff$V3

write.csv(df,'genemap.csv')

genemap <- read.csv("./genemap.csv", row.names=2)
ID2gene<-tapply(genemap$symbol,rownames(genemap),paste)
gene2ID<-tapply(rownames(genemap),genemap$symbol,paste,collapse="; ")

########### update genemap with Refseq and GO term ###########
GO_list <- read.delim("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/LiLiang_neisseria/GO_list.txt")

GO_list<-GO_list[!duplicated(GO_list$ID),]
write.csv(GO_list,'GO_list.dedup.csv')


gff.sorted <- read.delim("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/LiLiang_neisseria/gff.sorted.txt")
gff.sorted$ID<-str_replace(gff.sorted$ID,'_','-')
gff.sorted<-arrange(gff.sorted,ID)

all(gff.sorted$ID==rownames(genemap))
all(gff.sorted$symbols==genemap$symbol)


########### check previous gene ###########
feature<-c('QFG95-10175','QFG95-10180', #comP
           'QFG95-01345', #nodQ
           'QFG95-03860',#cysP
           'QFG95-03270',
           'grpE',
           'QFG95-09300', #CT0075
           'QFG95-00250'
           )

comp1<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-10175')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('comP.QFG95-10175.pdf')
comp1
dev.off()



comp2<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-10180')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('comP.QFG95-10180.pdf')
comp2
dev.off()


nodQ<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-01345')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('nodQ.QFG95-01345.pdf')
nodQ
dev.off()


cysP<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-03860')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('cysP.QFG95-03860.pdf')
cysP
dev.off()


calcium.binding.protein<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-03270')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('calcium.pro.QFG95-03270.pdf')
calcium.binding.protein
dev.off()


grpE<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='grpE')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('grpE.pdf')
grpE
dev.off()



CT0075<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-09300')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('CT0075.QFG95-09300.pdf')
CT0075
dev.off()


QFG95.00250<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-00250')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('hypothetical.QFG95-00250.pdf')
QFG95.00250
dev.off()

########### 08 Nov 2024 check QFG95-10090 ###########

Neiss<-readRDS("Neiss.RDS")

QFG95.10090<-VlnPlot(Neiss,split.plot = T,split.by = 'orig.ident',feature='QFG95-10090')+
  scale_fill_manual(values = c('black','red'))

dev.new()
pdf('QFG95.10090.vln.pdf')
QFG95.10090
dev.off()


#heatmap



QFG95.10090.feature<-FeaturePlot(Neiss,features = 'QFG95-10090',pt.size = 1.5,split.by = 'orig.ident',slot='scale.data',
                                max.cutoff = '0.99',min.cutoff = '0.5',order=T,alpha = 0.9,
                                cols = c("grey80", "red"))

QFG95.10090.feature<-QFG95.10090.feature+
  labs(y = 'QFG95.10090')

dev.new()
pdf('QFG95.10090.feature.pdf',width=8,height = 5)
QFG95.10090.feature
dev.off()


Neiss.marker<- FindAllMarkers(Neiss, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
Neiss.marker$symbols<-ID2gene[Neiss.marker$gene]
Neiss.marker <- Neiss.marker %>%
  mutate(symbols=case_when(!is.na(symbols)~symbols,
                           T~gene))
Neiss.marker$symbols<-str_split_i(Neiss.marker$symbols,'\\.',1)
Neiss.marker$ID<-gene2ID[Neiss.marker$symbols]


Neiss.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5

hm.marker<-  DoHeatmap(subset(Neiss,downsample=500), features = top5$gene)

dev.new()
pdf("hm.marker2.pdf")
hm.marker
dev.off()

  
