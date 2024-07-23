 setwd('~/')
library(AnnoProbe)
library(Cairo)
library(CellChat)
library(celldex)
library(circlize)
library(ClusterGVis)
library(clusterProfiler)
library(clustree)
library(cols4all)
library(ComplexHeatmap)
library(cowplot)
library(CytoTRACE)
library(data.table)
library(DESeq2) 
library(DoubletFinder)
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggsci)
library(ggVolcano)
library(gplots)
library(GseaVis)
library(GSVA)
library(harmony)
library(HGNChelper)
library(igraph)
library(infercnv)
library(limma)
library(monocle)
library(msigdbr)
library(NMF)
library(openxlsx)
library(org.Hs.eg.db)
library(paletteer)
library(patchwork)
library(pheatmap)
library(preprocessCore)
library(RColorBrewer)
library(ReactomeGSA)
library(readr)
library(reshape2)
library(rsvd)
library(scAB)
library(Scissor)
library(scRNAtoolVis)
library(Seurat)
library(singleseqgset)
library(stringr)
library(tidydr)
library(tidyverse)
library(VennDiagram)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

matrix_data <- read.table("GSE149614_HCC.scRNAseq.S71915.count.txt", sep="\t", header=T, row.names=1)
dim(matrix_data)

scRNA <- CreateSeuratObject(counts = matrix_data,project="os")
metadata<- read.table("HCC.metadata.txt", sep="\t", header=T)

T = metadata %>% filter(site == 'Tumor') %>% select (Cell) %>% unlist() %>% as.vector()
N = metadata %>% filter(site =='Normal') %>% select (Cell) %>% unlist() %>% as.vector()
P= metadata %>% filter(site =='PVTT') %>% select (Cell) %>% unlist() %>% as.vector()
sce.all <- scRNA[,c(T,N,P)]
dim(sce.all)

metadata1<- subset(metadata,metadata$site == 'Tumor')
metadata2<- subset(metadata,metadata$site == 'Normal')
metadata3<- subset(metadata,metadata$site == 'PVTT')
metadata=rbind(metadata1,metadata2,metadata3)
sce.all@meta.data$group=as.factor(metadata$site)
sce.all@meta.data$orig.ident=as.factor(metadata$sample)
########################################  1-QC   #################################################

dir.create("./1-QC")
setwd("./1-QC")

mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] 
mito_genes 
sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")

ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)

rownames(sce.all)[grep("^HB[^(p)]", rownames(sce.all))]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(p)]", col.name = "percent_hb")

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident",pt.size = 0,  features = feats, ncol = 2) + 
    NoLegend()
p1
ggsave(filename="Vlnplot1.pdf",plot=p1,width = 10,height = 5)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
p2 
ggsave(filename="Vlnplot2.pdf",plot=p2, width = 18,height = 5)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="Scatterplot.pdf",plot=p3,width = 6,height = 5)

selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 300 & nCount_RNA>1000)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]
sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all.filt) 
table(sce.all.filt@meta.data$orig.ident) 

selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 20)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 0.1)
length(selected_hb)
length(selected_ribo)
length(selected_mito)

sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)
table(sce.all.filt$orig.ident) 

feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
p1_filtered
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = 10,height = 5)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
    NoLegend()
p2_filtered
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered, width = 18,height = 5)

sce.all.filt <- NormalizeData(sce.all.filt,verbose=FALSE)
sce.all.filt <- FindVariableFeatures(sce.all.filt)
sce.all.filt <- ScaleData(sce.all.filt,verbose=FALSE)
sce.all.filt <- RunPCA(sce.all.filt,verbose=FALSE)
ElbowPlot(sce.all.filt , ndims=50)
sce.all.filt <- RunUMAP(sce.all.filt, dims = 1:30,verbose=FALSE)
sce.all.filt <- RunTSNE(sce.all.filt, dims = 1:30,verbose=FALSE)

sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
   theme_minimal()
ggsave(filename="cycle_details.pdf", width = 5,height = 4 )

dim.usage <- 30
sce.all.list <- SplitObject(sce.all.filt, split.by = "orig.ident")
sce.all.list
phe_lt <- lapply(names(sce.all.list),function(x){
sce.all.filt=sce.all.list[[x]]
sce.all.filt = NormalizeData(sce.all.filt)
sce.all.filt = FindVariableFeatures(sce.all.filt)
sce.all.filt = ScaleData(sce.all.filt, 
                           vars.to.regress = c("nFeature_RNA", "percent_mito"))
sce.all.filt = RunPCA(sce.all.filt)
sce.all.filt<- FindNeighbors(sce.all.filt, dims = 1:dim.usage) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:dim.usage) %>% 
  RunTSNE(dims = 1:dim.usage)
  sweep.res.list <- paramSweep_v3(sce.all.filt, PCs = 1:dim.usage, sct = FALSE) # 若使用SCT方法 标准化则'sct=T'
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]$pK))
  DoubletRate = ncol(sce.all.filt)*8*1e-6 
  nExp_poi <- round(DoubletRate*ncol(sce.all.filt))
  homotypic.prop <- modelHomotypic(sce.all.filt@meta.data$seurat_clusters) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

 sce.all.filt <- doubletFinder_v3(sce.all.filt, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
 
  DF.name = colnames(sce.all.filt@meta.data)[grepl("DF.classification", 
                                                   colnames(sce.all.filt@meta.data))]
  p5.dimplot=cowplot::plot_grid(ncol = 2, DimPlot(sce.all.filt, group.by = "orig.ident") + NoAxes(), 
                                DimPlot(sce.all.filt, group.by = DF.name) + NoAxes())
  p5.dimplot
  ggsave(filename=paste0("doublet_dimplot_",x,".pdf"),
         plot=p5.dimplot)
  
  p5.vlnplot=VlnPlot(sce.all.filt, features = "nFeature_RNA", 
                     group.by = DF.name, pt.size = 0.1)
  p5.vlnplot
  ggsave(paste0("doublet_vlnplot_",x,".pdf"),
         plot=p5.vlnplot)
  print(table(sce.all.filt@meta.data[, DF.name] ))
  phe=sce.all.filt@meta.data
  phe
})

kpCells=unlist(lapply(phe_lt, function(x){
  table(x[,ncol(x)])
  rownames(x[ x[,ncol(x)]=='Singlet', ])
}))
kp = colnames(sce.all.filt) %in% kpCells
table(kp)
sce.all.filt=sce.all.filt[,kp]

sce.all.filt <- NormalizeData(sce.all.filt,verbose=FALSE)
sce.all.filt <- FindVariableFeatures(sce.all.filt, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
sce.all.filt <- ScaleData(sce.all.filt,verbose=FALSE)
sce.all.filt <- RunPCA(sce.all.filt,verbose=FALSE)
ElbowPlot(sce.all.filt , ndims=50)
sce.all.filt <- RunUMAP(sce.all.filt, dims = 1:30,verbose=FALSE)
sce.all.filt <- RunTSNE(sce.all.filt, dims = 1:30,verbose=FALSE)

save(sce.all.filt,file = 'sce.all.filt2.Rdata') 


