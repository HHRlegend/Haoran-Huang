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

sce.all=run_seurat(sce.all.filt)
set.seed(4)
sce.all <- RunHarmony(sce.all,c( "orig.ident" ))
harmony_embeddings <- Embeddings(sce.all, 'harmony')
sce.all=RunUMAP(sce.all,reduction = "harmony",dims = 1:20)

mycolor1 <- getPalette(15)
pcompareumap=plot_grid(ncol = 1,
                     DimPlot(sce.all.filt, reduction = "umap", group.by = "orig.ident",cols=mycolor1)+NoAxes()+ggtitle("UMAP raw_data"),
                     DimPlot(sce.all, reduction = "umap", group.by = "orig.ident",cols=mycolor1)+NoAxes()+ggtitle("UMAP integrated")
)
pcompareumap
ggsave(filename = "pcompareumap.pdf", width = 6, height = 11)

sce.all=FindNeighbors(sce.all, reduction = "harmony",  dims = 1:20 )
for (res in c(0.1, 0.3,0.5,0.6,0.7,0.8,0.9,1,1.2,1.5)) {
  sce.all=FindClusters(sce.all,  
                       resolution = res, 
                       algorithm = 1)}
apply(sce.all@meta.data[,grep("RNA_snn_res",
                              colnames(sce.all@meta.data))],2,table)
ptree=clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
ptree
ggsave(plot=ptree, filename="Tree_diff_resolution.pdf",width = 10, height = 10)

mycolor1 <- getPalette(38)
DimPlot(sce.all, reduction = "umap",  group.by = "RNA_snn_res.1",cols=mycolor1,pt.size = 0.5,label = TRUE)+NoLegend()
ggsave(filename = "umap_res.pdf", width =5, height = 5)##############

pdf('orig.ident-vs-RNA_snn_res.pdf', width =10, height = 10 )
gplots::balloonplot(table(sce.all$RNA_snn_res.1,sce.all$orig.ident))
dev.off()

table(sce.all$RNA_snn_res.1)
My_levels<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26",'27','28','29','30','31')
sce.all$RNA_snn_res.1<-factor(sce.all$RNA_snn_res.1,levels= My_levels)
table(sce.all$RNA_snn_res.1)
sce.all@meta.data$seurat_clusters=sce.all@meta.data$RNA_snn_res.1

Epi =c('KRT18','PROX1','ALDH1A1')
Mac =c('CD163','CD68','VSIG4')  
Mono=c('S100A8', 'FCN1')  
DC=c("CLEC9A",'IDO1') 
B =c('MS4A1','CD79A')
Pla=c('MZB1','TNFRSF17')
T_NK  =c('CD2', 'CD3D', 'CD7')
Fibro=c('COL3A1','COL1A2','DCN')
Endo =c('VWF','CD34','CDH5')
Cyc=c( "CDK1",'MKI67','TOP2A')
genes_to_checkall =list(
Mac =Mac ,
Mono=Mono,
DC=DC,
B =B ,
Pla=Pla,
T_NK=T_NK,
Epi=Epi,
Fibro=Fibro,
Endo =Endo
)
genes_to_checkall  = lapply(genes_to_checkall , str_to_upper)

p_all_markers=DotPlot(sce.all , 
                      features = genes_to_checkall ,
                      assay='RNA',group.by = "RNA_snn_res.1" )+
  theme_bw()+
  scale_color_continuous(low="#A6CEE3FF",high =  "#FD6467")+
   theme(legend.position = "right",legend.box = "vertical",
   legend.margin=margin(t= 0, unit='cm'),
  axis.text.x  = element_text(color="black",size=10,angle = 45,vjust = 0.5, hjust=0.5),
    axis.text.y  = element_text(color="black",size=10),
   legend.text = element_text(size =10,color="black"),
    legend.title = element_text(size =10,color="black")
  ) ;p_all_markers
ggsave('allmarkers.pdf',height = 9,width = 9)

sce.all= sce.all[,sce.all@meta.data[["seurat_clusters"]]%in% c(0:20,22:24,26:27,29,31)]

celltype=data.frame(ClusterID=0:31,celltype= 0:31)   
celltype[celltype$ClusterID %in% c(1,14),2] <- 'T/NK cells'
celltype[celltype$ClusterID %in% c(12),2] <- 'Monocytes'
celltype[celltype$ClusterID %in% c(0,5,10,16,22),2] <- 'Macrophages'
celltype[celltype$ClusterID %in% c(18),2] <- 'Dendritic cells'
celltype[celltype$ClusterID %in% c(19),2] <- 'B cells'
celltype[celltype$ClusterID %in% c(13,26,31),2] <- 'Plasma cells'
celltype[celltype$ClusterID%in% c(9,23),2] <- 'Fibroblasts'
celltype[celltype$ClusterID %in% c(8,15,20,29),2] <- 'Endothelial cells'
celltype[celltype$ClusterID %in% c(2,3,4,6,7,11,17,24,27),2] <- 'Epithelial cells'

table(celltype$celltype)
sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.1== celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)

mycolor3 <- c('#FB9A99FF','#46ACC8','#FD6467','#FDBF6FFF','#CAB2D6FF','#B2DF8AFF','#7294D4','#D870A9B2','DEEPSKYBLUE2','#A73030B2','#725663B2','#8DBC80FF')###
dfall=sce.all@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=sce.all@meta.data$celltype)

celltype_positionall <- dfall%>%
  group_by(celltype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
   umap_2 = median(umap_2))

p_UMAPall=ggplot(dfall, aes(umap_1,umap_2, color=celltype))+
  geom_point(size = 0.01) + 
  theme(panel.border = element_blank(), 
   axis.title = element_blank(),  
   axis.text = element_blank(), 
   axis.ticks = element_blank(),
   panel.background = element_rect(fill = 'white'), 
   plot.background=element_rect(fill="white"),
   legend.title = element_blank(),
   legend.key=element_rect(fill='white'), 
   legend.text = element_text(size=25), 
   legend.key.size=unit(10,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 10)))+ 
  scale_color_manual(values=mycolor3)+theme_dr(xlength = 0.22, ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())+ 
theme(legend.position = "bottom")+
guides(color = guide_legend(override.aes = list(size = 5))) 
p_UMAPall
ggsave(filename = "UMAP1.pdf", width =7, height =7.5)
save(sce.all,file = 'sce.all1.Rdata') 

sceepi = sce.all[,sce.all@meta.data$celltype%in% 'Epithelial cells']###########
sceepi<- NormalizeData(sceepi, normalization.method = "LogNormalize", scale.factor = 1e4) 
sceepi<- FindVariableFeatures(sceepi, selection.method = 'vst', nfeatures = 2000)
sceepi <- ScaleData(sceepi)
sceepi<- RunPCA(sceepi, features = VariableFeatures(object = sceepi)) 
ElbowPlot(sceepi) 
sceepi  <- FindNeighbors(sceepi , dims = 1:20)
mycolor5 <- getPalette(9)
  
sceepi<- FindClusters(sceepi, resolution = 0.1)
table(sceepi@meta.data$seurat_clusters)

table(sce.all$celltype)
sce1=sce.all[,sce.all@meta.data$celltype=='Macrophages']
  tam.cells  <-  row.names(sce1@meta.data) 
 tam.cells=sample(tam.cells,500)
  tamMat=as.data.frame(GetAssayData(subset(sce1, cells=tam.cells),
                                     slot='counts',assay='RNA'))

sce2=sce.all[,sce.all@meta.data$celltype=='T/NK cells']
  ecs.cells  <-  row.names(sce2@meta.data) 
  ecs.cells=sample(ecs.cells,500)
  ecsMat=as.data.frame(GetAssayData(subset(sce2, cells=ecs.cells),
                                     slot='counts',assay='RNA'))

 Macrophages = tamMat
  T_NK_cells= ecsMat

epiMat=as.data.frame(GetAssayData(sceepi,
                                  slot='counts',assay='RNA'))
ids = intersect(rownames(epiMat),rownames(  Macrophages))
this_dat=cbind(epiMat[ids,],  Macrophages[ids,],  T_NK_cells[ids,])
groupinfo=data.frame(v1=colnames(this_dat),
                     v2=c( sceepi@meta.data$seurat_clusters,
                           rep('Macrophages',500),
                           rep('T_NK_cells',500)))
head(groupinfo) 
groupFiles='groupFiles.txt'
write.table(groupinfo,file = groupFiles,
            sep = '\t',quote = F,col.names = F,row.names = F)
print(dim(this_dat))

  geneInfor=annoGene(rownames(this_dat),"SYMBOL",'human')
  colnames(geneInfor)
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  length(unique(geneInfor[,1]))
  head(geneInfor)

  dat=this_dat[rownames(this_dat) %in% geneInfor[,1],]
  dat=dat[match( geneInfor[,1], rownames(dat) ),]
  dim(dat)
  expFile='expFile.txt'
  library(data.table)
  fwrite(dat, file = expFile,  row.names = TRUE,sep = "\t", quote = FALSE)
  
  head(geneInfor)
  geneFile='geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = FALSE,col.names =FALSE,row.names = FALSE)
  
  expFile='expFile.txt'
  groupFiles='groupFiles.txt'
  geneFile='geneFile.txt'
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c('Macrophages',
                                                        'T_NK_cells'))  
  
  infercnv_obj2 = infercnv::run(infercnv_obj,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                                cluster_by_groups=TRUE ,   # cluster
                                denoise=TRUE)

infer_CNV_obj<-readRDS('./infercnv_output/run.final.infercnv_obj')
expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)

if(TRUE){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$'Macrophages']
  tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$'T_NK_cells']
  tmp= cbind(tmp1,tmp2)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
   
  cnv_score_table_pts[1:4,1:4]
  str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 
}
head(cell_scores_CNV) 
score=cell_scores_CNV
head(score)

My_level<-c("0","1","2","3","4","5","6","7","8","9","10","11")
sceepi$seurat_clusters<-factor(sceepi$seurat_clusters,levels= My_level)
sceepi@meta.data$Cells<- sceepi@meta.data$seurat_clusters
sce1@meta.data$Cells<- sce1@meta.data$celltype
sce2@meta.data$Cells<- sce2@meta.data$celltype
sce3=merge(sceepi,y=c(sce2,sce1))
meta3 <- sce3@meta.data
meta3$TotalCNV = score[match(colnames(sce3),
                            rownames(score)),1] 
meta3$Cells=factor(meta3$Cells,levels=c('Macrophages',My_level,'T/NK cells'))
mycolor10 <- getPalette(25)
p1=ggplot(meta3, aes(x=Cells , y=TotalCNV, fill=Cells )) +  geom_violin(color="white") +  geom_boxplot(width=0.1,color='white',outlier.shape = NA)  +
  scale_fill_manual(values = mycolor10)+labs(x=NULL,title=NULL)+theme_classic()
p1
ggsave('Infercnv_score.pdf',width = 10,height =4)

 celltypeepi1=data.frame(ClusterID=0:11,celltypeepi1=0:11)    ####
  celltypeepi1[celltypeepi1$ClusterID %in% c(0:11),2]='Malignant cells'  #####

  table(celltypeepi1$celltypeepi1)
  sceepi@meta.data$celltypeepi1= "NA"  
  for(i in 1:nrow(celltypeepi1)){
    sceepi@meta.data[which(sceepi@meta.data$RNA_snn_res.0.1== celltypeepi1$ClusterID[i]),'celltypeepi1'] <- celltypeepi1$celltypeepi1[i]}
  table(sceepi@meta.data$celltypeepi1)
  
metaepi1<- sce.all@meta.data
metaepi1$cell_typeepi1 <- NA
metaepi1$cell_typeepi1 <- sceepi@meta.data[match(rownames(sce.all@meta.data), rownames(sceepi@meta.data)), 'celltypeepi1']
print(table(metaepi1$cell_typeepi1))
DT::datatable(metaepi1)
celltype_namesepi1<- NULL
for(i in 1:dim(metaepi1)[1]){
  sub_dataepi1<- metaepi1[i,]
  if(is.na(sub_dataepi1$cell_typeepi1)){
    sub_dataepi1$cell_typeepi1<- sub_dataepi1$celltype
    celltype_namesepi1<- c(celltype_namesepi1, sub_dataepi1$celltype)
  }else{
    celltype_namesepi1<- c(celltype_namesepi1, sub_dataepi1$cell_typeepi1)
  }
}
print(table(celltype_namesepi1))
metaepi1$celltype <- celltype_namesepi1
sce.all@meta.data <- metaepi1

mycolor3 <- c('#FB9A99FF','#46ACC8','#FD6467','#FDBF6FFF','#CAB2D6FF','#B2DF8AFF','#7294D4','#D870A9B2','DEEPSKYBLUE2','#A73030B2','#725663B2','#8DBC80FF')
dfall=sce.all@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=sce.all@meta.data$celltype)

celltype_positionall <- dfall%>%
  group_by(celltype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
   umap_2 = median(umap_2))

p_UMAPall1=ggplot(dfall, aes(umap_1,umap_2, color=celltype))+
  geom_point(size = 0.01) +
  theme(panel.border = element_blank(), 
   axis.title = element_blank(),  
   axis.text = element_blank(), 
   axis.ticks = element_blank(),
   panel.background = element_rect(fill = 'white'), 
   plot.background=element_rect(fill="white"),
   legend.title = element_blank(), 
   legend.key=element_rect(fill='white'), 
   legend.text = element_text(size=25), 
   legend.key.size=unit(10,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 10)))+
  scale_color_manual(values=mycolor3)+theme_dr(xlength = 0.22, ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())+ 
theme(legend.position = "bottom")+
guides(color = guide_legend(override.aes = list(size = 5))) 
p_UMAPall1
ggsave(filename = "UMAP2.pdf", width =8, height =8.5)

rt=read_table("allexp.txt")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
bulk_dataset=avereps(data)

phenotype=read.table("phenotype_sur.txt",row.names=1,header=T)
colnames(bulk_dataset) == row.names(phenotype)
phenotype=phenotype[,1:2]

sc_dataset=sce.all
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.02, 
       family = "cox", Save_file ='Scissor_survival.RData')
Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "Scissor")

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "222"
Scissor_select[infos1$Scissor_neg] <- "111"
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "Scissor1")
UMAP_scissor <- DimPlot(sc_dataset, reduction = 'umap', 
        group.by = 'Scissor',
        cols = c('#88c4e8','#93cc82','#ea9c9d'), 
        pt.size = 0.001, order = c("Scissor+","Scissor-")
)+theme(legend.position=c(.7,.2))
UMAP_scissor
ggsave(filename = "UMAP_scissor.pdf", width =5, height =5)
sce.all=sc_dataset
save(sce.all,file = 'sce.all2.Rdata') 

pdf('celltypeall-vs-scissor.pdf', width =8, height =8 )
gplots::balloonplot(table(sce.all$Scissor,sce.all$celltype))
dev.off()

load('Scissor_survival.RData')
numbers <- length(infos1$Scissor_pos) + length(infos1$Scissor_neg)
result1 <- reliability.test(X, Y, network, alpha = 0.02, family = "cox", cell_num = numbers, n = 10, nfold = 10)

source("./custom_seurat_functions_scissor.R")
per1 = plot.clusters.group(data =sc_dataset,clusters = "Scissor", xlab = "Scissor",log =F, group = "celltype",legend.title = "Celltype",widths = c(3,1),color = 1)
per1
ggsave(filename = "percent_scissor.pdf", width =10, height =5)

sceepi = sce.all[,sce.all@meta.data$celltype%in% 'Malignant cells']###########
sceepi<- NormalizeData(sceepi, normalization.method = "LogNormalize", scale.factor = 1e4) 
sceepi<- FindVariableFeatures(sceepi, selection.method = 'vst', nfeatures = 2000)
sceepi <- ScaleData(sceepi)
sceepi<- RunPCA(sceepi, features = VariableFeatures(object = sceepi)) 
sceepi  <- FindNeighbors(sceepi , dims = 1:20)
sceepi=RunUMAP(sceepi,reduction = "harmony",dims = 1:20)

sc_dataset1=sceepi
rt=read_table("TCGAMVIexp.txt")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
bulk_dataset1=avereps(data)

phenotype1=read.table("phenotype_mvi.txt",row.names=1,header=T)
colnames(bulk_dataset1) == row.names(phenotype1)
phenotype1=phenotype1[,1]
tag=c('MVI_negative','MVI_positive')

infos2 <- Scissor(bulk_dataset1, sc_dataset1, phenotype1, alpha = 0.05, tag=tag, family = "binomial", Save_file ='Scissor_MVI.RData')
Scissor_select1 <- rep(0, ncol(sc_dataset1))
names(Scissor_select1) <- colnames(sc_dataset1)
Scissor_select1[infos2$Scissor_pos] <- "MVI_Scissor+"
Scissor_select1[infos2$Scissor_neg] <- "MVI_Scissor-"
sc_dataset1 <- AddMetaData(sc_dataset1, metadata = Scissor_select1, col.name = "MVI_Scissor")
Scissor_select1 <- rep(0, ncol(sc_dataset1))
names(Scissor_select1) <- colnames(sc_dataset1)
Scissor_select1[infos2$Scissor_pos] <- "22"
Scissor_select1[infos2$Scissor_neg] <- "11"
sc_dataset1 <- AddMetaData(sc_dataset1, metadata = Scissor_select1, col.name = "MVI_Scissor1")
sceepi=sc_dataset1

UMAP_scissor1 <- DimPlot(sc_dataset1, reduction = 'umap', 
        group.by = 'MVI_Scissor',
        cols = c('#ffff99ff','#7294d4','#db6968'), 
        pt.size = 0.001, order = c("MVI_Scissor+","MVI_Scissor-"))+theme(legend.position=c(.6,.8))
UMAP_scissor1
ggsave(filename = "UMAP_MVI+scissor.pdf", width =5, height =5)

load('Scissor_MVI.RData')
numbers1 <- length(infos2$Scissor_pos) + length(infos2$Scissor_neg)
result2 <- reliability.test(X, Y, network, alpha = 0.05, family = "binomial", cell_num = numbers1, n = 10, nfold = 10)
save(sceepi,file = 'sceepi1.Rdata') 

markers_scissormvi_HCC <- FindMarkers(sceepi, ident.1 = "MVI_Scissor+", 
                       group.by = 'MVI_Scissor', 
                       logfc.threshold = 0.585)
markers_scissormvi_HCC1 <- markers_scissormvi_HCC[which(markers_scissormvi_HCC$p_val_adj<0.05),]
write.table(markers_scissormvi_HCC1,file="markers_scissormvi_HCC.txt",sep="\t",row.names=T,quote=F)

markers_scissor_HCC <- FindMarkers(sceepi, ident.1 = "Scissor+", 
                       group.by = 'Scissor', 
                       logfc.threshold = 0.585)
markers_scissor_HCC1 <- markers_scissor_HCC[which(markers_scissor_HCC$p_val_adj<0.05),]
write.table(markers_scissor_HCC1,file="markers_scissor_HCC.txt",sep="\t",row.names=T,quote=F)

markers_SCISSOR_HCC <- FindMarkers(sceepi, ident.1 = "Scissor+", ident.2= "Scissor-", 
                       group.by = 'Scissor', 
                       logfc.threshold = 0)
markers_SCISSOR_HCC <- markers_SCISSOR_HCC [which(markers_SCISSOR_HCC $p_val_adj<0.05),] 
write.table(markers_SCISSOR_HCC ,file="markers_SCISSOR_HCC.txt",sep="\t",row.names=T,quote=F)

sparse_data <- as(as.matrix(sceepi@assays$RNA@counts),'sparseMatrix')
mdata <- new('AnnotatedDataFrame',data=sceepi@meta.data)
fData <- data.frame(gene_short_name=row.names(sparse_data),row.names = row.names(sparse_data))
fd <- new('AnnotatedDataFrame',data=fData)
cds <- newCellDataSet(cellData = sparse_data,
                              phenoData = mdata,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds,min_expr = 0.1)

express_genes <- VariableFeatures(sceepi)
 cds1<- setOrderingFilter(cds, express_genes)
 plot_ordering_genes(cds1)

cds <- reduceDimension(cds1, #residualModelFormulaStr = "~study",
reduction_method = 'DDRTree')

source('order_cells.R')
cds <- orderCells(cds,root_state=3)

colour=c('#88c4e8','#93cc82','#ea9c9d')
colour1=c("#9370DB","#FDBF6FFF","#8DBC80FF","#7CFC00","#FFFF00",  
         "#808000","#FA8072","#7B68EE","#9400D3","#800080")

data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% 
  select_(Component_1 = 1, Component_2 = 2) %>% 
  rownames_to_column("cells") %>% 
  mutate(pData(cds)$State) %>% 
  mutate(pData(cds)$Pseudotime, 
         pData(cds)$group, pData(cds)$Scissor,pData(cds)$MVI_Scissor,
         pData(cds)$celltype)
colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","Group",'Scissor','MVI_Scissor',"celltype")

dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))

edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% select_(source = "sample_name", 
                                     source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
  left_join(ica_space_df %>% select_(target = "sample_name", 
                                     target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

 Cellratio1 <- prop.table(table(data_df$Group,data_df$State), margin = 2)
 Cellratio1 <- as.data.frame(Cellratio1)
 colnames(Cellratio1) <- c("Group",'State',"Freq1")
Cellratio2 <- prop.table(table(data_df$Scissor,data_df$State ), margin = 2)
 Cellratio2 <- as.data.frame(Cellratio2)
 colnames(Cellratio2) <- c('Group','State',"Freq2")
Cellratio3 <- prop.table(table(data_df$State,data_df$MVI_Scissor), margin = 2)
 Cellratio3 <- as.data.frame(Cellratio3)
 colnames(Cellratio3) <- c('State',"Group","Freq3")

g1 <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                 y = Component_2,
                                 color =Pseudotime)) + 
  scale_color_viridis()+
 geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
 theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g1
ggsave("Pseudotime1.pdf",width =6,height =4 )

mycol <- c4a('pastel',7)
mycol2 <- c4a('superfishel_stone',7)
g2 <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                 y = Component_2,
                                 color =State)) +
scale_colour_manual(values = c('#88c4e8','#a1d5b9','#ea9c9d',"#9370DB","#FDBF6FFF","#8DBC80FF","#7CFC00"))+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g2
ggsave("Pseudotime2.pdf",width =5.5,height =4 )

pp1=ggplot(Cellratio1) +
  geom_col(aes(x = 3,
               y = Freq1,
               fill = Group),
           width = 1.5,
           color = 'white') +
facet_grid(.~State)+
 coord_polar(theta = "y") +
  xlim(c(0.2, 3.8))+
scale_fill_manual(values = mycol) +
  theme_void()
pp1
ggsave("Pseudotime2131.pdf",width =3,height =1)

pp2=ggplot(Cellratio2) +
  geom_col(aes(x = 3,
               y = Freq2,
               fill = Group),
           width = 1.5,
           color = 'white') +
facet_grid(.~State)+
 coord_polar(theta = "y") +
  xlim(c(0.2, 3.8))+
scale_fill_manual(values = mycol2) +
  theme_void()
pp2
ggsave("Pseudotime213.pdf",width =3,height =1)

g3<- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                 y = Component_2,
                                 color =MVI_Scissor)) + 
scale_colour_manual(values = c('#a3d393','#f8984e','#db6968',"#9370DB","#FDBF6FFF","#8DBC80FF","#7CFC00"))+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g3
ggsave("Pseudotime3.pdf",width =6,height =4 )

pp3=ggplot(Cellratio3) +
  geom_col(aes(x = 3,
               y = Freq3,
               fill =State ),
           width = 1.5,
           color = 'white') +
facet_grid(.~Group)+
 coord_polar(theta = "y") +
  xlim(c(0.2, 3.8))+
scale_fill_manual(values = c('#88c4e8','#93cc82','#ea9c9d',"#9370DB","#FDBF6FFF","#8DBC80FF","#7CFC00")) +
  theme_void()
pp3
ggsave("Pseudotime21.pdf",width =3,height =1)

#CytoTRACE
monocle_meta <- data.frame(t(cds @reducedDimS), 
                         cds $Pseudotime, 
                           cds $State, 
                           cds $MVI_Scissor)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "MVI_Scissor")
phenot1 <- monocle_meta$MVI_Scissor
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle <- monocle_meta[,1:2]
mat<-as.matrix(sceepi@assays$RNA@counts)
results <- CytoTRACE(mat = mat)
getwd()
plotCytoTRACE(results, phenotype = phenot1, emb = emb_monocle)

BEAM_res <- BEAM(cds, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
df <-plot_pseudotime_heatmap2(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                             num_clusters = 3,
                             cores = 8,
                             use_gene_short_name = T,
                             show_rownames = T)
enrich <- enrichCluster(object = df,
                        OrgDb = org.Hs.eg.db,
                        type = "CC",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 10,
                        seed = 5201314)
markGenes = c('HLA-DRA1','HLA-DPA1','HLA-DQA1','HLA-DRA','HLA-DPB1','HLA-DQB1',
'CD2BP2','SF3B1','SRSF2','PCBP1','PCBP2','SFPQ',
'RPL12','RPL3','RPL36','RPL4','RPL31','RPL21','RPL5','RPL11','RPS5','RPS27A','RPS3A','RPS8','RPS15',
'CD63','CD9','GAPDH','A2M','TUBB','TUBA1B','HSPA8','RAC1','ACTB','ANXA2','SDCBP','TSG101','LDHA'
,'MT-CO1','MT-ATP6','VEGFA','ITGB1','TMEM45A','LGALS3','SERPINA12','NDUFA6','S100A10','STMN1','RPL8','BHMT','CYP3A5','CFHR3'
)
pdf('branch-enrich.pdf',height = 9,width = 16,onefile = F)
visCluster(object = df,
           plot.type = "both",
           cluster.order = c(1,3,2),
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           go.col = rep(jjAnno::useMyCol("calm",n =3),each =10),
           add.bar = T,
           line.side = "left")
dev.off()
saveRDS(cds,file = 'monocle.rds')

scemac = sce.all[,sce.all@meta.data$celltype%in% 'Macrophages']###########
scemac <- NormalizeData(scemac , normalization.method = "LogNormalize", scale.factor = 1e4) 
scemac <- FindVariableFeatures(scemac , selection.method = 'vst', nfeatures = 2000)
scemac <- ScaleData(scemac )
scemac <- RunPCA(scemac , features = VariableFeatures(object = scemac )) 
ElbowPlot(scemac ) 
scemac   <- FindNeighbors(scemac , dims = 1:20)

markers_SCISSOR_MAC <- FindMarkers(scemac, ident.1 = "Scissor+", ident.2= "Scissor-", 
                       group.by = 'Scissor', 
                       logfc.threshold = 0)
markers_SCISSOR_MAC <- markers_SCISSOR_MAC[which(markers_SCISSOR_MAC $p_val_adj<0.05),] 
write.table(markers_SCISSOR_MAC,file="markers_SCISSOR_MAC.txt",sep="\t",row.names=T,quote=F)

sceTNK = sce.all[,sce.all@meta.data$celltype%in% 'T/NK cells']###########
sceTNK <- NormalizeData(sceTNK  , normalization.method = "LogNormalize", scale.factor = 1e4) 
sceTNK<- FindVariableFeatures(sceTNK  , selection.method = 'vst', nfeatures = 2000)
sceTNK  <- ScaleData(sceTNK)
sceTNK <- RunPCA(sceTNK  , features = VariableFeatures(object = sceTNK )) 
ElbowPlot(sceTNK) 
sceTNK<- FindNeighbors(sceTNK  , dims = 1:20)

markers_SCISSOR_TNK<- FindMarkers(sceTNK, ident.1 = "Scissor+", ident.2= "Scissor-", 
                       group.by = 'Scissor', 
                       logfc.threshold =0)
markers_SCISSOR_TNK<- markers_SCISSOR_TNK[which(markers_SCISSOR_TNK$p_val_adj<0.05),] 
write.table(markers_SCISSOR_TNK,file="markers_SCISSOR_TNK.txt",sep="\t",row.names=T,quote=F)

sceFIB = sce.all[,sce.all@meta.data$celltype%in% 'Fibroblasts']###########
sceFIB<- NormalizeData(sceFIB  , normalization.method = "LogNormalize", scale.factor = 1e4) 
sceFIB<- FindVariableFeatures(sceFIB , selection.method = 'vst', nfeatures = 2000)
sceFIB  <- ScaleData(sceFIB)
sceFIB <- RunPCA(sceFIB , features = VariableFeatures(object = sceFIB)) 
ElbowPlot(sceFIB) 
sceFIB<- FindNeighbors(sceFIB  , dims = 1:20)

markers_SCISSOR_FIB<- FindMarkers(sceFIB, ident.1 = "Scissor+", ident.2= "Scissor-", 
                       group.by = 'Scissor', 
                       logfc.threshold = 0)
markers_SCISSOR_FIB<- markers_SCISSOR_FIB[which(markers_SCISSOR_FIB$p_val_adj<0.05),] 
write.table(markers_SCISSOR_FIB,file="markers_SCISSOR_FIB.txt",sep="\t",row.names=T,quote=F)

#ENDO
sceENDO = sce.all[,sce.all@meta.data$celltype%in% 'Endothelial cells']###########
sceENDO<- NormalizeData(sceENDO  , normalization.method = "LogNormalize", scale.factor = 1e4) 
sceENDO<- FindVariableFeatures(sceENDO , selection.method = 'vst', nfeatures = 2000)
sceENDO  <- ScaleData(sceENDO)
sceENDO <- RunPCA(sceENDO, features = VariableFeatures(object = sceENDO)) 
ElbowPlot(sceENDO) 
sceENDO<- FindNeighbors(sceENDO  , dims = 1:20)

markers_SCISSOR_ENDO<- FindMarkers(sceENDO, ident.1 = "Scissor+", ident.2= "Scissor-", 
                       group.by = 'Scissor', 
                       logfc.threshold = 0)
markers_SCISSOR_ENDO<- markers_SCISSOR_ENDO[which(markers_SCISSOR_ENDO$p_val_adj<0.05),] 
write.table(markers_SCISSOR_ENDO,file="markers_SCISSOR_ENDO.txt",sep="\t",row.names=T,quote=F)

VOL=read.table('markers_SCISSOR_ALL.txt',sep='\t',header=T)
mycolor <- getPalette(5)
jjVolcano(diffData =VOL,log2FC.cutoff=1,
          tile.col =mycolor  ,aesCol = c('#7294D4','#FD6467'))
ggsave(filename = "VOL-SCISSOR.pdf", width =8, height =8)

scemvi= sceepi[,sceepi@meta.data$MVI_Scissor %in% 'MVI_Scissor+']###########
scemvi<- NormalizeData(scemvi, normalization.method = "LogNormalize", scale.factor = 1e4) 
scemvi<- FindVariableFeatures(scemvi, selection.method = 'vst', nfeatures = 2000)
scemvi <- ScaleData(scemvi)
scemvi<- RunPCA(scemvi, features = VariableFeatures(object = scemvi)) 
scemvi <- FindNeighbors(scemvi, dims = 1:20)
 scemvi@meta.data$celltypeepi= "NA"  
 scemvi@meta.data[which( scemvi@meta.data$MVI_Scissor== 'MVI_Scissor+'),'celltypeepi'] <-  'MVI_Scissor+ malignant cells'
 table( scemvi@meta.data$celltypeepi)

sceother= sce.all[,sce.all@meta.data$celltype %in% c('Macrophages','B cells','Plasma cells','Dendritic cells','T/NK cells','Endothelial cells','Fibroblasts','Monocytes')]###########
sceother <- NormalizeData(sceother , normalization.method = "LogNormalize", scale.factor = 1e4) 
sceother<- FindVariableFeatures(sceother, selection.method = 'vst', nfeatures = 2000)
sceother <- ScaleData(sceother )
sceother <- RunPCA(sceother, features = VariableFeatures(object = sceother )) 
ElbowPlot(sceother) 
sceother <- FindNeighbors(sceother , dims = 1:20)

scecellchat = merge(scemvi,sceother)#####
scecellchat<- NormalizeData(scecellchat, normalization.method = "LogNormalize", scale.factor = 1e4) 
scecellchat <- FindVariableFeatures(scecellchat, selection.method = 'vst', nfeatures = 2000)
scecellchat <- ScaleData(scecellchat )
scecellchat <- RunPCA(scecellchat, features = VariableFeatures(object = scecellchat )) 
ElbowPlot(scecellchat ) 
scecellchat   <- FindNeighbors(scecellchat , dims = 1:20)

metacellchat  <- scecellchat @meta.data
metacellchat  $celltypecellchat   <- NA
metacellchat  $celltypecellchat <- scecellchat@meta.data[match(rownames(scecellchat@meta.data), rownames(scemvi@meta.data)), 'celltypeepi']
print(table(metacellchat$celltypecellchat))
#DT::datatable(metacellchat )
celltype_namescellchat  <- NULL
for(i in 1:dim(metacellchat )[1]){
  sub_datacellchat  <- metacellchat [i,]
  if(is.na(sub_datacellchat $celltypecellchat )){
    sub_datacellchat $celltypecellchat <- sub_datacellchat $celltypecellchat 
    celltype_namescellchat  <- c(celltype_namescellchat , sub_datacellchat $celltypecellchat )
  }else{
    celltype_namescellchat <- c(celltype_namescellchat , sub_datacellchat $celltypecellchat )
  }
}
print(table(celltype_namescellchat ))
metacellchat $celltypecellchat  <- celltype_namescellchat 
scecellchat @meta.data <- metacellchat 

metaother <- scecellchat @meta.data
metaother   $celltypeother   <- NA
metaother   $celltypeother   <- sceother  @meta.data[match(rownames(scecellchat@meta.data), rownames(sceother  @meta.data)), 'celltype']
print(table(metaother  $celltypeother  ))
celltype_namesother  <- NULL
for(i in 1:dim(metaother  )[1]){
  sub_dataother  <- metaother  [i,]
  if(is.na(sub_dataother  $celltypeother  )){
    sub_dataother  $celltypecellchat<- sub_dataother  $celltypecellchat
    celltype_namesother  <- c(celltype_namesother  , sub_dataother  $celltypecellchat)
  }else{
    celltype_namesother  <- c(celltype_namesother  , sub_dataother  $celltype)
  }
}
print(table(celltype_namesother  ))
metaother  $celltypecellchat<- celltype_namesother  
scecellchat@meta.data <- metaother  

data.input <- scecellchat@assays$RNA@data
meta.data <- scecellchat@meta.data

cellchat <- createCellChat(object=data.input,
                           meta = meta.data,
                           group.by='celltypecellchat')
cellchat <- addMeta(cellchat,meta = meta.data)

cellchatDB <- CellChatDB.human
cellchat@DB <- cellchatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat,min.cells = 10)

df.net <- subsetCommunication(cellchat)
df.pathway <- subsetCommunication(cellchat,slot.name = 'netP')
write.csv(df.net, "cell-cell_communications.all.csv")
write.csv(df.pathway , "cell-cell_communicationspath.all.csv")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,
                 weight.scale = TRUE,label.edge = FALSE,
                 title.name = 'Number of Interaction')
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,
                 weight.scale = TRUE,label.edge = FALSE,
                 title.name = 'Interaction Weight')

mat <- cellchat@net$weight
par(mfrow = c(2,3),xpd=TRUE)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,
                   weight.scale = TRUE,edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}

cellchat@netP$pathways
pathway.show <- c('MIF')#####

par(mfrow=c(1,1),xpd=TRUE)
netVisual_aggregate(cellchat,signaling = pathway.show,
                    layout = 'circle')

par(mfrow=c(1,1))
netVisual_heatmap(cellchat,
                  signaling = pathway.show,
                  color.heatmap = c("white", "#b2182b"))

plotGeneExpression(cellchat, signaling = 'MIF', type = 'violin')

levels(cellchat@idents)
netVisual_bubble(cellchat,
                 sources.use = 7,
                 targets.use = 3, angle.x = 0,
                 signaling = c("MK","VEGF",'SPP1','CD99','FN1','GDF'), 
                 remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:9), signaling = c("MIF"), remove.isolate = FALSE)

sub.markers <- FindMarkers(sceepi,group.by = 'MVI_Scissor',
                           ident.1 = 'MVI_Scissor+',logfc.threshold = 0.01)
sub.markers.sig <- subset(sub.markers, p_val_adj<0.05 & avg_log2FC >0.5)
mydata <- data.frame(Gene=rownames(sub.markers),logFC=sub.markers$avg_log2FC)%>%
    arrange(desc(logFC))
genelist <- bitr(mydata$Gene,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")%>%
    inner_join(mydata,by= c('SYMBOL'='Gene'))%>%
   arrange(desc(logFC))
gsea_input <- genelist$logFC
names(gsea_input) <- genelist$ENTREZID
geneset <- read.gmt("h.all.v2023.2.Hs.entrez.gmt")  
gg <- GSEA(gsea_input, TERM2GENE=geneset,verbose=F,
            pvalueCutoff=0.05, pAdjustMethod = "BH")
sortgg<- gg[order(gg$NES, decreasing = T),]
sortgg<- sortgg[sortgg$p.adjust <0.05,]
write.csv(sortgg , "gsea.csv")
mycolor <- getPalette(9)
pdf('gsea.pdf',height = 4,width = 11)
dotplotGsea(data = gg,
                             topn=8,
                             order.by = "NES",
                              add.seg = T,
                              line.col =mycolor ,
                              line.type = 1
                  )
 dev.off()


