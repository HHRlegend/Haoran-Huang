library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
setwd("./")
img2L = Seurat::Read10X_Image('HCC-2L/spatial', image.name = 'tissue_lowres_image.png')
spatial2L = Seurat::Load10X_Spatial(data.dir='HCC-2L',
                                filename = 'filtered_feature_bc_matrix.h5',
                                assay='Spatial',
                                slice='slice2L',
                                image = img2L
)
spatial2L <- SCTransform(spatial2L, assay = "Spatial", verbose = FALSE)

spatial2L <- RunPCA(spatial2L, assay = "SCT", verbose = FALSE)
spatial2L <- FindNeighbors(spatial2L, reduction = "pca", dims = 1:30)
spatial2L <- FindClusters(spatial2L, resolution = 1, verbose = FALSE)
spatial2L <- RunUMAP(spatial2L, reduction = "pca", dims = 1:30)

p1 <- DimPlot(spatial2L, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(spatial2L, label = TRUE, label.size = 3, group.by = 'seurat_clusters')
p1 + p2
ggsave(filename = "spacluster2L.pdf", width =10, height =5)

spa2L <- spatial2L
sce_reference <- SCTransform(sce.all, #ncells = 3000, 
verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
anchors2L <- FindTransferAnchors(reference = sce_reference, query =spa2L ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2L <- TransferData(anchorset = anchors2L, refdata = sce_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction =spa2L[["pca"]], dims = 1:30)
spa2L[["predictions"]] <- predictions.assay2L
DefaultAssay(spa2L ) <- "predictions"
SpatialFeaturePlot(spa2L , features = c("T/NK cells", "Malignant cells",'B cells','Plasma cells','Mast cells','Monocytes','Macrophages','Fibroblasts','Dendritic cells','Endothelial cells'), pt.size.factor =3.5,ncol = 5, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spacelltype2L.pdf", width =20, height =8)

SpatialDimPlot(spatial2L, cells.highlight = CellsByIdentities(object = spatial2L,
                                                         idents = c(0, 1, 2,3,4, 6, 5, 7,8,9)), facet.highlight = TRUE, ncol = 3)
spaepi2L <-subset(spatial2L, idents = c(0, 3, 9))
sce_reference1 <- SCTransform(sceepispa, #ncells = 3000, 
verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
anchors2L1 <- FindTransferAnchors(reference = sce_reference1, query =spaepi2L ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2L1 <- TransferData(anchorset = anchors2L1, refdata = sce_reference1$MVI_Scissor, prediction.assay = TRUE,
                                  weight.reduction =spaepi2L[["pca"]], dims = 1:30)
spaepi2L[["predictions"]] <- predictions.assay2L1
DefaultAssay(spaepi2L ) <- "predictions"
SpatialFeaturePlot(spaepi2L , features = c("0", "22",'11'), pt.size.factor =3,ncol = 3, crop = F, alpha = c(0.1, 1))
ggsave(filename = "spamvi2L.pdf", width =15, height =5)

img2T = Seurat::Read10X_Image('HCC-2T/spatial', image.name = 'tissue_lowres_image.png')
spatial2T= Seurat::Load10X_Spatial(data.dir='HCC-2T',
                                filename = 'filtered_feature_bc_matrix.h5',
                                assay='Spatial',
                                slice='slice2T',
                                image = img2T
)
spatial2T<- SCTransform(spatial2T, assay = "Spatial", verbose = FALSE)

spatial2T <- RunPCA(spatial2T, assay = "SCT", verbose = FALSE)
spatial2T<- FindNeighbors(spatial2T, reduction = "pca", dims = 1:30)
spatial2T<- FindClusters(spatial2T, resolution = 1, verbose = FALSE)
spatial2T <- RunUMAP(spatial2T, reduction = "pca", dims = 1:30)

p3 <- DimPlot(spatial2T, reduction = "umap", label = TRUE)
p4 <- SpatialDimPlot(spatial2T, label = TRUE, label.size = 3, group.by = 'seurat_clusters')
p3 + p4
ggsave(filename = "spacluster2T.pdf", width =10, height =5)

spa2T <- spatial2T
anchors2T <- FindTransferAnchors(reference = sce_reference, query =spa2T ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2T <- TransferData(anchorset = anchors2T, refdata = sce_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction =spa2T[["pca"]], dims = 1:30)
spa2T[["predictions"]] <- predictions.assay2T
DefaultAssay(spa2T ) <- "predictions"
write.csv(spa2T@assays[["predictions"]]@data , "predictions2T.csv")
celltype_2T=read.table('celltype2T.txt',header = F,row.names = 1,sep = '\t')
spa2T<- AddMetaData(spa2T, metadata = celltype_2T, col.name = "celltype_2T")
SpatialFeaturePlot(spa2T , features = c("T/NK cells", "Malignant cells",'B cells','Plasma cells','Mast cells','Monocytes','Macrophages','Fibroblasts','Dendritic cells','Endothelial cells'), pt.size.factor =3.5,ncol = 5, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spacelltype2T.pdf", width =20, height =8)

SpatialDimPlot(spatial2T, cells.highlight = CellsByIdentities(object = spatial2T,
                                                         idents = c(0, 1, 2,3,4, 6, 5, 7,8,9,10,11,12,13,14)), facet.highlight = TRUE, ncol = 3)
spaepi2T <- subset(spatial2T, idents = c(1, 2, 3, 5,6,8,9,11,12,13))
anchors2T1 <- FindTransferAnchors(reference = sce_reference1, query =spaepi2T ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2T1 <- TransferData(anchorset = anchors2T1, refdata = sce_reference1$MVI_Scissor, prediction.assay = TRUE,
                                  weight.reduction =spaepi2T[["pca"]], dims = 1:30)
spaepi2T[["predictions"]] <- predictions.assay2T1
DefaultAssay(spaepi2T ) <- "predictions"
write.csv(spaepi2T@assays[["predictions"]]@data , "predictions2TEPI.csv")
SpatialFeaturePlot(spaepi2T , features = c("0", "22",'11'), pt.size.factor =3,ncol = 3, crop =T, alpha = c(0.1, 1))
ggsave(filename = "spamvi2T.pdf", width =15, height =5)

img2P = Seurat::Read10X_Image('HCC-2P/spatial', image.name = 'tissue_lowres_image.png')
spatial2P = Seurat::Load10X_Spatial(data.dir='HCC-2P',
                                filename = 'filtered_feature_bc_matrix.h5',
                                assay='Spatial',
                                slice='slice2P',
                                image = img2P
)
spatial2P<- SCTransform(spatial2P, assay = "Spatial", verbose = FALSE)

spatial2P <- RunPCA(spatial2P, assay = "SCT", verbose = FALSE)
spatial2P<- FindNeighbors(spatial2P, reduction = "pca", dims = 1:30)
spatial2P<- FindClusters(spatial2P, resolution = 1, verbose = FALSE)
spatial2P <- RunUMAP(spatial2P, reduction = "pca", dims = 1:30)

p5 <- DimPlot(spatial2P, reduction = "umap", label = TRUE)
p6<- SpatialDimPlot(spatial2P, label = TRUE, label.size = 3, group.by = 'seurat_clusters')
p5 + p6
ggsave(filename = "spacluster2P.pdf", width =10, height =5)

spa2P <- spatial2P
anchors2P <- FindTransferAnchors(reference = sce_reference, query =spa2P ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2P <- TransferData(anchorset = anchors2P, refdata = sce_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction =spa2P[["pca"]], dims = 1:30)
spa2P[["predictions"]] <- predictions.assay2P
DefaultAssay(spa2P ) <- "predictions"
write.csv(spa2P@assays[["predictions"]]@data , "predictions2P.csv")
celltype_2P=read.table('celltype2P.txt',header = F,row.names = 1,sep = '\t')
spa2P<- AddMetaData(spa2P, metadata = celltype_2P, col.name = "celltype_2P")
SpatialFeaturePlot(spa2P , features = c("T/NK cells", "Malignant cells",'B cells','Plasma cells','Mast cells','Monocytes','Macrophages','Fibroblasts','Dendritic cells','Endothelial cells'), pt.size.factor =3.5,ncol = 5, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spacelltype2P.pdf", width =20, height =8)

SpatialDimPlot(spatial2P, cells.highlight = CellsByIdentities(object = spatial2P,
                                                         idents = c(0, 1, 2,3,4, 6, 5, 7,8,9,10)), facet.highlight = TRUE, ncol = 3)
spaepi2P <- subset(spatial2P, idents = c(4, 5, 6, 9))
anchors2P1 <- FindTransferAnchors(reference = sce_reference1, query =spaepi2P ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2P1 <- TransferData(anchorset = anchors2P1, refdata = sce_reference1$MVI_Scissor, prediction.assay = TRUE,
                                  weight.reduction =spaepi2P[["pca"]], dims = 1:30)
spaepi2P[["predictions"]] <- predictions.assay2P1
DefaultAssay(spaepi2P ) <- "predictions"
SpatialFeaturePlot(spaepi2P , features = c("0", "22",'11'), pt.size.factor =3,ncol = 3, crop =T, alpha = c(0.1, 1))
ggsave(filename = "spamvi2P.pdf", width =15, height =5)

sce_referencecellchat <- SCTransform(scecellchat, #ncells = 3000, 
verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
anchors2Tcellchat <- FindTransferAnchors(reference = sce_referencecellchat, query =spa2T,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay2Tcellchat <- TransferData(anchorset = anchors2Tcellchat, refdata = sce_referencecellchat$celltypecellchat, prediction.assay = TRUE,
                                  weight.reduction =spa2T[["pca"]], dims = 1:30)
spa2T[["predictionscellchat"]] <- predictions.assay2Tcellchat
DefaultAssay(spa2T) <- "predictionscellchat"

data.input = Seurat::GetAssayData(spa2T, slot = "data", assay = "SCT") 

meta =  celltype_2T
unique(meta$V2)

spatial.locs = Seurat::GetTissueCoordinates(spa2T, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 

scale.factors = jsonlite::fromJSON(txt = 
                                   file.path("./", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)
cellchatST <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "V2", 
                           datatype = "spatial", ###
                           coordinates = spatial.locs, 
                           scale.factors = scale.factors)
CellChatDB <- CellChatDB.human
cellchatST@DB=CellChatDB
cellchatST <- subsetData(cellchatST) 
future::plan("multisession", workers = 8) 
cellchatST <- identifyOverExpressedGenes(cellchatST)
cellchatST <- identifyOverExpressedInteractions(cellchatST)

cellchatST <- computeCommunProb(cellchatST, 
                              type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, 
                              scale.distance = 0.01)
cellchatST <- filterCommunication(cellchatST, min.cells = 10)
cellchatST <- computeCommunProbPathway(cellchatST)

cellchatST <- aggregateNet(cellchatST)
cellchatST@netP$pathways
df.netST <- subsetCommunication(cellchatST)
df.pathwayST <- subsetCommunication(cellchatST,slot.name = 'netP')
write.csv(df.netST, "cell-cell_communications.allST.csv")
write.csv(df.pathwayST , "cell-cell_communicationspath.allST.csv")
pathways.show <- c("MIF")
levels(cellchatST@idents)  
vertex.receiver = c(3,2,4,5,8,9)
netVisual_aggregate(cellchatST, signaling = pathways.show,                      
vertex.receiver = vertex.receiver,layout = "hierarchy")
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchatST, signaling = pathways.show, layout = "spatial", 
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)


