library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
setwd("./")
img1L = Seurat::Read10X_Image('HCC-1L/spatial', image.name = 'tissue_lowres_image.png')
spatial1L = Seurat::Load10X_Spatial(data.dir='HCC-1L',
                                filename = 'filtered_feature_bc_matrix.h5',
                                assay='Spatial',
                                slice='slice1L',
                                image = img1L
)
spatial1L <- SCTransform(spatial1L, assay = "Spatial", verbose = FALSE)

spatial1L <- RunPCA(spatial1L, assay = "SCT", verbose = FALSE)
spatial1L <- FindNeighbors(spatial1L, reduction = "pca", dims = 1:30)
spatial1L <- FindClusters(spatial1L, resolution = 1, verbose = FALSE)
spatial1L <- RunUMAP(spatial1L, reduction = "pca", dims = 1:30)

p7 <- DimPlot(spatial1L, reduction = "umap", label = TRUE)
p8 <- SpatialDimPlot(spatial1L, label = TRUE, label.size = 3, group.by = 'seurat_clusters')
p7+ p8
ggsave(filename = "spacluster1L.pdf", width =10, height =5)

sce_reference <- SCTransform(sce.all, 
verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
sce_reference1 <- SCTransform(sceepi, 
verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

spa1L <- spatial1L
anchors1L <- FindTransferAnchors(reference = sce_reference, query =spa1L ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay1L <- TransferData(anchorset = anchors1L, refdata = sce_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction =spa1L[["pca"]], dims = 1:30)
spa1L[["predictions"]] <- predictions.assay1L
DefaultAssay(spa1L ) <- "predictions"
SpatialFeaturePlot(spa1L , features = c("T/NK cells", "Malignant cells",'B cells','Plasma cells','Mast cells','Monocytes','Macrophages','Fibroblasts','Dendritic cells','Endothelial cells'),  pt.size.factor =3.5, ncol = 5, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spacelltype1L.pdf", width =20, height =8)

SpatialDimPlot(spatial1L, cells.highlight = CellsByIdentities(object = spatial1L,
                                                         idents = c(0, 1, 2,3,4, 6, 5, 7,8,9,10,11)), facet.highlight = TRUE, ncol = 3)
spaepi1L <-subset(spatial1L, idents = c(0, 2, 4))
anchors1L1 <- FindTransferAnchors(reference = sce_reference1, query =spaepi1L ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay1L1 <- TransferData(anchorset = anchors1L1, refdata = sce_reference1$MVI_Scissor, prediction.assay = TRUE,
                                  weight.reduction =spaepi1L[["pca"]], dims = 1:30)
spaepi1L[["predictions"]] <- predictions.assay1L1
DefaultAssay(spaepi1L ) <- "predictions"
SpatialFeaturePlot(spaepi1L , features = c("0", "22",'11'),  pt.size.factor =4, ncol = 3, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spamvi1L.pdf", width =15, height =5)

img1T = Seurat::Read10X_Image('HCC-1T/spatial', image.name = 'tissue_lowres_image.png')
spatial1T= Seurat::Load10X_Spatial(data.dir='HCC-1T',
                                filename = 'filtered_feature_bc_matrix.h5',
                                assay='Spatial',
                                slice='slice1T',
                                image = img1T
)
spatial1T<- SCTransform(spatial1T, assay = "Spatial", verbose = FALSE)

spatial1T <- RunPCA(spatial1T, assay = "SCT", verbose = FALSE)
spatial1T<- FindNeighbors(spatial1T, reduction = "pca", dims = 1:30)
spatial1T<- FindClusters(spatial1T, resolution = 1, verbose = FALSE)
spatial1T <- RunUMAP(spatial1T, reduction = "pca", dims = 1:30)

p9 <- DimPlot(spatial1T, reduction = "umap", label = TRUE)
p10 <- SpatialDimPlot(spatial1T, label = TRUE, label.size = 3, group.by = 'seurat_clusters')
p9 + p10
ggsave(filename = "spacluster1T.pdf", width =10, height =5)

spa1T <- spatial1T
anchors1T <- FindTransferAnchors(reference = sce_reference, query =spa1T ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay1T <- TransferData(anchorset = anchors1T, refdata = sce_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction =spa1T[["pca"]], dims = 1:30)
spa1T[["predictions"]] <- predictions.assay1T
DefaultAssay(spa1T ) <- "predictions"
SpatialFeaturePlot(spa1T , features = c("T/NK cells", "Malignant cells",'B cells','Plasma cells','Mast cells','Monocytes','Macrophages','Fibroblasts','Dendritic cells','Endothelial cells'), pt.size.factor =3.5, ncol = 5, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spacelltype1T.pdf", width =20, height =8)

SpatialDimPlot(spatial1T, cells.highlight = CellsByIdentities(object = spatial1T,
                                                         idents = c(0, 1, 2,3,4, 6, 5, 7,8,9)), facet.highlight = TRUE, ncol = 3)
spaepi1T <-subset(spatial1T, idents = c(0,1,2,6,7, 8))
anchors1T1 <- FindTransferAnchors(reference = sce_reference1, query =spaepi1T ,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay1T1 <- TransferData(anchorset = anchors1T1, refdata = sce_reference1$MVI_Scissor, prediction.assay = TRUE,
                                  weight.reduction =spaepi1T[["pca"]], dims = 1:30)
spaepi1T[["predictions"]] <- predictions.assay1T1
DefaultAssay(spaepi1T ) <- "predictions"
SpatialFeaturePlot(spaepi1T , features = c("0", "22",'11'),  pt.size.factor =3.5,ncol = 3, crop = T, alpha = c(0.1, 1))
ggsave(filename = "spamvi1T.pdf", width =15, height =5)






