## Visium for cell type annotation based on scRNA-seq clustering
s_df = readRDS("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/s_df_filtered.rds")
s_df_for_visium = s_df
Idents(s_df_for_visium) <- s_df_for_visium$cell_types
#For speed, and for a more fair comparison of the celltypes, we will subsample all celltypes to a maximum of 200 cells per class (subclass)
s_df_for_visium <- subset(s_df_for_visium, cells = WhichCells(s_df_for_visium, downsample = 200))
s_df_for_visium <- SCTransform(s_df_for_visium, ncells = 3000, verbose = FALSE, method = "poisson") %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:15)
DimPlot(s_df_for_visium, label = TRUE)

#54A
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/54A/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U054", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE, resolution = 0)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
DimPlot(visium, reduction = "umap", group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(visium)
visium_sub <- subset(visium, idents = c(0))
SpatialDimPlot(visium_sub, crop = TRUE, label = F, cols = "orchid2")

anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE, weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
visium_54A = visium

#54B
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/54B/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U054", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))
anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE, weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
visium_54B = visium

#55A
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/55A/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U055", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))
anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
visium_55A = visium

#55B
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/55B/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U055", filter.matrix = T, to.upper = F)
# store mitochondrial percentage in object meta data
#visium[["percent.mito"]] = PercentageFeatureSet(visium, pattern = "^mt-")
#visium[["percent.ribo"]] = PercentageFeatureSet(visium, pattern = "^Rp[sl]")
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))
anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
visium_55B = visium


#55C2
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/55C2/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U055", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)

## visium region demarcation based on Cd68/Adgre1
p = SpatialFeaturePlot(object = visium, features = "Grn", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 2.2)
p
p = SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 0.6)
p
p = SpatialFeaturePlot(object = visium, features = "Egfr", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 1.5)
p
p = SpatialFeaturePlot(object = visium, features = "Adgre1", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 0.7)
p

p + ggplot2::scale_fill_continuous(limits = c(0,2), breaks = c(0, 1, 2), type = "viridis")
p + ggplot2::scale_colour_gradient(low = "white", high = "black")
p + ggplot2::scale_fill_continuous(limits = c(0,2), breaks = c(0, 1, 2), type = "gradient")

visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE, resolution = 0.2)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
DimPlot(visium, reduction = "umap", group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(visium)
visium_sub <- subset(visium, idents = c(3))
SpatialDimPlot(visium_sub, crop = TRUE, label = F, cols = "orchid2")

anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
SpatialFeaturePlot(object = visium, features = "Hepatocyte", pt.size.factor = 2, alpha = c(0.1, 1))
SpatialFeaturePlot(object = visium, features = "Macrophage", pt.size.factor = 2, alpha = c(0.1, 1))
visium_55C2 = visium


#57A
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/57A/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U057", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE, resolution = 0)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
DimPlot(visium, reduction = "umap", group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(visium)
visium_sub <- subset(visium, idents = c(0))
SpatialDimPlot(visium_sub, crop = TRUE, label = F, cols = "orchid2")

anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
SpatialFeaturePlot(object = visium, features = "Hepatocyte", pt.size.factor = 2, alpha = c(0.1, 1))
visium_57A = visium

#57B
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/57B/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U057", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE, resolution = 0)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
DimPlot(visium, reduction = "umap", group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(visium)
visium_sub <- subset(visium, idents = c(0))
SpatialDimPlot(visium_sub, crop = TRUE, label = F, cols = "orchid2")

anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
SpatialFeaturePlot(object = visium, features = "Hepatocyte", pt.size.factor = 2, alpha = c(0.1, 1))
visium_57B = visium


#57C1
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/57C1/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U057", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)


#DefaultAssay(visium) <- "SCT"
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))

## visium region demarcation based on Cd68/Adgre1
p = SpatialFeaturePlot(object = visium, features = "Grn", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 2)
p
p = SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 0.6)
p
p = SpatialFeaturePlot(object = visium, features = "Egfr", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 1.5)
p
p = SpatialFeaturePlot(object = visium, features = "Adgre1", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 0.7)
p


visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE, resolution = 0.3)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
DimPlot(visium, reduction = "umap", group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(visium)
visium_sub <- subset(visium, idents = c(3))
SpatialDimPlot(visium_sub, crop = TRUE, label = F, cols = "orchid2")

anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% dplyr::select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
SpatialFeaturePlot(object = visium, features = "Hepatocyte", pt.size.factor = 2, alpha = c(0.1, 1))
visium_57C1 = visium


#57C2
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/57C2/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U057", filter.matrix = T, to.upper = F)
visium <- PercentageFeatureSet(visium, "^mt-", col.name = "percent_mito")
visium <- PercentageFeatureSet(visium, "^Hb.*-", col.name = "percent_hb")
#filter spots
visium = visium[, visium$nFeature_Spatial > 500 & visium$percent_mito <= 10 & visium$percent_hb <= 10]

#filter genes
# Filter Mitocondrial
visium <- visium[!grepl("^mt-", rownames(visium)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
visium <- visium[!grepl("^Hb.*-", rownames(visium)), ]

# Filter genes starting with "Gm"
#visium <- visium[!grepl("^Gm", rownames(visium)), ]

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.1, 1))


## visium region demarcation based on Cd68/Adgre1
p = SpatialFeaturePlot(object = visium, features = "Grn", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 2.2)
p
p = SpatialFeaturePlot(object = visium, features = "Cd68", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 0.7)
p
p = SpatialFeaturePlot(object = visium, features = "Egfr", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 1.5)
p
p = SpatialFeaturePlot(object = visium, features = "Adgre1", pt.size.factor = 2, alpha = c(0.5, 1), min.cutoff = 0.7)
p

#visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE, resolution = 0.4)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
DimPlot(visium, reduction = "umap", group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(visium)
visium_sub <- subset(visium, idents = c(3))
SpatialDimPlot(visium_sub, crop = TRUE, label = F, cols = "orchid2")

anchors <- FindTransferAnchors(reference = s_df_for_visium, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = s_df_for_visium$cell_types, prediction.assay = TRUE, weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#SpatialFeaturePlot(visium, features = c("Hepatocyte", "CD3+ T Cell", "Endothelial Cell", "M2 Macrophage", "B Cell", "Macrophage", "Germinal Center B Cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, combine = F)

visium <- FindSpatiallyVariableFeatures(visium, assay = "predictions", selection.method = "markvariogram",
                                        features = rownames(visium), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(visium), 4)
SpatialPlot(object = visium, features = top.clusters, ncol = 2)


df = visium@assays$predictions[1,] %>% t()
df2 = visium@assays$predictions@data %>% as.data.frame()
df_hepatocytes = df[rowSums(df) > 0.7, ] %>% as.data.frame() %>% rownames()
df_remaining = df2 %>% .[2:9, which(!(names(df2) %in% df_hepatocytes))] %>% as.data.frame()


#Let soln be the solution vector that stores the corresponding maximum value of each column              
soln=c()
column_names = c()
#Traverse the matrix column-wise
for (i in 1:ncol(df_remaining))
{
  #Extract all rows of the ith column and find the maxiumum value in the same column
  soln[i] = rownames(df_remaining) %>% .[which.max(df_remaining[,i])]
  column_names[i] = colnames(df_remaining)[i]
} 
tmp = data.frame(df_hepatocytes, "Hepatocyte")
colnames(tmp) = c("column_names", "soln")
df_merged = rbind(data.frame(column_names, soln), tmp)
df_without_rownames = tibble::rownames_to_column(as.data.frame(df), var = "column_names")
df_final = left_join(df_without_rownames, df_merged, by = "column_names") %>% dplyr::select(-Hepatocyte)
visium@meta.data$cell_type_scRNAseq = df_final$soln
SpatialDimPlot(visium, label = F, label.size = 3, group.by = "cell_type_scRNAseq", pt.size.factor = 1)
p = SpatialFeaturePlot(object = visium, features = "Hepatocyte", pt.size.factor = 2, alpha = c(0.1, 1))
p + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")
p = SpatialFeaturePlot(object = visium, features = "Macrophage", pt.size.factor = 2, alpha = c(0.1, 1))
p + ggplot2::scale_fill_continuous(limits = c(0.0,0.5), breaks = c(0.0, 0.25, 0.5), type = "viridis")

##1## ISG score for visium (deprecated, please use ##3## for this analysis)
visium_55A@meta.data$orig.ident = "55A"
visium_57B@meta.data$orig.ident = "57B"
visium_57C1@meta.data$orig.ident = "57C1"
samples = c("visium_54A","visium_54B","visium_55A","visium_55B","visium_55C2","visium_57A","visium_57B","visium_57C1","visium_57C2")
for (val in samples) {
  vals = get(val)
  genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
  vals = ScaleData(vals, assay = "Spatial", features = genes_isg, do.center = T, do.scale = T, model.use = "poisson")
  vals@meta.data$ISG_score = colMeans(vals[["Spatial"]]@scale.data)
  p0 = SpatialFeaturePlot(object = vals, features = "ISG_score", pt.size.factor = 2, alpha = c(0.1, 1))
  p = p0 + ggplot2::scale_fill_continuous(limits = c(0,2), breaks = c(0, 1, 2), type = "viridis")
  ggsave(filename = paste0("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/visium_ISG_plot.", val, ".pdf"), plot = p, width = 6, height = 6)
}

##2## Basic box plot
visium_55A@meta.data$orig.ident = "55A"
visium_57B@meta.data$orig.ident = "57B"
visium_57C1@meta.data$orig.ident = "57C1"
library(ggpubr)
genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
visium_55A = ScaleData(visium_55A, assay = "SCT", features = genes_isg, do.center = F, do.scale = T)
visium_55A@meta.data$ISG_score = colMeans(visium_55A[["SCT"]]@scale.data)
df1 = visium_55A@meta.data %>% select(orig.ident, ISG_score) %>% mutate(log_ISG= log10(ISG_score + 1)) 
visium_57B = ScaleData(visium_57B, assay = "SCT", features = genes_isg, do.center = F, do.scale = T)
visium_57B@meta.data$ISG_score = colMeans(visium_57B[["SCT"]]@scale.data)
df2 = visium_57B@meta.data %>% select(orig.ident, ISG_score) %>% mutate(log_ISG= log10(ISG_score + 1)) 
visium_57C1 = ScaleData(visium_57C1, assay = "SCT", features = genes_isg, do.center = F, do.scale = T)
visium_57C1@meta.data$ISG_score = colMeans(visium_57C1[["SCT"]]@scale.data)
df3 = visium_57C1@meta.data %>% select(orig.ident, ISG_score) %>% mutate(log_ISG= log10(ISG_score + 1)) 
df = rbind(df1,df2,df3)
my_comparisons <- list( c("55A", "57B"), c("57B", "57C1"), c("55A", "57C1"))
ggboxplot(df, x = "orig.ident", y = "ISG_score", #ylim = c(0,2),
          color = "orig.ident", palette = c("#868686FF", "#EFC000FF", "#0073C2FF"))+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')

boxplot(log_ISG ~ orig.ident, data = df, #ylim=c(-0.5,0.6), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#868686FF", "#EFC000FF", "#0073C2FF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'ISG score of the three samples')


##3## Plot visium Spatial featureplot for ISG score:
samples = c("visium_55A","visium_57B","visium_57C1")
genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
visium = ScaleData(visium_55A, assay = "SCT", features = genes_isg, do.center = F, do.scale = T)
visium@meta.data$ISG_score = colMeans(visium[["SCT"]]@scale.data)
#visium = visium %>% subset(visium %in% "Hepatocyte") 
p1 = SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.4, 1))
p1 + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")
visium = ScaleData(visium_57B, assay = "SCT", features = genes_isg, do.center = F, do.scale = T)
visium@meta.data$ISG_score = colMeans(visium[["SCT"]]@scale.data)
#visium = visium %>% subset(visium %in% "Hepatocyte") 
p1 = SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.4, 1))
p1 + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")
visium = ScaleData(visium_57C1, assay = "SCT", features = genes_isg, do.center = F, do.scale = T)
visium@meta.data$ISG_score = colMeans(visium[["SCT"]]@scale.data)
#visium = visium %>% subset(visium %in% "Hepatocyte") 
p1 = SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.4, 1))
p1 + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")

visium_scRNAseq_celltrek_57C1 = SetIdent(visium_scRNAseq_celltrek_57C1, value = visium_scRNAseq_celltrek_57C1[["cell_types"]])
p = SpatialFeaturePlot(object = visium_scRNAseq_celltrek_57C1, features = "Hepatocyte", pt.size.factor = 2, alpha = c(0.1, 1))
p + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")
p = SpatialFeaturePlot(object = visium_scRNAseq_celltrek_57C1, features = "Macrophage", pt.size.factor = 2, alpha = c(0.1, 1))
p + ggplot2::scale_fill_continuous(limits = c(0.0,0.5), breaks = c(0.0, 0.25, 0.5), type = "viridis")

##4#### Plot two features simutaneously in Visium ST
test = visium_55A@images$V10U055@coordinates[,2:3] %>% as.matrix()
colnames(test) = c("PC_1", "PC_2")
visium_55A@reductions$pca@cell.embeddings = test
#FeaturePlot(object = visium2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 1)
p = FeaturePlot(object = visium_55A, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 0)
suppressMessages(ggsave(paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/Grn_Egfr_co-expression_55A.light_grey_background.pdf"), plot = p, width = 17.95, height = 4.7))

test = visium_57B@images$V10U057@coordinates[,2:3] %>% as.matrix()
colnames(test) = c("PC_1", "PC_2")
visium_57B@reductions$pca@cell.embeddings = test
#FeaturePlot(object = visium2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 1)
p = FeaturePlot(object = visium_57B, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 0)
suppressMessages(ggsave(paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/Grn_Egfr_co-expression_57B.light_grey_background.pdf"), plot = p, width = 17.95, height = 4.7))

test = visium_57C1@images$V10U057@coordinates[,2:3] %>% as.matrix()
colnames(test) = c("PC_1", "PC_2")
visium_57C1@reductions$pca@cell.embeddings = test
#FeaturePlot(object = visium2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 1)
p = FeaturePlot(object = visium_57C1, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 0)
suppressMessages(ggsave(paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/Grn_Egfr_co-expression_57C1.light_grey_background.pdf"), plot = p, width = 17.95, height = 4.7))

test = visium_57C2@images$V10U057@coordinates[,2:3] %>% as.matrix()
colnames(test) = c("PC_1", "PC_2")
visium_57C2@reductions$pca@cell.embeddings = test
#FeaturePlot(object = visium2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 1)
p = FeaturePlot(object = visium_57C2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 0)
suppressMessages(ggsave(paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/Grn_Egfr_co-expression_57C2.light_grey_background.pdf"), plot = p, width = 17.95, height = 4.7))

test = visium_55C2@images$V10U055@coordinates[,2:3] %>% as.matrix()
colnames(test) = c("PC_1", "PC_2")
visium_55C2@reductions$pca@cell.embeddings = test
#FeaturePlot(object = visium2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 1)
p = FeaturePlot(object = visium_55C2, features = c('Grn', 'Egfr'), blend = TRUE, cols = c("light grey", "orange", "blue"), pt.size = 1.5, order = T, min.cutoff = 0)
suppressMessages(ggsave(paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/Grn_Egfr_co-expression_55C2.light_grey_background.pdf"), plot = p, width = 17.95, height = 4.7))




visium_57C2 = visium

visium_54A@meta.data$orig.ident = "sampleA"
visium_55A@meta.data$orig.ident = "sampleA"
visium_57A@meta.data$orig.ident = "sampleA"
visium_54B@meta.data$orig.ident = "sampleB"
visium_55B@meta.data$orig.ident = "sampleB"
visium_57B@meta.data$orig.ident = "sampleB"
visium_55C2@meta.data$orig.ident = "sampleC"
visium_57C1@meta.data$orig.ident = "sampleC"
visium_57C2@meta.data$orig.ident = "sampleC"
#visium_merged <- merge(merge(merge(merge(merge(visium_54A, visium_55A), visium_55C2), visium_57A), visium_57C1), visium_57C2)
visium_merged <- merge(merge(merge(merge(merge(merge(merge(merge(visium_54A, visium_54B), visium_55A), visium_55B), visium_55C2), visium_57A), visium_57B), visium_57C1), visium_57_C2)

##3## Find the DEGs between samples within one cluster
# Save current cluster identites in object@meta.data under 'cell_types'
# Run only if Seurat::FindClusters() was executed
s_df_new = SetIdent(visium_merged, value = visium_merged[["cell_type_scRNAseq"]])
object <- Seurat::StashIdent(object = s_df_new, save.name = "cell_types")
object_sub = subset(object, idents = "Hepatocyte")

# Set cell identity to sample identity
Idents(object_sub = object_sub) <- object_sub@meta.data$orig.ident

# Find all sample specific marker genes 
#markers <- Seurat::FindMarkers(object = object,ident.1 = "Hepatocyte_C_Hom", ident.2 = c("Hepatocyte_A_WT", "Hepatocyte_B_Het"))

## Do the above command in loop:
#for each cluster name: load file and create outline
markers <- Seurat::FindMarkers(object = object_sub,ident.1 = "sampleC", ident.2 = "sampleA", assay = "Spatial")
markers = markers %>% tibble::rownames_to_column(var = "gene")
write.delim(markers, "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/DEG_clusters_from_C_compared_to_A.scCATCH_annotated.visium_only_hepatocyte.txt")
















samples = c("visium_54A","visium_54B","visium_55A","visium_55B","visium_55C2","visium_57A","visium_57B","visium_57C1","visium_57C2")
for (val in samples) {
## get the DEG genes based on ISG score (ISG high vs ISG low) of hepatocytes using Visium spatial transcriptomics
genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
vals = get(val)
visium = ScaleData(vals, assay = "Spatial", features = genes_isg, do.center = T, do.scale = T)
visium@meta.data$ISG_score = colMeans(visium[["Spatial"]]@scale.data)
ISG = visium@meta.data %>% subset(cell_type_scRNAseq %in% "Hepatocyte") %>% select(ISG_score) %>% t() %>% as.data.frame()
quantiles = ISG %>% quantile(seq(0, 1, 1/4)) %>% as.numeric()
high_aids = ISG[which(ISG <= quantiles[5] & ISG >= quantiles[4])] %>% colnames()
low_aids = ISG[which(ISG <= quantiles[2] & ISG >= quantiles[1])] %>% colnames()
markers = FindMarkers(visium, ident.1 = high_aids, ident.2 = low_aids, assay = "Spatial", min.pct = 0.1)
markers = markers %>% tibble::rownames_to_column(var = "gene")
markers %>% write.delim(paste0("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/DE_genes_based_on_ISG_score_of_hepatocytes_minpct_0.1_using_Visium.",val,".txt"))
}














## DEG using DESeq2
library(DESeq2)
cts = visium@assays$Spatial@counts %>% as.data.frame()
high = data.frame(name=high_aids)
high$condition = "high"
low = data.frame(name=low_aids)
low$condition = "low"
coldata = rbind(high, low) %>% tibble::column_to_rownames(var = "name")
ids=c(high_aids,low_aids)
cts = cts%>% select(ids)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds, test="LRT", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
resultsNames(dds) # lists the coefficients
res <- results(dds)
res = res[complete.cases(res), ] %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene")

write.delim(res, "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/DESeq2_genes_based_on_ISG_score_of_hepatocytes.using_Visium.57C2.txt")




visium@meta.data$macrophage_score = as.data.frame(visium@assays$predictions@data)[rownames(as.data.frame(visium@assays$predictions@data)) %in% "Macrophage",] %>% t()
SpatialPlot(object = visium, features = "macrophage_score", ncol = 1)
