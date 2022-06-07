options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")

## 1. read visium ST data -- 57C2
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

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE, method = "poisson") %>%
  RunPCA(assay = "SCT",verbose = FALSE)
#visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
visium@meta.data$id = rownames(visium@meta.data)
s_df@meta.data$id = rownames(s_df@meta.data)
SpatialDimPlot(visium, group.by = c("seurat_clusters"), label = F)
## Visualize the ST data
#SpatialDimPlot(visium)

## Rename the cells/spots with syntactically valid names
visium <- RenameCells(visium, new.names=make.names(Cells(visium)))
#s_df_for_visium <- RenameCells(s_df_for_visium, new.names=make.names(Cells(s_df_for_visium)))


## 2. Visualize the scRNA-seq data
s_df_for_visium = readRDS("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.rds")
s_df_for_visium <- RenameCells(s_df_for_visium, new.names=make.names(Cells(s_df_for_visium)))


## 3. Cell charting using CellTrek
visium_traint <- CellTrek::traint(st_data=visium, sc_data=s_df_for_visium, sc_assay='RNA', cell_names='cell_types')

## We can check the co-embedding result to see if there is overlap between these two data modalities
DimPlot(visium_traint, group.by = "type")

## 4. Chart single cells to their spatial locations
visium_scRNAseq_celltrek <- CellTrek::celltrek(st_sc_int=visium_traint, int_assay='traint', sc_data=s_df_for_visium, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000,  dist_thresh=0.55, top_spot=5, spot_n=10, repel_r=20, repel_iter=20, keep_model=T)$celltrek

## 5. Interactively visualize the CellTrek result using celltrek_vis
visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_54A.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_54A = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_55A.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_55A = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57A.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57A = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_54B.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_54B = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_55B.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_55B = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57B.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57B = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_55C2.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_55C2 = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57C1.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57C1 = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57C2.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57C2 = visium_scRNAseq_celltrek@meta.data %>% as.data.frame()

id_list=list(visium_scRNAseq_celltrek_54A,visium_scRNAseq_celltrek_55A,visium_scRNAseq_celltrek_57A,visium_scRNAseq_celltrek_54B,visium_scRNAseq_celltrek_55B,visium_scRNAseq_celltrek_57B,visium_scRNAseq_celltrek_55C2,visium_scRNAseq_celltrek_57C1,visium_scRNAseq_celltrek_57C2)
cor_counts = do.call("cbind", lapply(id_list, function(id) {
  #id = as.data.frame(id)
  #colnames(id)
  list1=table(id[['cell_types']])
  list2=length(id[['cell_types']])
  #list1
  mapply(FUN = `/`, list1, list2, SIMPLIFY = FALSE)
}))
write.table(cor_counts, "C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/cell_type_population_by_cellTrek.txt", row.names = T, col.names = T, sep = '\t')

CellTrek::celltrek_vis(visium_scRNAseq_celltrek_54A@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek_54A@images$V10U054@image, visium_scRNAseq_celltrek_54A@images$V10U054@scale.factors$lowres)

##6## Visualization of ISG score based on CellTrek result
visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_54A.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_54A = visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_55A.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_55A = visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57A.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57A = visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_54B.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_54B = visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_55B.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_55B= visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57B.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57B= visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_55C2.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_55C2= visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57C1.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57C1 = visium_scRNAseq_celltrek

visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_57C2.rds")
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
visium_scRNAseq_celltrek_57C2 = visium_scRNAseq_celltrek


##Plot visium Spatial featureplot for ISG score:
samples = c("visium_55A","visium_57B","visium_57C1")
for (val in samples) {
  ## get the DEG genes based on ISG score (ISG high vs ISG low) of hepatocytes using Visium spatial transcriptomics
  genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
  vals = get(val)
  visium = ScaleData(visium_55A, assay = "Spatial", features = genes_isg, do.center = F, do.scale = T)
  visium@meta.data$ISG_score = colMeans(visium[["Spatial"]]@scale.data)
  #visium = visium@meta.data %>% subset(visium %in% "Hepatocyte") 
  plot=SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.1, 1))
  plot %>% write.delim(paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/visium_ISG_plot_",val,".mouse_scRNAseq.txt"))
}

genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
visium = ScaleData(visium_scRNAseq_celltrek_55A, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
visium = visium %>% subset(cell_types %in% "Hepatocyte") 
visium@meta.data$ISG_score = colMeans(visium[["RNA"]]@scale.data)
p1 = SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.01, 1))
p1 + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")

visium = ScaleData(visium_scRNAseq_celltrek_57B, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
visium = visium %>% subset(cell_types %in% "Hepatocyte") 
visium@meta.data$ISG_score = colMeans(visium[["RNA"]]@scale.data)
p1 = SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.01, 1))
p1 + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")

visium = ScaleData(visium_scRNAseq_celltrek_57C1, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
visium = visium %>% subset(cell_types %in% "Hepatocyte") 
visium@meta.data$ISG_score = colMeans(visium[["RNA"]]@scale.data)
p1 = SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.01, 1))
p1 + ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0), type = "viridis")


visium = ScaleData(visium_scRNAseq_celltrek_57C1, assay = "Spatial", features = genes_isg, do.center = F, do.scale = T)
visium@meta.data$ISG_score = colMeans(visium[["Spatial"]]@scale.data)
#visium = visium@meta.data %>% subset(visium %in% "Hepatocyte") 
SpatialFeaturePlot(object = visium, features = "ISG_score", pt.size.factor = 2, alpha = c(0.1, 1))


## 1. read visium ST data -- 54A
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/54A/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U054", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_54A.rds")
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/55A/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U055", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_55A.rds")
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/57A/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U057", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_57A.rds")


## 5. Interactively visualize the CellTrek result using celltrek_vis
visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/visium_scRNAseq_celltrek.55A.rds")
visium_scRNAseq_celltrek = visium_scRNAseq_celltrek_55A
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek@images$V10U055@image, visium_scRNAseq_celltrek@images$V10U055@scale.factors$lowres)


visium_scRNAseq_celltrek = visium_scRNAseq_celltrek_57B
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek@images$V10U057@image, visium_scRNAseq_celltrek@images$V10U057@scale.factors$lowres)


visium_scRNAseq_celltrek = visium_scRNAseq_celltrek_57C1
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek@images$V10U057@image, visium_scRNAseq_celltrek@images$V10U057@scale.factors$lowres)


## 1. read visium ST data -- 54B
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/54B/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U054", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_54B.rds")

visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/55B/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U055", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_55B.rds")

visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/57B/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U057", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_57B.rds")


visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/visium_scRNAseq_celltrek.55B.rds")
## 5. Interactively visualize the CellTrek result using celltrek_vis
visium_scRNAseq_celltrek$celltrek$cell_types <- factor(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek$celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek$celltrek@images$V10U055@image, visium_scRNAseq_celltrek$celltrek@images$V10U055@scale.factors$lowres)


visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/visium_scRNAseq_celltrek.54B.rds")
## 5. Interactively visualize the CellTrek result using celltrek_vis
visium_scRNAseq_celltrek$celltrek$cell_types <- factor(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek$celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek$celltrek@images$V10U054@image, visium_scRNAseq_celltrek$celltrek@images$V10U054@scale.factors$lowres)


visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/visium_scRNAseq_celltrek.57B.rds")
## 5. Interactively visualize the CellTrek result using celltrek_vis
visium_scRNAseq_celltrek$celltrek$cell_types <- factor(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek$celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek$celltrek@images$V10U057@image, visium_scRNAseq_celltrek$celltrek@images$V10U057@scale.factors$lowres)

## 1. read visium ST data -- 55C2
visium = Load10X_Spatial(data.dir = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/output/55C2/", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "V10U055", filter.matrix = T, to.upper = F)
saveRDS(visium, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/visium_object_55C2.rds")

## Then go to renxi@SpongeBob:/media/ext6/renxi/projects/WL_scRNAseq/WL_visium_analysis$ Rscript CellTrek_for_ST_data_annotation.54A.R


visium_scRNAseq_celltrek = readRDS("C:/Users/renxi/Desktop/WL analysis/CellTrek/visium_scRNAseq_celltrek.55C2.rds")
## 5. Interactively visualize the CellTrek result using celltrek_vis
visium_scRNAseq_celltrek$celltrek$cell_types <- factor(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek$celltrek@meta.data$cell_types)))

CellTrek::celltrek_vis(visium_scRNAseq_celltrek$celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek$celltrek@images$V10U055@image, visium_scRNAseq_celltrek$celltrek@images$V10U055@scale.factors$lowres)




##6## spatialfeatureplot one predicted cell type per time
#visium_scRNAseq_celltrek_57C1 = SetIdent(visium_scRNAseq_celltrek_57C1, value = visium_scRNAseq_celltrek_57C1[["cell_types"]])
#marco_cell = WhichCells(object = visium_scRNAseq_celltrek_57C1, idents = "Macrophage")
#SpatialDimPlot(object = visium_scRNAseq_celltrek_57C1, cells.highlight = marco_cell, pt.size.factor = 1.5, alpha = c(0.3, 1), cols.highlight = c("green", "white"), image.alpha = 0)
#SpatialDimPlot(object = visium_scRNAseq_celltrek_57C1, cells.highlight = marco_cell, pt.size.factor = 1.5, alpha = c(0.3, 1), cols.highlight = c("green", "white"), image.alpha = 1)
#plot all cell types in visium
for (id in c("54A", "54B", "55A", "57B", "57A", "55B", "55C2", "57C2", "57C1")){
  celltrek_name = get(paste0("visium_scRNAseq_celltrek_",id))
  visium_scRNAseq_celltrek = readRDS(paste0("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_", id, ".rds"))
  visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
  assign(paste("visium_scRNAseq_celltrek_", id, sep=""), visium_scRNAseq_celltrek) 
  
  ## plot cell type in boolean at Visium ST
  celltrek_name = SetIdent(celltrek_name, value = celltrek_name[["cell_types"]])
  p = SpatialDimPlot(object = celltrek_name, pt.size.factor = 1.2, alpha = c(0.5, 1), image.alpha = 0, cols = c("yellow4", "#00B9E3", "purple", "#00BA38", "#F8766D", "yellow", "springgreen3", "#DB72FB", "#FF61C3"))
  ggsave(filename = paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/visium/All cell types spots predicted by cellTrek",id,".pdf"), plot = p)
}

#borrowing FeaturePlot function
for (id in c("54A", "54B", "55A", "57A", "57B", "55B", "55C2", "57C2", "57C1")){
  celltrek_name = get(paste0("visium_scRNAseq_celltrek_",id))
  slide_no = substr(id, 2, 2)
  image_name = paste0("V10U05",slide_no)
  
visium_scRNAseq_celltrek = readRDS(paste0("C:/Users/renxi/Desktop/WL analysis/CellTrek/newest results/output/visium_scRNAseq_celltrek_", id, ".rds"))
visium_scRNAseq_celltrek$cell_types <- factor(visium_scRNAseq_celltrek@meta.data$cell_types, levels=sort(unique(visium_scRNAseq_celltrek@meta.data$cell_types)))
assign(paste("visium_scRNAseq_celltrek_", id, sep=""), visium_scRNAseq_celltrek) 

## plot cell type in boolean at Visium ST
celltrek_name = SetIdent(celltrek_name, value = celltrek_name[["cell_types"]])
cells.ident <- FetchData(celltrek_name, vars="ident")
hepotacyte_ids = grep("Hepatocyte", cells.ident$ident, ignore.case=T)
macrophage_ids = grep("Macrophage", cells.ident$ident, ignore.case=T)
CD3Tcell_ids = grep("CD3", cells.ident$ident, ignore.case=T)
CD8Tcell_ids = grep("CD8", cells.ident$ident, ignore.case=T)

celltrek_name@reductions$umap <- NULL
test = celltrek_name@reductions$celltrek_raw@cell.embeddings[,1:2] %>% as.matrix()
colnames(test) = c("pca_1", "pca_2")
test = test*celltrek_name@images[[image_name]]@scale.factors$lowres
celltrek_name@reductions$pca@cell.embeddings = test

celltrek_name@meta.data$CD8Tcell_boolean = 0
for(i in 1:nrow(celltrek_name@meta.data)) {
  if(celltrek_name@meta.data[i,8] != "CD8+ T Cell") {
    celltrek_name@meta.data[i,14] = 0
  } else {
    celltrek_name@meta.data[i,14] = 1
  }
}
p = FeaturePlot(object = celltrek_name, features = c('CD8Tcell_boolean'), cols = c("light grey", "#00B9E3"), pt.size = 1.5, order = T, min.cutoff = 0, cells = CD8Tcell_ids)
#p = FeaturePlot(object = celltrek_name, features = c('CD3Tcell_boolean'), cols = c("light grey", "yellow4"), pt.size = 1.5, order = T, min.cutoff = 0, cells = CD8Tcell_ids)
#p = FeaturePlot(object = celltrek_name, features = c('Hepatocyte_boolean'), cols = c("light grey", "#F8766D"), pt.size = 1.5, order = T, min.cutoff = 0, cells = hepotacyte_ids)
#p = FeaturePlot(object = celltrek_name, features = c('Macrophage_boolean'), cols = c("light grey", "springgreen3"), pt.size = 1.5, order = T, min.cutoff = 0, cells = macrophage_ids)
ggsave(filename = paste0("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/visium/CD8+ T cell spots predicted by cellTrek",id,".pdf"), plot = p)
}

##default ggplot2 color scheme for Seurat Dimplot function when ploting different cell types
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n] }
gg_color_hue(9)
#[1] "#F8766D" "#D39200" "#93AA00" "#00BA38" "#00C19F" "#00B9E3" "#619CFF" "#DB72FB" "#FF61C3"