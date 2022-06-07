options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")

## 1. read visium ST data -- 57C2 as an example
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

CellTrek::celltrek_vis(visium_scRNAseq_celltrek_57C2@meta.data %>% dplyr::select(coord_x, coord_y, cell_types:id_new), visium_scRNAseq_celltrek_57C2@images$V10U057@image, visium_scRNAseq_celltrek_57C2@images$V10U057@scale.factors$lowres)
