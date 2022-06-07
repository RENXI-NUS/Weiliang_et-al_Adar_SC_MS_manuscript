## Integrating different scRNAseq datasets
library(Seurat)
library(SeuratData)
library(caroline)
library(magrittr)
library(tidyverse)
library(viridis)
library(patchwork)
library(scCATCH)
library(rafalib)
library(clustree)
library(caroline)
library(magrittr)
library(tidyr)
library(dplyr)
library(gplots)
library(ggplot2)

project_name = "mouse_liver"
min_cells = 3 #Include features detected in at least this many cells. Will subset the counts matrix as well. 
min_features = 100 #Include cells where at least this many features are detected.
max_features = 3100
min_nFeature_RNA = 100 #nFeature_RNA is the number of genes detected in each cell.
max_nFeature_RNA = 3100
#min_genes = 100
#max_genes = 3100
min_UMIs = 500
max_UMIs = 10000
percent_mito = 20
percent_ribo = 15
top_n_features = 20
n_dims = 15

s_data_list = Read10X(data.dir = paste0("C:/Users/renxi/Downloads/WL_scRNAseq_batch2_enforcing10k/aggregated_filtered_feature_bc_matrix/"))
s_df <- CreateSeuratObject(counts = s_data_list, project = "mouse_liver")
s_df[["percent.mito"]] = PercentageFeatureSet(s_df, pattern = "^mt-")
s_df[["percent.ribo"]] = PercentageFeatureSet(s_df, pattern = "^Rp[sl]")
s_df = subset(s_df, subset = nFeature_RNA >= min_nFeature_RNA & percent.mito <= percent_mito )
#s_df <- s_df[!grepl("^mt-", rownames(s_df)), ] #When plotting the heatmap of DEG, remove the ^mt genes

orig.ident = str_split_fixed(s_df@meta.data %>% as.data.frame() %>% tibble::rownames_to_column(var = "barcodes") %>% .$barcodes, "-", 2) %>% as.data.frame() %>% .$V2
s_df@meta.data$orig.ident = orig.ident

s_df.list <- SplitObject(s_df, split.by = "orig.ident")

s_df.list <- lapply(X = s_df.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
hepa.anchors <- FindIntegrationAnchors(object.list = s_df.list, dims = 1:20)
s_df <- IntegrateData(anchorset = hepa.anchors, dims = 1:20)

##Perform an integrated analysis
DefaultAssay(s_df) <- "integrated"

# Run the standard workflow for visualization and clustering
s_df <- ScaleData(s_df, verbose = FALSE)
s_df <- RunPCA(s_df, verbose = FALSE)
s_df <- RunUMAP(s_df, reduction = "pca", dims = 1:15, min.dist = 0.3, n.neighbors = 30L)
s_df <- FindNeighbors(s_df, reduction = "pca", dims = 1:15)
s_df <- FindClusters(s_df, resolution = 0.2, algorithm = 1, random.seed = 42)


##Identify cell type markers
s_df = readRDS("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/s_df_filtered.rds")
s_df = SetIdent(s_df, value = s_df[["seurat_clusters"]])
DefaultAssay(s_df) <- "RNA"
clu_markers = FindAllMarkers(s_df, assay = "RNA") %T>% write.delim("DEG_clusters.txt")

UMAPPlot(s_df, label = F, split.by = "orig.ident") 
#UMAPPlot(s_df, label = T, group.by = "seurat_clusters") 
#clu_markers = read.delim("DEG_clusters.txt")

#DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 10 markers (or all markers if less than 10) for each cluster.
clu_markers$cluster = gsub(pattern = "Hepatocyte", replacement = "0", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "CD3\\+ T cell", replacement = "1", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "CD8\\+ T Cell", replacement = "2", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "Germinal Center B Cell", replacement = "3", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "Hematopoietic Stem Cell", replacement = "4", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "Lymphocyte", replacement = "5", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "M2 Macrophage", replacement = "6", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "Macrophage", replacement = "7", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "Regulatory T Cell", replacement = "8", clu_markers$cluster)
clu_markers$cluster = gsub(pattern = "Stromal Cell", replacement = "9", clu_markers$cluster)


library(scCATCH) ## Using version 2.1
clu_anno = scCATCH::scCATCH(clu_markers,
                            species = 'Mouse',
                            cancer = NULL,
                            tissue = c('Liver', 'Fetal liver', 'Blood', 'Peripheral blood', 'Serum', 'Umbilical cord blood'))

s_df$cell_types = s_df@meta.data$seurat_clusters %>% as.character() %>% data.frame() %>%
  left_join(clu_anno %>% dplyr::select(cluster, cell_type), by = c("." = "cluster")) %>% pull(cell_type)
#write.delim(clu_anno, "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/scCATCH_marker_gene_annotation.txt")
# saveRDS(object = s_df, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/s_df_filtered.new.rds")

cell_types_counts = table(s_df$cell_types, s_df$orig.ident) %>%
  format() %>% as.data.frame() %>% rownames_to_column(var = "Cell_Type") %T>%
  write.delim("cell_types_counts.txt")


UMAPPlot(s_df, label = T, group.by = "cell_types")
FeaturePlot(s_df, features = "Cd68", label = F, reduction = "umap", order = T)


##### If read seurat object from RDS #####

##Identify conserved cell type markers
s_df = readRDS("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/s_df_filtered.rds")

s_df = SetIdent(s_df, value = s_df[["seurat_clusters"]])
s_df <- s_df[!grepl("^mt-", rownames(s_df)), ]
DefaultAssay(s_df) <- "RNA"
s_df = ScaleData(s_df)
clu_markers = FindAllMarkers(s_df, assay = "RNA") %T>% write.delim("DEG_clusters.txt")

UMAPPlot(s_df, label = F, split.by = "orig.ident") 
DimPlot(s_df, label = F, split.by = "orig.ident", pt.size = 0.1) 
#UMAPPlot(s_df, label = T, group.by = "seurat_clusters") 
#clu_markers = read.delim("DEG_clusters.txt")

#DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 10 markers (or all markers if less than 10) for each cluster.


library(scCATCH)
clu_anno = scCATCH::scCATCH(clu_markers,
                            species = 'Mouse',
                            cancer = NULL,
                            tissue = c('Liver', 'Fetal liver', 'Blood', 'Peripheral blood', 'Serum', 'Umbilical cord blood'))

s_df$cell_types = s_df@meta.data$seurat_clusters %>% as.character() %>% data.frame() %>%
  left_join(clu_anno %>% dplyr::select(cluster, cell_type), by = c("." = "cluster")) %>% pull(cell_type)
#write.delim(clu_anno, "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/scCATCH_marker_gene_annotation.txt")
# saveRDS(object = s_df, file = "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/s_df_filtered.new.rds")

cell_types_counts = table(s_df$cell_types, s_df$orig.ident) %>%
  format() %>% as.data.frame() %>% rownames_to_column(var = "Cell_Type") %T>%
  write.delim("cell_types_counts.txt")


UMAPPlot(s_df, label = F, group.by = "cell_types", split.by = "orig.ident")


##1## plot heatmap for top 5 DEGs of each clusters
#s_df = readRDS(file = "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.rds")
s_df = SetIdent(s_df, value = s_df[["seurat_clusters"]])
clu_markers_merged = left_join(clu_markers, clu_anno, by = "cluster")
clu_markers_merged %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5 
#s_df = ScaleData(s_df, assay = "RNA")
DoHeatmap(s_df, features = top5$gene, assay = "RNA")
pdf(file = paste0("DoHeatmap_plot_using_scCATCH_annotation.pdf"))
DoHeatmap_plot=DoHeatmap(s_df, features = top5$gene, assay = "RNA")
dev.off()

##2## Louvain clustering of ~30,000 cells from three samples - reassign cell type names based on scCATCH annotation, DEGs and markers
s_df = SetIdent(s_df, value = "seurat_clusters")
s_df.new.ids<-c("Hepatocyte","Hepatocyte","Hepatocyte","Hepatocyte","CD3+ T cell", "Endothelial Cell", "Macrophage", "Germinal Center B Cell", "Macrophage", "Hepatocyte", "CD8+ T Cell", "Lymphocyte", "Stromal Cell", "Regulatory T Cell")
names(s_df.new.ids)<-levels(s_df)
s_df<-RenameIdents(s_df, s_df.new.ids)
s_df$cell_types = s_df@active.ident

pdf(file = paste0("plots/UMAP_raw2.pdf"))
DimPlot(s_df,reduction = "umap",label = F, pt.size = 0.5, split.by = "orig.ident" )
dev.off()
saveRDS(object = s_df, file = "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.rds")

##3## ISG score
s_df = readRDS(file = "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.rds")

genes_isg = c("Irf7","Ifit1","Isg15","Stat1","Ifitm3","Ifi44","Tgtp1","Ifnb1")
s_df = ScaleData(s_df, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
s_df@meta.data$ISG_score = colMeans(s_df[["RNA"]]@scale.data)
p = FeaturePlot(s_df, features = "ISG_score", label = F, reduction = "umap", split.by = "orig.ident", order = T, raster = F, pt.size = 0.7, cols = c("grey90", viridis::inferno(10000, begin = 0, end = 1, direction = -1)))

for(i in 1:nrow(s_df@meta.data)) {
  if(s_df@meta.data[i,8] != "Hepatocyte") {
    s_df@meta.data[i,9] = 0
  }
}
FeaturePlot(s_df, features = "ISG_score", label = F, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"), pt.size = 0.7)


#ggsave(filename = paste0("plots/umap_genes_isg_collapsed_sample.pdf"), plot = plot_genes_isg_collapsed_sample, width = 12, height = 6)

# ISG score boxplot comparing three samples
s_df = ScaleData(s_df, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
s_df@meta.data$ISG_score = colMeans(s_df[["RNA"]]@scale.data)
df = s_df@meta.data %>% select(orig.ident, ISG_score) %>% mutate(log_ISG= log10(ISG_score + 1)) 

my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
ggboxplot(df, x = "orig.ident", y = "ISG_score", ylab = "Normalized ISG expression level (log10(x+1))",
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"))+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

boxplot(log_ISG ~ orig.ident, data = df, #ylim=c(-0.15,0.55), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#F8766D", "#00BA38", "#619CFF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'ISG score of the three samples')


## cytotoxic signature score
genes_cytotox = c("Prf1", "Gzmb")
s_df = ScaleData(s_df, assay = "RNA", features = genes_cytotox, do.center = T, do.scale = T)
s_df@meta.data$cytotox_score = colMeans(s_df[["RNA"]]@scale.data)
FeaturePlot(s_df, features = "cytotox_score", label = T, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"), pt.size = 1)

# cytotoxic score boxplot comparing three samples
s_df = ScaleData(s_df, assay = "RNA", features = genes_cytotox, do.center = F, do.scale = T)
s_df@meta.data$cytotox_score = colMeans(s_df[["RNA"]]@scale.data)
df = s_df@meta.data %>% select(orig.ident, cytotox_score) %>% mutate(log_cytotox= log10(cytotox_score + 1)) 

my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
ggboxplot(df, x = "orig.ident", y = "cytotox_score", ylab = "Normalized cytotoxic score", 
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"))+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

boxplot(log_cytotox ~ orig.ident, data = df, #ylim=c(-0.15,0.55), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#F8766D", "#00BA38", "#619CFF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Cytotoxic score of the three samples')

## exhaustion signature score
genes_exhaustion = c("Pdcd1", "Ctla4", "Lag3", "Cd160", "Havcr2", "Cd244")
s_df = ScaleData(s_df, assay = "RNA", features = genes_exhaustion, do.center = T, do.scale = T)
s_df@meta.data$exhaustion_score = colMeans(s_df[["RNA"]]@scale.data)
FeaturePlot(s_df, features = "exhaustion_score", label = T, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"), pt.size = 0.6)

# exhaustion score boxplot comparing three samples
s_df = ScaleData(s_df, assay = "RNA", features = genes_exhaustion, do.center = F, do.scale = T)
s_df@meta.data$exhaustion_score = colMeans(s_df[["RNA"]]@scale.data)
df = s_df@meta.data %>% select(orig.ident, exhaustion_score) %>% mutate(log_exhaustion= log10(exhaustion_score + 1)) 

my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
ggboxplot(df, x = "orig.ident", y = "exhaustion_score", ylab = "Normalized exhaustion score", 
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"))+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

boxplot(log_exhaustion ~ orig.ident, data = df, #ylim=c(-0.15,0.55), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#F8766D", "#00BA38", "#619CFF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Exhaustion score of the three samples')


##4## Cd68 gene plotting colored by samples in a merged plot
s_df_new = SetIdent(s_df, value = s_df[["orig.ident"]])
Cd68_expressed_A = WhichCells(object = s_df_new, expression = Cd68 > 0, idents = "1")
Cd68_expressed_B = WhichCells(object = s_df_new, expression = Cd68 > 0, idents = "2")
Cd68_expressed_C = WhichCells(object = s_df_new, expression = Cd68 > 0, idents = "3")

Cd68_unexpressed_A = WhichCells(object = s_df_new, expression = Cd68 == 0, idents = "1")
Cd68_unexpressed_B = WhichCells(object = s_df_new, expression = Cd68 == 0, idents = "2")
Cd68_unexpressed_C = WhichCells(object = s_df_new, expression = Cd68 == 0, idents = "3")

DimPlot(s_df_new, reduction = "umap", pt.size = 0.7, label = F, group.by = "orig.ident", label.size = 5, cells.highlight= list(Cd68_expressed_C, Cd68_expressed_B, Cd68_expressed_A), cols.highlight = c("#F8766D", "#00BA38", "#619CFF"), cols= "grey")


#pdf(file = paste0("plots/Cd68_plot_merged.pdf"), width = 8, height = 7)
#DP
#dev.off()
#DimPlot(s_df_new, reduction = "umap", pt.size = 0.7, label = F, group.by = "orig.ident", label.size = 5, cells.highlight= list(Cd68_expressed_C, Cd68_unexpressed_C), cols.highlight = c("grey", "#619CFF"), cols= "grey")


##5## Find the DEGs between samples within one cluster
# Save current cluster identites in object@meta.data under 'cell_types'
# Run only if Seurat::FindClusters() was executed
s_df_new = SetIdent(s_df, value = s_df[["cell_types"]])
object <- Seurat::StashIdent(object = s_df_new, save.name = "cell_types")

# Attribute sample identity to each cell type
sampleID <- paste(object@meta.data$cell_types, object@meta.data$orig.ident, sep = "_")

# Add sample identity to object@meta.data under 'sampleID'
object@meta.data$sampleID <- sampleID

# Set cell identity to sample identity
Idents(object = object) <- object@meta.data$sampleID

# Find all sample specific marker genes 
#markers <- Seurat::FindMarkers(object = object,ident.1 = "Hepatocyte_C_Hom", ident.2 = c("Hepatocyte_A_WT", "Hepatocyte_B_Het"))

## Do the above command in loop:
#for each cluster name: load file and create outline
allmarkers_from_C_compared_to_A_B <- lapply(levels(s_df_new), function(cluster_name) {
  markers <- Seurat::FindMarkers(object = object,ident.1 = paste0(cluster_name, "_3"), ident.2 = c(paste0(cluster_name, "_1"), paste0(cluster_name, "_2")), assay = "RNA") %>% mutate(cluster = cluster_name) %>% tibble::rownames_to_column(var = "geneName") 
  markers
})

# merge all outlines into one data frame (by appending them row-wise)
allmarkers_from_C_compared_to_A_B.merged <- do.call(rbind, allmarkers_from_C_compared_to_A_B) %T>% write.delim("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/DEG_clusters_from_C_compared_to_A_B.scCATCH_annotated.cluster11_change_to_lymphcytes.txt")

## Do the above command in loop:
#for each cluster name: load file and create outline
allmarkers_from_C_compared_to_A <- lapply(levels(s_df_new), function(cluster_name) {
  markers <- Seurat::FindMarkers(object = object,ident.1 = paste0(cluster_name, "_3"), ident.2 = paste0(cluster_name, "_1"), assay = "RNA") %>% mutate(cluster = cluster_name) %>% tibble::rownames_to_column(var = "geneName")
  markers
})

# merge all outlines into one data frame (by appending them row-wise)
allmarkers_from_C_compared_to_A.merged <- do.call(rbind, allmarkers_from_C_compared_to_A) %T>% write.delim("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/DEG_clusters_from_C_compared_to_A.scCATCH_annotated.cluster11_change_to_lymphcytes.txt")


##6## check the two sub-clusters of Macrophage using ToppGene
library(data.table)
s_df_new = s_df
Idents(s_df_new) = s_df_new$seurat_clusters
degs = FindMarkers(s_df_new, ident.1 = 11)
degs %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames_to_column('gene') %>% 
  fwrite('cl11_degs.macrophage.csv')

degs = FindMarkers(s_df_new, ident.1 = 8)
degs %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames_to_column('gene') %>% 
  fwrite('cl8_degs.macrophage.csv')

degs = FindMarkers(s_df_new, ident.1 = 6)
degs %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames_to_column('gene') %>% 
  fwrite('cl6_degs.macrophage.csv')


##7## featureplot for certain groups of cell types
s_df_new = SetIdent(s_df, value = s_df[["orig.ident"]])
highligted_cells1 = WhichCells(object = s_df_new, expression = cell_types == "Macrophage", idents = c("1", "2", "3"))
highligted_cells2 = WhichCells(object = s_df_new, expression = cell_types == "Hepatocyte", idents = c("1", "2", "3"))

#DimPlot(s_df, reduction = "umap", label = F, split.by = "orig.ident", label.size = 5, cells.highlight= list(highligted_cells1, highligted_cells2), cols= "grey")
FeaturePlot(object = s_df, features = c('Grn', 'Egfr'), cells = c(highligted_cells1, highligted_cells2), blend = TRUE)

p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells1, split.by = "orig.ident", pt.size = 1.3, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
fix.sc <- scale_color_gradientn( colours = c('gray100', 'red'),  limits = c(1, 3))
p2 <- lapply(p1, function (x) x + fix.sc + xlim(c(8, 11.5)) + ylim(c(-7, -1)))
CombinePlots(p2, ncol = 3)

DefaultAssay(s_df) <- "RNA"
p1 = FeaturePlot(object = s_df, features = c('Grn'), cells = highligted_cells2, split.by = "orig.ident", pt.size = 0.5, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
fix.sc <- scale_color_gradientn( colours = c('gray100', 'red'),  limits = c(1, 3))
p2 <- lapply(p1, function (x) x + fix.sc + xlim(c(-10, 8)) + ylim(c(-10, 15)))
CombinePlots(p2, ncol = 3)

p1 = FeaturePlot(object = s_df, features = c('Tnfrsf1b'), cells = highligted_cells1, split.by = "orig.ident", pt.size = 1.3, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
fix.sc <- scale_color_gradientn( colours = c('gray100', 'blue'),  limits = c(0.2, 4.5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2, ncol = 3)

p1 = FeaturePlot(object = s_df, features = c('Cd209a'), cells = highligted_cells1, split.by = "orig.ident", pt.size = 1.3, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
fix.sc <- scale_color_gradientn( colours = c('gray100', 'blue'),  limits = c(1, 3))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2, ncol = 3)

highligted_cells3 = WhichCells(object = s_df_new, expression = cell_types == "CD3+ T cell", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells3, split.by = "orig.ident", pt.size = 0.8, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

highligted_cells4 = WhichCells(object = s_df_new, expression = cell_types == "CD8+ T Cell", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells4, split.by = "orig.ident", pt.size = 1.3, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

highligted_cells5 = WhichCells(object = s_df_new, expression = cell_types == "Regulatory T Cell", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells5, split.by = "orig.ident", pt.size = 1.5, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
fix.sc <- scale_color_gradientn( colours = c('gray100', 'blue'),  limits = c(0.1, 2.5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2, ncol = 3)

p1 = FeaturePlot(object = s_df, features = c('Itga4'), cells = highligted_cells5, split.by = "orig.ident", pt.size = 1.5, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

highligted_cells6 = WhichCells(object = s_df_new, expression = cell_types == "Germinal Center B Cell", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells6, split.by = "orig.ident", pt.size = 1.3, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

highligted_cells7 = WhichCells(object = s_df_new, expression = cell_types == "Lymphocyte", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells7, split.by = "orig.ident", pt.size = 1.8, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

highligted_cells8 = WhichCells(object = s_df_new, expression = cell_types == "Stromal Cell", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells8, split.by = "orig.ident", pt.size = 1.2, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)
p1 = FeaturePlot(object = s_df, features = c('Tnfrsf1a'), cells = highligted_cells8, split.by = "orig.ident", pt.size = 1.2, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

highligted_cells9 = WhichCells(object = s_df_new, expression = cell_types == "Hematopoietic Stem Cell", idents = c("1", "2", "3"))
p1 = FeaturePlot(object = s_df, features = c('Egfr'), cells = highligted_cells9, split.by = "orig.ident", pt.size = 1.2, combine = FALSE, order = TRUE)# & theme(legend.position = c(1,0.2))
CombinePlots(p1, ncol = 3)

##8## Basic box plot
library(ggpubr)
s_df = SetIdent(s_df, value = "cell_types")
s_df@meta.data$Grn_exp = s_df@assays$RNA@counts["Grn", ]
s_df@meta.data$Egfr_exp = s_df@assays$RNA@counts["Egfr", ]
#s_df@meta.data$Grn_exp = s_df@assays$integrated@data["Grn", ]
#s_df@meta.data$Egfr_exp = s_df@assays$integrated@data["Egfr", ]
s_df_Mac = subset(s_df, idents = c("Macrophage"))
s_df_Hepa = subset(s_df, idents = c("Hepatocyte"))

df_Mac = s_df_Mac@meta.data %>% dplyr::select(orig.ident, Egfr_exp) %>% mutate(log_Egfr_exp = log10(Egfr_exp + 1)) #%>% mutate(Egfr_exp_1 = Egfr_exp +1)
df_Hepa = s_df_Hepa@meta.data %>% dplyr::select(orig.ident, Grn_exp) %>% mutate(log_Grn_exp = log10(Grn_exp + 1)) #%>% mutate(Grn_exp_1 = Grn_exp +1)
write.table(df_Hepa, file = "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/table for plotting Grn expression level in hepatocytes.txt", row.names = F, col.names = T, sep = '\t')
write.table(df_Mac, file = "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/table for plotting Egfr expression level in Macrophage.txt", row.names = F, col.names = T, sep = '\t')

my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
p = ggboxplot(df_Mac, x = "orig.ident", y = "Egfr_exp",
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"), add = "jitter", shape = "orig.ident")+
  stat_compare_means(comparisons = my_comparisons)+
  scale_y_log10(breaks=c(0,1,2,3,4,5))+
  ylab("Egfr expression level in logarithmic y-axis")
p

p = ggboxplot(df_Hepa, "orig.ident", "Grn_exp_1", 
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"), add = "jitter", shape = "orig.ident")+
  stat_compare_means(comparisons = my_comparisons)+
  scale_y_log10(breaks=c(5,4,3,2,1,0))+
  ylab("Grn expression level in logarithmic y-axis")
p

breaks = seq(0, 2, 2)
labels <- seq(0, 2 ,2)
p <- p + facet_grid(scales="free", space="free")
p <- p + scale_y_continuous(breaks=breaks, labels=labels, expand=c(0.075,0))
p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())
p

boxplot(Grn_exp ~ orig.ident, data = df_Hepa,
        main="Different boxplots for each month",
        xlab="Month Number",
        ylab="Degree Fahrenheit",
        border=c("#F8766D", "#00BA38", "#619CFF")
)
library(beeswarm)
beeswarm(Grn_exp ~ orig.ident, data = df_Hepa, #ylim=c(0.69,18),  #at=c(1.5,2.5), xlim=c(0.8,3.2),
         col=c("grey3", "blue2", "purple3"), pch = 16, add = TRUE, spacing = 0.75)

##9## violin plot
VlnPlot(s_df, features = c("Grn", "Egfr"), split.by = "orig.ident", group.by = "cell_types", pt.size = 0.2, combine = FALSE)

##10## Plot a series of feature plots without x or y axis
p = FeaturePlot(s_df, features = c("Alb","Cd68", "S100a4","Itgam","Cd8b1","Cd3d","Cd3e","Cd3g"), label = F, reduction = "umap", order = T, raster = F, combine = FALSE, min.cutoff = 5)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
p2 = cowplot::plot_grid(plotlist = p)

ggsave("C:/Users/renxi/Desktop/WL analysis/eight_feature_plots.pdf", p2, width = 12, height = 2)

# for albumin gene and other genes
FeaturePlot(s_df,"Alb",order = T, min.cutoff = 5)
FeaturePlot(s_df, features = c("Cd68", "S100a4", "Itgam", "Cd8b1", "Cd3d", "Cd3e", "Cd3g"), label = F, reduction = "umap", order = T, raster = F)

##11## featureplot for co-expression of Grn and Egfr using mouse scRNAseq
FeaturePlot(s_df, features = c("Grn", "Egfr"), blend = TRUE, order = T, split.by = "orig.ident")
