# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

s_data_list = Read10X(data.dir = paste0("C:/Users/renxi/Desktop/WL analysis/in_house_HCC_and_normal_scRNASeq/aggregate/filtered_feature_bc_matrix_by_NGS_portal/"))
s_df <- CreateSeuratObject(counts = s_data_list, project = "HCC_normal")
s_df[["percent.mito"]] = PercentageFeatureSet(s_df, pattern = "^MT-")
#s_df[["percent.ribo"]] = PercentageFeatureSet(s_df, pattern = "^RP[SL]")
min_nFeature_RNA = 100
percent_mito = 20
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



# scale and run PCA
s_df <- ScaleData(s_df, features = rownames(s_df))
s_df <- RunPCA(s_df, features = VariableFeatures(object = s_df))

# Check number of PC components 
ElbowPlot(s_df)

# cluster and visualize
s_df <- FindNeighbors(s_df, dims = 1:15)
s_df <- FindClusters(s_df, resolution = 0.2, algorithm = 1, random.seed = 42)
s_df <- RunUMAP(s_df, reduction = "pca", dims = 1:15, min.dist = 0.3, n.neighbors = 30L)
DimPlot(s_df, reduction = "umap", group.by = "seurat_clusters", label = T)


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = c("Liver", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = s_df[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(s_df@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(s_df@meta.data[s_df@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(s_df@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay the identified cell types on UMAP plot:
s_df@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  s_df@meta.data$customclassif[s_df@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
#Plotting
DimPlot(s_df, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')  
DimPlot(s_df, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', split.by = "orig.ident")    

##1## Based on DEGs and marker genes, we re-annotate some cell types from sc-Type
s_df@meta.data$customclassif = gsub(pattern = "Kupffer cells", replacement = "Macrophage", s_df@meta.data$customclassif)
s_df@meta.data$customclassif = gsub(pattern = "Memory CD8\\+ T cells", replacement = "CD8+ T cells", s_df@meta.data$customclassif)
s_df@meta.data$customclassif = gsub(pattern = "Memory B cells", replacement = "B cells", s_df@meta.data$customclassif)
s_df@meta.data$customclassif = gsub(pattern = "Myeloid Dendritic cells", replacement = "Dendritic cells", s_df@meta.data$customclassif)

saveRDS(s_df, "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.sc_type_annotated.from_NGS_portal.rds")
#saveRDS(s_df, "C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.sc_type_annotated.using individual cell ranger result instead of using the aggregated result from cell ranger.rds")
s_df = readRDS("C:/Users/renxi/Desktop/WL analysis/Wrapping_up_manuscript_WL_requests/s_df.sc_type_annotated.from_NGS_portal.rds")

# calculate fold change in DEG lists of hepatocytes and macrophage
s_df = SetIdent(s_df, value = s_df[["customclassif"]])
s_df_hep = subset(s_df, idents = "Hepatocytes")
s_df_hep = SetIdent(s_df_hep, value = s_df_hep[["orig.ident"]])
DE_table_receiver = FindMarkers(object = s_df_hep, ident.1 = "2", ident.2 = "1", min.pct = 0.10, assay = "RNA") %>% rownames_to_column("gene")
DE_table_receiver[DE_table_receiver$gene %in% "GRN",]

s_df_mac = subset(s_df, idents = "Kupffer cells")
s_df_mac  = SetIdent(s_df_mac, value = s_df_mac[["orig.ident"]])
DE_table_receiver = FindMarkers(object = s_df_mac, ident.1 = "2", ident.2 = "1", min.pct = -Inf, logfc.threshold = -Inf, assay = "RNA") %>% rownames_to_column("gene")
DE_table_receiver[DE_table_receiver$gene %in% "EGFR",]




##2## ISG score
#genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1","IFI27","IIGP1","LPIN1","LY6A","LY6E")
genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")

s_df = ScaleData(s_df, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
s_df@meta.data$ISG_score = colMeans(s_df[["RNA"]]@scale.data)
FeaturePlot(s_df, features = "ISG_score", label = F, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"))

s_df_new = SetIdent(s_df, value = s_df[["customclassif"]])
hepotacyte_ids = WhichCells(object = s_df_new, idents = "Hepatocytes")

cells.ident <- FetchData(s_df, vars="ident")
hepotacyte_ids = grep("Hepatocytes", cells.ident$ident, ignore.case=T)
FeaturePlot(s_df_new, features = "ISG_score", label = F, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"), cells = hepotacyte_ids)


# ISG score boxplot comparing two samples
s_df_hep = subset(s_df, idents = "Hepatocytes")
s_df_hep = ScaleData(s_df_hep, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
s_df_hep@meta.data$ISG_score = colMeans(s_df_hep[["RNA"]]@scale.data)
df = s_df_hep@meta.data %>% select(orig.ident, ISG_score) %>% mutate(log_ISG= log10(ISG_score + 1)) 

my_comparisons <- list( c("1", "2") )
library(ggpubr)
ggboxplot(df, x = "orig.ident", y = "ISG_score", ylab = "Normalized ISG score",
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"))+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

boxplot(log_ISG ~ orig.ident, data = df, #ylim=c(-0.15,0.55), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#F8766D", "#00BA38", "#619CFF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'ISG score of the two samples')

# ISG score violin plot comparing two samples
s_df_hep = subset(s_df, idents = "Hepatocytes")
VlnPlot(s_df_hep, features = "ISG_score", split.by = "orig.ident", pt.size = 0)

s_df = SetIdent(s_df, value = "orig.ident")
VlnPlot(s_df, features = "ISG_score", split.by = "orig.ident", pt.size = 0.2)
s_df = SetIdent(s_df, value = s_df[["customclassif"]])
VlnPlot(s_df, features = "ISG_score", split.by = "orig.ident", pt.size = 0.2)

##2-2 cytotoxic signature score
genes_cytotox = c("PRF1", "GZMH", "GZMB")
s_df = ScaleData(s_df, assay = "RNA", features = genes_cytotox, do.center = T, do.scale = T)
s_df@meta.data$cytotox_score = colMeans(s_df[["RNA"]]@scale.data)
FeaturePlot(s_df, features = "cytotox_score", label = T, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"))

# cytotoxic score boxplot comparing three samples
s_df = ScaleData(s_df, assay = "RNA", features = genes_cytotox, do.center = F, do.scale = T)
s_df@meta.data$cytotox_score = colMeans(s_df[["RNA"]]@scale.data)
df = s_df@meta.data %>% select(orig.ident, cytotox_score) %>% mutate(log_cytotox= log10(cytotox_score + 1)) 

my_comparisons <- list( c("1", "2") )
ggboxplot(df, x = "orig.ident", y = "cytotox_score", ylab = "Normalized cytotoxic score", 
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"))+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

boxplot(log_cytotox ~ orig.ident, data = df, #ylim=c(-0.15,0.55), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#F8766D", "#00BA38", "#619CFF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Cytotoxic score of the two samples')

##2-3 exhaustion signature score
genes_exhaustion = c("PDCD1", "CTLA4", "LAG3", "CD160", "TIGHT", "HAVCR2", "CD244")
s_df = ScaleData(s_df, assay = "RNA", features = genes_exhaustion, do.center = T, do.scale = T)
s_df@meta.data$exhaustion_score = colMeans(s_df[["RNA"]]@scale.data)
FeaturePlot(s_df, features = "exhaustion_score", label = T, reduction = "umap", split.by = "orig.ident", order = T, raster = F, cols = c("honeydew3", "red4"))

# exhaustion score boxplot comparing three samples
s_df = ScaleData(s_df, assay = "RNA", features = genes_exhaustion, do.center = F, do.scale = T)
s_df@meta.data$exhaustion_score = colMeans(s_df[["RNA"]]@scale.data)
df = s_df@meta.data %>% select(orig.ident, exhaustion_score) %>% mutate(log_exhaustion= log10(exhaustion_score + 1)) 

my_comparisons <- list( c("1", "2") )
ggboxplot(df, x = "orig.ident", y = "exhaustion_score", ylab = "Normalized exhaustion score", 
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"))+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

boxplot(log_exhaustion ~ orig.ident, data = df, #ylim=c(-0.15,0.55), #at=c(1.5,2.5), xlim=c(0.8,3.2),
        col=c("#F8766D", "#00BA38", "#619CFF"),
        ylab = "Normalized ISG expression level (log10(x+1))",
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Exhaustion score of the two samples')

##3## Check on GRN expression split by HCC ISG high and HCC ISG low
s_df = SetIdent(s_df, value = s_df[["customclassif"]])
genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
s_df = ScaleData(s_df, assay = "RNA", features = genes_isg, do.center = F, do.scale = T)
s_df@meta.data$ISG_score = colMeans(s_df[["RNA"]]@scale.data)
s_df_hep = subset(s_df, idents = "Hepatocytes")
VlnPlot(s_df_hep, features = "ISG_score", split.by = "orig.ident", pt.size = 0.2, combine = FALSE, assay = "integrated")

ISG = s_df_hep@meta.data %>% subset(orig.ident %in% "2") %>% subset(customclassif %in% "Hepatocytes") %>% select(ISG_score) %>% t() %>% as.data.frame()
quantiles = ISG %>% as.matrix() %>% quantile(seq(0, 1, 1/4)) %>% as.numeric()

ids=s_df_hep$orig.ident
group_info= ifelse(ids == "1", "Normal", ifelse(s_df_hep$ISG_score >= quantiles[4], "HCC - high ISG score", ifelse(s_df_hep$ISG_score <= quantiles[2], "HCC - low ISG score", "HCC - medium ISG score")))
s_df_hep@meta.data$group = group_info
s_df_hep = SetIdent(s_df_hep, value = s_df_hep[["group"]])
VlnPlot(s_df_hep, features = "ISG_score", pt.size = 0.2, assay = "integrated")

s_df_hep@meta.data$GRN = s_df_hep@assays$integrated@data[rownames(s_df_hep@assays$integrated@data) %in% "GRN",]
VlnPlot(s_df_hep, features = "GRN", pt.size = 0, assay = "integrated")

##4## visualize a bubble plot showing all the cell types that were considered by ScType for cluster annotation ##
# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

scater::multiplot(DimPlot(s_df, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)


##5## Basic box plot for GRN and EGFR
library(ggpubr)
s_df = SetIdent(s_df, value = "customclassif")
s_df@meta.data$GRN_exp = s_df@assays$integrated@data["GRN", ]
s_df@meta.data$EGFR_exp = s_df@assays$integrated@data["EGFR", ]
s_df_Mac = subset(s_df, idents = c("Macrophage"))
s_df_Hepa = subset(s_df, idents = c("Hepatocytes"))

df_Mac = s_df_Mac@meta.data %>% dplyr::select(orig.ident, EGFR_exp)# %>% mutate(log_Egfr_exp = log10(EGFR_exp + 1)) 
df_Hepa = s_df_Hepa@meta.data %>% dplyr::select(orig.ident, GRN_exp)# %>% mutate(log_Grn_exp = log10(GRN_exp + 1)) 

my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
ggboxplot(df_Mac, x = "orig.ident", y = "EGFR_exp",
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"), add = "jitter", shape = "orig.ident")+
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(df_Hepa, "orig.ident", "GRN_exp",
          color = "orig.ident", palette = c("#F8766D", "#00BA38", "#619CFF"), add = "jitter", shape = "orig.ident")+
  stat_compare_means(comparisons = my_comparisons)

## Use our custom way to calculate log2FC and p value
s_df_hep = s_df_hep = subset(s_df, idents = "Hepatocytes")
df = s_df_hep@meta.data
x=df[df$orig.ident %in% "1",11]
y=df[df$orig.ident %in% "2",11]
t.test(y,x)
wilcox.test(y,x)
df %>% group_by(orig.ident) %>% summarise(mean=mean(GRN_exp))
log2(0.445)-log2(0.622)
s_df_hep = SetIdent(s_df_hep, value = s_df_hep[["orig.ident"]])
VlnPlot(s_df_hep, features = "GRN", pt.size = 0, assay = "integrated")

s_df_Mac = s_df_Mac = subset(s_df, idents = "Macrophage")
df = s_df_Mac@meta.data
x=df[df$orig.ident %in% "1",12]
y=df[df$orig.ident %in% "2",12]
t.test(y,x)
wilcox.test(y,x)
df %>% group_by(orig.ident) %>% summarise(mean=mean(EGFR_exp))
log2(0.0450)-log2(0.0764)
s_df_Mac = SetIdent(s_df_Mac, value = s_df_Mac[["orig.ident"]])
VlnPlot(s_df_Mac, features = "EGFR", pt.size = 0, assay = "integrated", y.max = 1)


##6## violin plot
VlnPlot(s_df, features = c("GRN", "EGFR"), split.by = "orig.ident", group.by = "customclassif", pt.size = 0.2, combine = FALSE)

