library(caroline)
library(magrittr)
library(tidyverse)
library(openxlsx)
library(Seurat)
library(viridis)
library(patchwork)
library(scCATCH)
library(rafalib)
library(clustree)
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)

#s_df = readRDS(url("https://zenodo.org/record/3531889/files/s_df.rds"))
#s_df@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513


#Read in NicheNet's ligand-target prior model, ligand-receptor network and weighted integrated networks:
ligand_target_matrix = readRDS("C:/Users/renxi/Desktop/e0149673/r_scripts/PLC lab/Nichenet_network_files/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network = readRDS("C:/Users/renxi/Desktop/e0149673/r_scripts/PLC lab/Nichenet_network_files/lr_network.rds")
head(lr_network)
## ???[38;5;246m# A tibble: 6 x 4???[39m
##   from  to    source         database
##   ???[3m???[38;5;246m<chr>???[39m???[23m ???[3m???[38;5;246m<chr>???[39m???[23m ???[3m???[38;5;246m<chr>???[39m???[23m          ???[3m???[38;5;246m<chr>???[39m???[23m   
## ???[38;5;250m1???[39m CXCL1 CXCR2 kegg_cytokines kegg    
## ???[38;5;250m2???[39m CXCL2 CXCR2 kegg_cytokines kegg    
## ???[38;5;250m3???[39m CXCL3 CXCR2 kegg_cytokines kegg    
## ???[38;5;250m4???[39m CXCL5 CXCR2 kegg_cytokines kegg    
## ???[38;5;250m5???[39m PPBP  CXCR2 kegg_cytokines kegg    
## ???[38;5;250m6???[39m CXCL6 CXCR2 kegg_cytokines kegg

weighted_networks = readRDS("C:/Users/renxi/Desktop/e0149673/r_scripts/PLC lab/Nichenet_network_files/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## ???[38;5;246m# A tibble: 6 x 3???[39m
##   from  to     weight
##   ???[3m???[38;5;246m<chr>???[39m???[23m ???[3m???[38;5;246m<chr>???[39m???[23m   ???[3m???[38;5;246m<dbl>???[39m???[23m
## ???[38;5;250m1???[39m A1BG  ABCC6  0.422 
## ???[38;5;250m2???[39m A1BG  ACE2   0.101 
## ???[38;5;250m3???[39m A1BG  ADAM10 0.097???[4m0???[24m
## ???[38;5;250m4???[39m A1BG  AGO1   0.052???[4m5???[24m
## ???[38;5;250m5???[39m A1BG  AKT1   0.085???[4m5???[24m
## ???[38;5;250m6???[39m A1BG  ANXA7  0.457
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## ???[38;5;246m# A tibble: 6 x 3???[39m
##   from  to     weight
##   ???[3m???[38;5;246m<chr>???[39m???[23m ???[3m???[38;5;246m<chr>???[39m???[23m   ???[3m???[38;5;246m<dbl>???[39m???[23m
## ???[38;5;250m1???[39m A1BG  A2M    0.029???[4m4???[24m
## ???[38;5;250m2???[39m AAAS  GFAP   0.029???[4m0???[24m
## ???[38;5;250m3???[39m AADAC CYP3A4 0.042???[4m2???[24m
## ???[38;5;250m4???[39m AADAC IRF8   0.027???[4m5???[24m
## ???[38;5;250m5???[39m AATF  ATM    0.033???[4m0???[24m
## ???[38;5;250m6???[39m AATF  ATR    0.035???[4m5???[24m

#convert the NicheNet network gene symbols from human to mouse based on one-to-one orthology:
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()





s_df_new = SetIdent(s_df, value = s_df[["cell_types"]])

## Perform the NicheNet analysis
#1. Define a "sender/niche" cell population and a "receiver/target" cell population present in your expression data and determine which genes are expressed in both populations
## receiver
receiver = c("Macrophage")
#receiver = c("Macrophage")
list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes, s_df_new, 0.10, assay_oi = "integrated")
expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Hepatocyte")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, s_df_new, 0.10, assay_oi = "integrated") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


#2. Define a gene set of interest: these are the genes in the "receiver/target" cell population that are potentially affected by ligands expressed by interacting cells (e.g. genes differentially expressed upon cell-cell interaction)
seurat_obj_receiver= subset(s_df_new, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "3" #sender
condition_reference = c("1", "2") #reference

#DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0, logfc.threshold = 0.1, assay = "RNA") %>% rownames_to_column("gene")
DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10, assay = "RNA") %>% rownames_to_column("gene")
potential_receptor = DE_table_receiver %>% filter((p_val_adj <= 0.05 | p_val < 0.005) & abs(avg_log2FC) >= 0.25) %>% pull(gene)
#potential_receptor = DE_table_receiver %>% pull(gene)
geneset_oi = potential_receptor %>% .[. %in% rownames(ligand_target_matrix)]

#3. Define a set of potential ligands: these are ligands that are expressed by the "sender/niche" cell population and bind a (putative) receptor expressed by the "receiver/target" population
seurat_obj_ligands= subset(s_df_new, idents = sender_celltypes)
seurat_obj_ligands = SetIdent(seurat_obj_ligands, value = seurat_obj_ligands[["orig.ident"]])

condition_oi = "3" #sender
condition_reference = c("1", "2") #reference

#DE_table_ligands = FindMarkers(object = seurat_obj_ligands, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0, logfc.threshold = 0.1, assay = "RNA") %>% rownames_to_column("gene")
DE_table_ligands = FindMarkers(object = seurat_obj_ligands, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10, assay = "RNA") %>% rownames_to_column("gene")
#DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

potential_ligand = DE_table_ligands %>% filter((p_val_adj <= 0.05 | p_val < 0.005) & abs(avg_log2FC) >= 0.25) %>% pull(gene)
#potential_ligand = DE_table_ligands %>% pull(gene)
potential_ligands = lr_network %>% filter(from %in% potential_ligand & to %in% potential_receptor) %>% pull(from) %>% unique()

#4) Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
## ???[38;5;246m# A tibble: 44 x 5???[39m
##    test_ligand auroc  aupr pearson  rank
##    ???[3m???[38;5;246m<chr>???[39m???[23m       ???[3m???[38;5;246m<dbl>???[39m???[23m ???[3m???[38;5;246m<dbl>???[39m???[23m   ???[3m???[38;5;246m<dbl>???[39m???[23m ???[3m???[38;5;246m<dbl>???[39m???[23m
## ???[38;5;250m 1???[39m Ebi3        0.638 0.234  0.197      1
## ???[38;5;250m 2???[39m Il15        0.582 0.163  0.096???[4m1???[24m     2
## ???[38;5;250m 3???[39m Crlf2       0.549 0.163  0.075???[4m8???[24m     3
## ???[38;5;250m 4???[39m App         0.499 0.141  0.065???[4m5???[24m     4
## ???[38;5;250m 5???[39m Tgfb1       0.494 0.140  0.055???[4m8???[24m     5
## ???[38;5;250m 6???[39m Ptprc       0.536 0.149  0.055???[4m4???[24m     6
## ???[38;5;250m 7???[39m H2-M3       0.525 0.157  0.052???[4m8???[24m     7
## ???[38;5;250m 8???[39m Icam1       0.543 0.142  0.048???[4m6???[24m     8
## ???[38;5;250m 9???[39m Cxcl10      0.531 0.141  0.040???[4m8???[24m     9
## ???[38;5;250m10???[39m Adam17      0.517 0.137  0.035???[4m9???[24m    10
## ???[38;5;246m# ... with 34 more rows???[39m

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
#best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

#DotPlot(s_df_new, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()


#5) Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
#pdf(file = paste0("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/NicheNet/p_ligand_target_network.in-house_mouse_scRNAseq.sender_Macrophage.receiver_Hepatocyte.pdf"))

write.delim(vis_ligand_target %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_name"), "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/NicheNet/NicheNet_regulatory_potential.Hepatocyte_Macrophage.C_vs_A_B.txt")

pdf(file = paste0("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/NicheNet/p_ligand_target_network.in-house_mouse_scRNAseq.sample1_compared_to_sample3.sender_Hepatocyte.receiver_Macrophage.pdf"))
p_ligand_target_network
dev.off()
p_ligand_target_network

##Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% potential_receptor) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()


p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network



#write.delim(vis_ligand_receptor_network %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene_name"), "C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/NicheNet/NicheNet_prior_interaction_potential.Hepatocyte_Macrophage.C_vs_A_B.txt")

vis_ligand_receptor_network_final = read.delim(paste0("C:/Users/renxi/Desktop/WL analysis/VisiumBatch1(19Nov2021)/analysis_output/NicheNet/NicheNet_prior_interaction_potential.Hepatocyte_Macrophage.C_vs_A_B.txt"), stringsAsFactors = F, header = T, check.names = F) %>% tibble::column_to_rownames(var = "gene_name") %>% as.matrix()

p_ligand_receptor_network = vis_ligand_receptor_network_final %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
