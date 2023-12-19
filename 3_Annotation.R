#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)

#setwd("/mnt/data/guseri/10x_analysi/All_sampless")
Project_name = "x_All"
FindMarker.flag = TRUE
Input.dir = "Output/1_Integrated/"
Norm.assay = "SCT"

if (dir.exists(path = paste0("Output/3_Labeling")) == FALSE) {
  print(paste0("Generating output directory Output/3_Labeling"))
  dir.create(path = paste0("Output/3_Labeling"), recursive = TRUE)
  Output.dir = paste0("Output/3_Labeling/")
} else if (dir.exists(path = paste0("Output/3_Labeling")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/3_Labeling/")
} else {
  print("Error with output directory")
}

# Plotting QC metrics 
Plotting_QC <- function(x, group = "") {
  
  # Generating ridgeplots
  RidgePlot(x.integrated, features = "percent.mt", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_mtDNA.pdf"), dpi = 700)
  
  RidgePlot(x.integrated, features = "percent.ribo", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_ribo.pdf"), dpi = 700)
  
  RidgePlot(x.integrated, features = "percent.hb", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_hb.pdf"), dpi = 700)
  
  RidgePlot(x.integrated, features = "nCount_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_nCount_RNA.pdf"), dpi = 700)
  
  RidgePlot(x.integrated, features = "nFeature_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_nFeature_RNA.pdf"), dpi = 700)
  
  RidgePlot(x.integrated, features = "S.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_S-Score.pdf"), dpi = 700)
  
  RidgePlot(x.integrated, features = "G2M.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_G2M-Score.pdf"), dpi = 700)
  
  # Generating violin plots
  VlnPlot(x.integrated, features = "percent.mt", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_mtDNA.pdf"), dpi = 700)
  
  VlnPlot(x.integrated, features = "percent.ribo", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_ribo.pdf"), dpi = 700)
  
  VlnPlot(x.integrated, features = "percent.hb", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_hb.pdf"), dpi = 700)
  
  VlnPlot(x.integrated, features = "nCount_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_nCount_RNA.pdf"), dpi = 700)
  
  VlnPlot(x.integrated, features = "nFeature_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_nFeature_RNA.pdf"), dpi = 700)
  
  VlnPlot(x.integrated, features = "S.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_S-Score.pdf"), dpi = 700)
  
  VlnPlot(x.integrated, features = "G2M.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_G2M-Score.pdf"), dpi = 700)
  
  # Make stacked barplot of group population per cell type to determine proportunality
  # between Control and PCOS in each cell type
  
  Idents(object = x.integrated) <- group
  pt <- table(Idents(x.integrated), x.integrated$Group_Stage)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Group") +
    ylab("Proportion") +
    theme(legend.title = element_blank()) +
    theme_cowplot()
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Stage.pdf"), dpi = 700)
  
  #Make stacked barplot of group population per cell type to determine proportunality
  # between Control and PCOS in each cell type
  
  Idents(object = x.integrated) <- group
  pt <- table(Idents(x.integrated), x.integrated$Group_Treatment)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Group") +
    ylab("Proportion") +
    theme(legend.title = element_blank()) +
    theme_cowplot()
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Treatment.pdf"), dpi = 700)
  
}

# Loading SCT Seurat
print("Seurat object loading")
x.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated.h5seurat"))
print("Setting DefaultAssay")

# Set the default assay to SCT
DefaultAssay(x.integrated) = Norm.assay #Alternative is RNA

# Setting cluster to be tested to SCT clustered
Idents(object = x.integrated) <- "seurat_clusters"

# Annotating the seurat object
# Renaming for groups and stage. 
print("Plotting Labelled_Clusters")
x.integrated$Labelled_Clusters = x.integrated$seurat_clusters
x.integrated = SetIdent(x.integrated, value = "Labelled_Clusters")

# Add the number of seurat clusters and rename them according
# to your finding
x.integrated = RenameIdents(x.integrated, "0" = "",
                               "1" = "",
                               "2" = "",
                               "3" = "",
                               "4" = "",
                               "5" = "",
                               "6" = "",
                               "7" = "",
                               "8" = "",
                               "9" = "",
                               "10" = "",
                               "11" = "",
                               "12" = "",
                               "13" = "",
                               "14" = "",
                               "15" = "",
                               "16" = "",
                               "17" = "",
                               "18" = "",
                               "19" = "",
                               "20" = "",
                               "21" = "",
                               "22" = "",
                               "23" = "",
                               "24" = "",
                               "25" = "",
                               "26" = "",
                               "27" = "",
                               "28" = "",
                               "29" = "")
x.integrated[["Labelled_Clusters"]] = Idents(object = x.integrated)
DimPlot(x.integrated, reduction = "umap", group.by = "Labelled_Clusters", 
        label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_", Norm.assay, "_UMAP_Labelled_Clusters.pdf"), dpi = 700)

# Copy the above section and rename in case additional labelling is required

# Save the labeled object
SaveH5Seurat(object = x.integrated, filename = paste0(Output.dir,Project_name, "_", Norm.assay,"_celltypes_integrated.h5seurat"), overwrite = TRUE)

# QC plotting of labelled clusters
Plot_CQ_labeled_clusters = Plotting_QC(x = x.integrated, group = "Labelled_Clusters")
Plot_CQ_labeled_clusters = Plotting_QC(x = x.integrated, group = "seurat_clusters")
Plot_CQ_labeled_clusters = Plotting_QC(x = x.integrated, group = "Group_Treatment")

# Setting cluster to be tested to SCT clustered
Idents(object = x.integrated) <- "Labelled_Clusters"

# Log2 normalising of new labels
DefaultAssay(x.integrated) = "RNA"
x.integrated = NormalizeData(x.integrated)
x.integrated = FindVariableFeatures(x.integrated)
all.genes = rownames(x.integrated)
x.integrated = ScaleData(x.integrated, features = all.genes)

# Save labeled and log2 normalised object
SaveH5Seurat(object = x.integrated, 
             filename = paste0(Output.dir,Project_name, "_log2_celltypes_integrated.h5seurat"), 
             overwrite = TRUE)

# Find Markers on seurat clusters and labelled clusters
if (FindMarker.flag == TRUE) {
  
  # Setting cluster to be tested to SCT clustered
  Idents(object = x.integrated) <- "seurat_clusters"
  
  print("Finding markers")
  x.markers <-FindAllMarkers(x.integrated, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.5)
  print("Saving markers as .csv")
  write.csv(x.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(x.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_x_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top5 = x.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  x.heatmap = DoHeatmap(x.integrated, features = top5$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top5_Genes_Heatmap.pdf"),
          plot = x.heatmap,
          dpi = 700)

  # Setting cluster to be tested to SCT clustered
  Idents(object = x.integrated) <- "Labelled_Clusters"

  print("Finding markers")
  x.markers <-FindAllMarkers(x.integrated, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.5)
  print("Saving markers as .csv")
  write.csv(x.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers_labelled.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(x.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_x_markers_labelled.rds"))
  print("Extracting top 10 markers per cluster")
  top5 = x.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


  print("Generating heatmap")
  x.heatmap = DoHeatmap(x.integrated, features = top5$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top5_Genes_Heatmap_labelled.pdf"),
          plot = x.heatmap,
          dpi = 700)

}
