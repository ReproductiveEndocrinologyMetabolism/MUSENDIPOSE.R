#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)

# Setting variables and flags
Project_name = ""
FindMarker.flag = TRUE
Norm.assay = "SCT"
Input.dir = "Output/1_Integrated/"
log2.clustering = FALSE
ScaleData.flag = TRUE # Rescale data after subsetting it
Ident.group = "seurat_clusters"
SingleR.flag == FALSE

# Marker list
# Updated cell cluster markers
Epithelium_markers = c("SCGB2A2", "CPM", "ESR1", "HMGB2", "AR", "PGR", "EPCAM", "KRT8", "LGR5")
Cell_markers = c("")

# Setting directories
if (dir.exists(path = paste0("Output/2_Clustering")) == FALSE) {
  print(paste0("Generating output directory Output/2_Clustering"))
  dir.create(path = paste0("Output/2_Clustering"), recursive = TRUE)
  Output.dir = paste0("Output/2_Clustering/")
} else if (dir.exists(path = paste0("Output/2_Clustering")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/2_Clustering/")
} else {
  print("Error with output directory")
}

# Plotting function
Marker_plotting <- function(marker, name, group = Ident.group, selected.assay = Norm.assay) {
  
  print(paste("Plotting", name, "with", selected.assay))
  
  DefaultAssay(x.integrated) = selected.assay
  VlnPlot(x.integrated, features = marker, pt.size = 0, sort = "increasing", group.by = group)
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_ViolinPlot.pdf"), dpi = 700)
  
  FeaturePlot(x.integrated, features = marker)
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_FeaturePlot.pdf"), dpi = 700)
  
  DotPlot(object = x.integrated, features = marker, group.by = group) + RotatedAxis()
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_Dotplot.pdf"), dpi = 700)
}

# Loading Seurat
print("Seurat object loading")
x.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated.h5seurat"))

# Setting cluster to be tested to SCT clustered
DefaultAssay(x.integrated) = Norm.assay
Idents(object = x.integrated) <- Ident.group

#Plotting SCT clustering UMAP
print("Plotting of UMAP")
UMAP_Ident = DimPlot(x.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Phase = DimPlot(x.integrated, reduction = "umap", group.by = "Phase")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Phase.pdf"), dpi = 700)
UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Clusters.pdf"), dpi = 700)
UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = FALSE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Clusters.pdf"), dpi = 700)
UMAP_Group = DimPlot(x.integrated, reduction = "umap", group.by = "Group_Stage")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Group_Stage.pdf"), dpi = 700)
print("Done with basic UMAP of log2 normalised object")

# Running marker plotting. Add run commands based on imputed cell markers above
run = Marker_plotting(marker = Epithelium_markers, name = "Epithelium_markers")
run = Marker_plotting(marker = Cell_markers, name = "Cell_markers")


# If true, perform log2 normalisation on the integrated object and do log2 clustering
if (log2.clustering == TRUE) {
  
  Norm.assay = "Log2"
  
  DefaultAssay(x.integrated) = "RNA"
  x.integrated = NormalizeData(object = x.integrated, normalization.method = "LogNormalize")
  x.integrated <- FindVariableFeatures(x.integrated, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(x.integrated), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(x.integrated)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  var_plot = plot1 + plot2
  ggsave2(plot = var_plot, filename = paste0(Output.dir,Project_name, "_Log2_var_plot_combined.pdf"), dpi = 700)
  
  # Scaling the data
  all.genes <- rownames(x.integrated)
  x.integrated <- ScaleData(x.integrated, features = all.genes)
  
  # Perform linear dimensional reduction
  x.integrated <- RunPCA(x.integrated, features = VariableFeatures(object = x.integrated))
  
  # Plotting the PCA 
  DimPlot(x.integrated, reduction = "pca")
  ggsave2(paste0(Output.dir,Project_name, "_Log2_PCA_plot.pdf"), dpi = 700)
  
  DimHeatmap(x.integrated, dims = 1:15, cells = 500, balanced = TRUE)
  ggsave2(paste0(Output.dir,Project_name, "_Log2_PCA_heatmap.pdf"), dpi = 700)
  
  ElbowPlot(x.integrated)
  ggsave2(paste0(Output.dir,Project_name, "_Log2_ElbowPlot.pdf"), dpi = 700)
  
  # Clustering the data 
  x.integrated = FindNeighbors(x.integrated, dims = 1:10)
  x.integrated = FindClusters(x.integrated, resolution = 0.7) # Resolution increases with more cells
  
  # Run non-linear dimensional reduction with UMAP
  x.integrated = RunUMAP(x.integrated, dims = 1:10)
  
  # Save the log2 normalised Seurat object 
  print("Saving the log2 normalised object")
  SaveH5Seurat(x.integrated, paste0(Output.dir, Project_name, "_integrated_log2.h5seurat"), overwrite = TRUE)
  
  # Setting cluster to be tested to SCT clustered
  Idents(object = x.integrated) <- "seurat_clusters"
  
  #Plotting SCT clustering UMAP
  print("Plotting of UMAP")
  UMAP_Ident = DimPlot(x.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Ident.pdf"), dpi = 700)
  UMAP_Phase = DimPlot(x.integrated, reduction = "umap", group.by = "Phase")
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Phase.pdf"), dpi = 700)
  UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Clusters.pdf"), dpi = 700)
  UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = FALSE)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Clusters.pdf"), dpi = 700)
  UMAP_Group = DimPlot(x.integrated, reduction = "umap", group.by = "Group_Stage")
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Group_Stage.pdf"), dpi = 700)
  print("Done with basic UMAP of log2 normalised object")
  
  # Save the log2 normalised Seurat object 
  print("Saving the log2 normalised object")
  SaveH5Seurat(x.integrated, paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated_log2.h5seurat"), overwrite = TRUE)
  
} 


# Run FindAllMarkers to identify markers for clusters
if (FindMarker.flag == TRUE) {
  
  # Data is re-scaled after subsetting as the mean and SD will have changed 
  if (ScaleData.flag == TRUE) {
    Norm.assay = "RNA"
    print("Setting DefaultAssay to RNA and log2 normalising it")
    DefaultAssay(x.integrated) = "RNA"
    x.integrated = NormalizeData(x.integrated)
    x.integrated = FindVariableFeatures(x.integrated)
    all.genes = rownames(x.integrated)
    x.integrated = ScaleData(x.integrated, features = all.genes)
    
  } else if (ScaleData.flag == FALSE) {
    print("Setting DefaultAssay to SCT")
    DefaultAssay(x.integrated) = "SCT" #Alternative is RNA
  }
  
  SaveH5Seurat(x.integrated, paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated.h5seurat"), overwrite = TRUE)

  print("Finding markers")
  x.markers <-FindAllMarkers(x.integrated, assay = Norm.assay, 
                                logfc.threshold = 0.5, min.pct = 0.25)
  print("Saving markers as .csv")
  write.csv(x.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(x.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_x_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top10 = x.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) # Maybe change to p-value?
  
  print("Generating heatmap")
  x.heatmap = DoHeatmap(x.integrated, features = top10$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top10_Genes_Heatmap.pdf"),
          plot = x.heatmap,
          dpi = 700)
}

if (SingleR.flag == TRUE) {
  
  library(celldex)
  library(SingleR)
  
  # Setting annotation for SingleR annotation
  SingleR.annotation.ref.Monaco = celldex::MonacoImmuneData()
  SingleR.annotation.ref.DbImmune = celldex::DatabaseImmuneCellExpressionData(cell.ont = "nonna")
  SingleR.annotation.ref.Blueprint = celldex::BlueprintEncodeData(cell.ont = "nonna")
  
  x.sce = as.SingleCellExperiment(DietSeurat(x.integrated))
  
  # Run automated annotation using several different references
  singleR.main.Monaco = SingleR(test = x.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Monaco,labels = SingleR.annotation.ref.Monaco$label.main)
  singleR.main.DbImmune = SingleR(test = x.sce, assay.type.test = 1, ref = SingleR.annotation.ref.DbImmune,labels = SingleR.annotation.ref.DbImmune$label.main)
  singleR.main.Blueprint = SingleR(test = x.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Blueprint,labels = SingleR.annotation.ref.Blueprint$label.main)
  
  singleR.fine.Monaco = SingleR(test = x.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Monaco,labels = SingleR.annotation.ref.Monaco$label.fine)
  singleR.fine.DbImmune = SingleR(test = x.sce, assay.type.test = 1, ref = SingleR.annotation.ref.DbImmune,labels = SingleR.annotation.ref.DbImmune$label.fine)
  singleR.fine.Blueprint = SingleR(test = x.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Blueprint,labels = SingleR.annotation.ref.Blueprint$label.fine)
  
  
  x.SingleR = x.integrated
  
  # Testing and plotting Monaco main
  x.SingleR$immune_labels_monaco_main = singleR.main.Monaco$pruned.labels
  Idents(x.SingleR) = "immune_labels_monaco_main"
  DimPlot(x.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_SingleR_Monaco_main.pdf"), dpi = 700)
  
  # Testing and plotting DbImmune main
  x.SingleR$immune_labels_DbImmune_main = singleR.main.DbImmune$pruned.labels
  Idents(x.SingleR) = "immune_labels_DbImmune_main"
  DimPlot(x.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_SingleR_DbImmune_main.pdf"), dpi = 700)
  
  # Testing and plotting Blueprint main
  x.SingleR$immune_labels_Blueprint_main = singleR.main.Blueprint$pruned.labels
  Idents(x.SingleR) = "immune_labels_Blueprint_main"
  DimPlot(x.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_SingleR_Blueprint_main.pdf"), dpi = 700)
}

print("Done")
