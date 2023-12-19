#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)

# Setting variables and flag
Project_name = "x_All"
FindMarker.flag = TRUE
#Integration.reference = c(1:5)
Input.dir = "Output/3_Labeling/"
Norm.assay = "SCT"
ScaleData.flag = TRUE
integration.feat = 3000


if (dir.exists(path = paste0("Output/4_Subsetting")) == FALSE) {
  print(paste0("Generating output directory Output/3_Labeling"))
  dir.create(path = paste0("Output/4_Subsetting"), recursive = TRUE)
  Output.dir = paste0("Output/4_Subsetting/")
} else if (dir.exists(path = paste0("Output/4_Subsetting")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/4_Subsetting/")
} else {
  print("Error with output directory")
}

print(paste0("Output saved in ", Output.dir))

# Loading Seurat
print("Seurat object loading")
x.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_", Norm.assay, "_celltypes_integrated.h5seurat"))

print("Setting DefaultAssay to SCT")
# Set the default assay to SCT
DefaultAssay(x.integrated) = Norm.assay #Alternative is RNA

# Setting cluster to be tested to SCT clustered
Idents(object = x.integrated) <- "seurat_clusters"

# Function to split and re-integrate subsetted datasets
Subset.Reintegration <- function(x, celltype = "", reference = NULL, read.anchors = FALSE) {
  
  print("Seurat object split")
  # If there is a sample with very few cells, then merge and split by group
  if (min(table(x$orig.ident)) > 30) {
    seurat.list <- SplitObject(x, split.by = "orig.ident")
  } else if (min(table(x$orig.ident)) < 30) {
    seurat.list <- SplitObject(x, split.by = "Group_Stage")
  }
  
  print("Processing seurat object")
  seurat.list <- lapply(seurat.list, SCTransform)
  features  <- SelectIntegrationFeatures(seurat.list, nfeatures = integration.feat)
  seurat.list <- PrepSCTIntegration(seurat.list, anchor.features = features)
  
  if (read.anchors == TRUE) {
    print("Loading anchors")
    anchors = readRDS(file = paste0(Output.dir, Project_name, "_", celltype,"_reclustered_anchors.rds"))
  } else if (read.anchors == FALSE) {
    if (is.null(reference) == TRUE) {
      print("Finding anchors without reference")
      anchors <- FindIntegrationAnchors(seurat.list, normalization.method = Norm.assay, anchor.features = features, dims = 1:30)
      saveRDS(anchors, paste0(Output.dir, Project_name, "_", celltype,"_reclustered_anchors.rds"))
    } else if (is.null(reference) == FALSE) {
      print("Finding anchors with integration reference")
      anchors <- FindIntegrationAnchors(seurat.list, normalization.method = Norm.assay, anchor.features = features, dims = 1:30, 
                                        reference = Integration.reference)
      saveRDS(anchors, paste0(Output.dir, Project_name, "_", celltype,"_reclustered_anchors.rds"))
    }
  }
  
  print("Integrating object")
  x.integrated <- IntegrateData(anchorset = anchors, normalization.method = Norm.assay)
  print("Running PCA")
  x.integrated <- RunPCA(x.integrated)
  print("Running UMAP")
  x.integrated <- RunUMAP(x.integrated, reduction = "pca", dims = 1:10)
  print("Finding Neighbords")
  x.integrated <- FindNeighbors(x.integrated, dims = 1:30)
  print("Finding clusters")
  x.integrated <- FindClusters(x.integrated)
  
  ###print("New Seurat object created")
  ####SaveH5Seurat(object = x.integrated, filename = paste0(Output.dir, Project_name, "_", celltype,"_reclustered.h5seurat"), overwrite = TRUE)
  
  # Plotting the integrated object
  print("Plotting of UMAP")
  UMAP_Ident = DimPlot(x.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir,Project_name, "_", celltype,"_reclustered_UMAP_Ident.pdf"), dpi = 700)
  UMAP_Phase = DimPlot(x.integrated, reduction = "umap", group.by = "Phase")
  ggsave2(paste0(Output.dir,Project_name, "_", celltype,"_reclustered_UMAP_Phase.pdf"), dpi = 700)
  UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir,Project_name, "_", celltype,"_reclustered_UMAP_Label.pdf"), dpi = 700)
  UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = FALSE, repel = FALSE)
  ggsave2(paste0(Output.dir,Project_name, "_", celltype,"_reclustered_UMAP.pdf"), dpi = 700)
  UMAP_Group = DimPlot(x.integrated, reduction = "umap", group.by = "Group_Stage", label = TRUE)
  ggsave2(paste0(Output.dir,Project_name, "_", celltype,"_reclustered_UMAP_Group.pdf"), dpi = 700)
  print("Done with basic UMAP")
  
  # Data is re-scaled after subsetting as the mean and SD will have changed 
  if (ScaleData.flag == TRUE) {
    print("Setting DefaultAssay to RNA and log2 normalising it")
    DefaultAssay(x.integrated) = "RNA"
    x.integrated = NormalizeData(x.integrated)
    x.integrated = FindVariableFeatures(x.integrated)
    all.genes = rownames(x.integrated)
    x.integrated = ScaleData(x.integrated, features = all.genes)
    
  } else if (ScaleData.flag == FALSE) {
    print("Setting DefaultAssay to SCT")
    DefaultAssay(x.integrated) = Norm.assay #Alternative is RNA
  }
  
  if (FindMarker.flag == TRUE) {
    
    # Switching to RNA assay to identify most variable genes and to
    # select markers
    
    DefaultAssay(x.integrated) = "RNA"
    x.integrated = NormalizeData(x.integrated)
    x.integrated = FindVariableFeatures(x.integrated)
    all.genes = rownames(x.integrated)
    x.integrated = ScaleData(x.integrated, features = all.genes)
    
    # Setting cluster to be tested to SCT clustered
    Idents(object = x.integrated) <- "seurat_clusters"
    
    # Running FindAllMarkers and saving output
    print("Finding markers")
    all.markers <- FindAllMarkers(x.integrated, min.pct = 0.25, logfc.threshold = 0.5) 
    
    print("Saving markers as .csv")
    write.csv(all.markers, paste0(Output.dir,Project_name, "_", celltype,"_FindAllMarkers.csv"), quote = F)
    print("Saving markers as .rds")
    saveRDS(all.markers, paste0(Output.dir,Project_name, "_", celltype,"_FindAllMarkers.rds"))
    
    # Extracting top 10 markers and saving the output
    print("Extracting top 10 markers per cluster")
    top10_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
    print("Saving top 10 markers as .csv")
    write.csv(top10_markers, paste0(Output.dir,Project_name, "_", celltype,"_FindAllMarkers_top10.csv"), quote = F)
    print("Saving top 10 markers as .rds")
    saveRDS(top10_markers, paste0(Output.dir,Project_name, "_", celltype,"_FindAllMarkers_top10.rds"))
    
  }
  
  print("New Seurat object created")
  SaveH5Seurat(object = x.integrated, filename = paste0(Output.dir, Project_name, "_", celltype,"_reclustered.h5seurat"), overwrite = TRUE)
  
}

# Set the ident to extract from
Idents(object = x.integrated) <- "Labelled_Clusters"

# Subset and reintregrate the selected cell type. Copy the section below and
# for each specific selected cell type
x.celltype = "Celltype" # Change celltype to your selected celltype
print(paste0("Subsetting and creating new object with ", x.celltype))
x.subset <- subset(x.integrated,  idents = x.celltype)
SaveH5Seurat(object = x.subset, filename = paste0(Output.dir, Project_name, "_", x.celltype, ".h5seurat"), overwrite = TRUE)
DimPlot(endo.subset, reduction = "umap", group.by = "seurat_clusters")
ggsave2(paste0(Output.dir, Project_name,"_", x.celltype, "_UMAP.pdf"), dpi = 700)
run.reintegration = Subset.Reintegration(x = x.subset, celltype = x.celltype, reference = NULL, read.anchors = FALSE)
