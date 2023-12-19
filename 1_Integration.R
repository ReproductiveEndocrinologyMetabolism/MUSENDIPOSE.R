#!/usr/bin/env Rscript

print("Loading packages")

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(future)

### Flags and directories to be set ###
LogNorm.flag = FALSE
SCT.flag = TRUE
sample.reference.flag = FALSE
do.integration.flag = TRUE
load.anchors = FALSE
integration.feat = 3000

#### Do not forget do change names in line 215 #####

Project_name = "Adipo_All"
Tissue.type = "Muscle" #Choose between Endometrium, Adipose or Muscle
sample.integration = list.files("Output/0_Pre-processed/")
sample.references = c(1:5) # These are the controls
Input.dir = "Output/0_Pre-processed/"
anchor.dir = paste0("Output/1_Integrated/", Project_name, "_anchors.RDS")

##########

if (dir.exists(path = paste0("Output/1_Integrated/")) == FALSE) {
  print(paste0("Generating output directory ", "1_Integrated"))
  dir.create(path = paste0("Output/1_Integrated/"), recursive = TRUE)
  Output.dir = paste0("Output/1_Integrated/")
} else if (dir.exists(path = paste0("Output/1_Integrated/")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/1_Integrated/")
} else {
  print("Error with output directory")
}

### Generating the Seurat object ###

# Load integrated seurat object if it exists
if (length(list.files(Output.dir, pattern = "integrated.h5seurat", full.names = TRUE)) == 1) {
  print(paste0("Loading ", list.files(Output.dir, pattern = "_integrated.h5seurat", full.names = FALSE)))
  x.integrated = LoadH5Seurat(file = list.files(Output.dir, pattern = "_integrated.h5seurat", full.names = TRUE))
} else if (length(list.files(Output.dir, pattern = "integrated.h5seurat", full.names = TRUE)) == 0 & load.anchors == FALSE) {
  # Perform integration if the integration has not been performed before
  print("Generating integrated seurat object")
  
  # Loop to load sample seurat objects
  seurat.dir = c()
  for (x in sample.integration) {
    sample.dir = x
    seurat.path = paste0(Input.dir, sample.dir, "/", sample.dir, "_Filtered_", Tissue.type,".h5seurat")
    seurat.dir = c(seurat.dir, seurat.path)
  }
  seurat.list = lapply(seurat.dir, function(x) LoadH5Seurat(file = x))
  
} else if (load.anchors == TRUE) {
  print("Anchors are to be loaded!")
}

if (LogNorm.flag == TRUE & load.anchors == FALSE) {
  print("Log Normalised")
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = TRUE)
    x <- FindVariableFeatures(x, verbose = TRUE)
  })
  # Finding integration anchors and integrating
  if (sample.reference.flag == TRUE) {
    print("Using reference")
    x.anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                        normalization.method = "LogNormalize",
                                        reference = sample.references, 
                                        reduction = "rpca", dims = 1:50)
    
    # The anchor dataset is saved so that it can be used on the cluster
    saveRDS(x.anchors, paste0(Output.dir,Project_name, "_anchors.RDS"))
    
    
  } else if (sample.reference.flag == FALSE) {
    print("Using reference")
    x.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                        reduction = "rpca", dims = 1:50)
    x.integrated <- IntegrateData(anchorset = x.anchors, dims = 1:50)
    x.integrated <- ScaleData(x.integrated, verbose = FALSE)
    x.integrated <- RunPCA(x.integrated, verbose = FALSE)
    x.integrated <- RunUMAP(x.integrated, dims = 1:50)
  }
  
  # The anchor dataset is saved
  saveRDS(x.anchors, paste0(Output.dir,Project_name, "_anchors.RDS"))
  
  # The data is integrated and saved
  if (do.integration.flag == TRUE) {
    x.integrated <- IntegrateData(anchorset = x.anchors, normalization.method = "LogNormalize", dims = 1:50)
    x.integrated <- ScaleData(x.integrated, verbose = FALSE)
    x.integrated <- RunPCA(x.integrated, verbose = FALSE)
    x.integrated <- RunUMAP(x.integrated, dims = 1:50)
    SaveH5Seurat(x.integrated, paste0(Output.dir,Project_name, "_integrated.h5seurat"), overwrite = TRUE)
  }
  
} else if (SCT.flag == TRUE & load.anchors == FALSE) {
  print("Integration of SCT transformed object")
  # Preparing SC transformed objects for integration
  #seurat.list <- lapply(X = seurat.list, FUN = SCTransform) ### They have already been normalised SCTransformed
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = integration.feat)
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
  seurat.list <- lapply(X = seurat.list, FUN = RunPCA, features = features)
  
  # Finding integration anchors and integrating
  if (sample.reference.flag == TRUE) {
    print("Using reference")
    x.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                        normalization.method = "SCT", 
                                        reference = sample.references, 
                                        anchor.features = features,
                                        dims = 1:30,
                                        reduction = "rpca", 
                                        k.anchor = 5)
  } else if (sample.reference.flag == FALSE) {
    x.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                        normalization.method = "SCT", 
                                        anchor.features = features,
                                        dims = 1:30,
                                        reduction = "rpca", 
                                        k.anchor = 5)
  }
  
  # The anchor dataset is saved
  saveRDS(x.anchors, paste0(Output.dir,Project_name, "_anchors.RDS"))
  
  # The data is integrated and saved
  if (do.integration.flag == TRUE) {
    x.integrated <- IntegrateData(anchorset = x.anchors, normalization.method = "SCT", dims = 1:30)
    x.integrated <- RunPCA(x.integrated, verbose = TRUE)
    x.integrated <- RunUMAP(x.integrated, reduction = "pca", dims = 1:30)
    x.integrated <- FindNeighbors(x.integrated, dims = 1:30)
    x.integrated <- FindClusters(x.integrated)
  }
  
} else if (load.anchors == TRUE & do.integration.flag == TRUE) {
  print("Loading the anchors!")
  x.anchors = readRDS(file = anchor.dir)
  
  # The data is integrated and saved
  if (do.integration.flag == TRUE) {
    
    if (SCT.flag == TRUE) {
      print("Integration with SCT and loaded anchors")
      x.integrated <- IntegrateData(anchorset = x.anchors, normalization.method = "SCT", dims = 1:30)
      x.integrated <- RunPCA(x.integrated, verbose = TRUE)
      x.integrated <- RunUMAP(x.integrated, reduction = "pca", dims = 1:30)
      x.integrated <- FindNeighbors(x.integrated, dims = 1:30)
      x.integrated <- FindClusters(x.integrated)
    } else if (LogNorm.flag == TRUE) {
      print("Integration with log2 normalisation and loaded anchors")
      x.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                          reduction = "rpca", dims = 1:50)
      x.integrated <- IntegrateData(anchorset = x.anchors, dims = 1:50)
      x.integrated <- ScaleData(x.integrated, verbose = FALSE)
      x.integrated <- RunPCA(x.integrated, verbose = FALSE)
      x.integrated <- RunUMAP(x.integrated, dims = 1:50)
    }
    SaveH5Seurat(object = x.integrated, filename = paste0(Output.dir,Project_name,"_integrated.h5seurat"), overwrite = TRUE)
  }
  
  # Doing log2 normalisation and clustering on SCT integrated object
  DefaultAssay(x.integrated) = "RNA"
  x.integrated = NormalizeData(object = x.integrated, normalization.method = "LogNormalize")
  x.integrated <- FindVariableFeatures(x.integrated, selection.method = "vst", nfeatures = 2000)
  
  # Scaling the data
  all.genes <- rownames(x.integrated)
  x.integrated <- ScaleData(x.integrated, features = all.genes)
  
  # Perform linear dimensional reduction
  x.integrated <- RunPCA(x.integrated, features = VariableFeatures(object = x.integrated))
  
  # Clustering the data 
  x.integrated = FindNeighbors(x.integrated, dims = 1:10)
  x.integrated = FindClusters(x.integrated, resolution = 0.7) # Resolution increases with more cells
  
  # Run non-linear dimensional reduction with UMAP
  x.integrated = RunUMAP(x.integrated, dims = 1:10)
  
  # Saving the seurat object
  SaveH5Seurat(object = x.integrated, filename = paste0(Output.dir,Project_name,"_integrated.h5seurat"), overwrite = TRUE)
  
}

if (SCT.flag == TRUE) {
  Norm.assay = "SCT"
} else if (LogNorm.flag == TRUE) {
  Norm.assay = "Log2"
}


# Plotting the integrated object
print("Plotting of UMAP")
UMAP_Ident = DimPlot(x.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Ident = DimPlot(x.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Phase = DimPlot(x.integrated, reduction = "umap", group.by = "Phase")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_UMAP_Phase.pdf"), dpi = 700)
UMAP_Label = DimPlot(x.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_labelled_UMAP_Clusters.pdf"), dpi = 700)
print("Done with basic UMAP")

# Renaming for groups and stage
x.integrated$Group_Stage = x.integrated$orig.ident
x.integrated = SetIdent(x.integrated, value = "Group_Stage")
x.integrated = RenameIdents(x.integrated, 'ctrl_204' = "Control",
                            "ctrl_209" = "Control",
                            "ctrl_210" = "Control",
                            "ctrl_211" = "Control",
                            "ls_003_w0" = "PCOS_W0",
                            "ls_003_w16" = "PCOS_W16_LS",
                            "ls_004_w0" = "PCOS_W0",
                            "ls_004_w16" = "PCOS_W16_LS",
                            "ls_027_w0" = "PCOS_W0",
                            "ls_027_w16" = "PCOS_W16_LS",
                            "met_005_w0" = "PCOS_W0",
                            "met_005_w16" = "PCOS_W16_Met",
                            "met_031_w0" = "PCOS_W0",
                            "met_031_w16" = "PCOS_W16_Met",
                            "met_047_w0" = "PCOS_W0",
                            "met_047_w16" = "PCOS_W16_Met",
                            "met_054_w0" = "PCOS_W0",
                            "met_054_w16" = "PCOS_W16_Met",
                            "met_061_w0" = "PCOS_W0",
                            "met_061_w16" = "PCOS_W16_Met",
                            "met_065_w0" = "PCOS_W0",
                            "met_073_w0" = "PCOS_W0",
                            "met_073_w16" = "PCOS_W16_Met")
x.integrated[["Group_Stage"]] = Idents(object = x.integrated)
DimPlot(x.integrated, reduction = "umap", group.by = "Group_Stage")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_UMAP_Group_Stage.pdf"), dpi = 700)

# Renaming for groups and stage
DimPlot(x.integrated, reduction = "umap", group.by = "orig.ident")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_UMAP_Pat-nr.pdf"), dpi = 700)

# Log2 normalising 
DefaultAssay(x.integrated) = "RNA"
x.integrated = NormalizeData(x.integrated)
x.integrated = FindVariableFeatures(x.integrated)
all.genes = rownames(x.integrated)
x.integrated = ScaleData(x.integrated, features = all.genes)

# Save labelled seurat object and overwrite old object
SaveH5Seurat(object = x.integrated, filename = paste0(Output.dir,Project_name, "_", Norm.assay,"_labeled_integrated.h5seurat"), overwrite = TRUE)

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
  
  # Make stacked barplot of group population per cell type to determine proportionality
  # between Control and PCOS in each cell type and after treatment
  
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
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                     vjust = 1, hjust = 1))
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Stage.pdf"), dpi = 700)
  
  # Make stacked barplot of group population per cell type to determine proportionality
  # between all different samples
  
  Idents(object = x.integrated) <- group
  pt <- table(Idents(x.integrated), x.integrated$orig.ident)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Sample") +
    ylab("Proportion") +
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                     vjust = 1, hjust = 1))
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Sample.pdf"), dpi = 700)
  
}

Plot_CQ_Group.Stage = Plotting_QC(x = x.integrated, group = "Group_Stage")
Plot_CQ_Pat_nr = Plotting_QC(x = x.integrated, group = "orig.ident")
Plot_CQ_seurat_clusters = Plotting_QC(x = x.integrated, group = "seurat_clusters")

# Make stacked barplot of group population per cell type to determine proportionality
# between Control and PCOS in each cell type and after treatment

Proportion_barplot <- function(x_group = "") {
  Idents(object = x.integrated) <- group
  pt <- table(Idents(x.integrated), x.integrated[[x_group]])
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Sample") +
    ylab("Proportion") +
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                     vjust = 1, hjust = 1))
  ggsave2(plot = barplot_prop, filename = paste0(Output.dir,Project_name, "_", Norm.assay, "_", "_barplot_proportions_", x_group, ".pdf"), dpi = 700)
  
}

barplot_prop_QC = Proportion_barplot(x_group = "Group_Stage")
barplot_prop_QC = Proportion_barplot(x_group = "orig.ident")
