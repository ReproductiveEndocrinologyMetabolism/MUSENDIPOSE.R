#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(openxlsx)
library(slingshot)
library(tradeSeq)

# Setting up script
# Manually changed the celltypes in the seurat object by editing line 49 and 62
Server.run = FALSE
Project_name = "Endo_All_Epithelium"
seurat.dir = "Output/5_Epithelium_Clustering/"
seurat.object = "Output/5_Epithelium_Clustering/Endo_All_Epithelium_reclustered_labelled.h5seurat"
select.idents = "Epithelium_labelled"
slingshot.groups = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")

##generate.rds_flag = FALSE
##load.rds_flag = TRUE

# Generating the output directory
if (dir.exists(path = paste0("Output/10_Slingshot_trajectory")) == FALSE) {
  print(paste0("Generating output directory Output/10_Slingshot_trajectory"))
  dir.create(path = paste0("Output/10_Slingshot_trajectory"), recursive = TRUE)
  Output.dir = paste0("Output/10_Slingshot_trajectory/")
} else if (dir.exists(path = paste0("Output/10_Slingshot_trajectory")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/10_Slingshot_trajectory/")
} else {
  print("Error with output directory")
}

# On the server, subset Slingshot input to later run locally for plotting
if (Server.run == TRUE) {
  
  # Load the seurat object from the seurat directory
  seurat.x = LoadH5Seurat(file = seurat.object)
  
  # Subset the data to run cell trajectory
  generate_slingshot_input <- function(seurat.x, subset.group = NULL, subset.flag = TRUE) {
    
    # If subset.flag is true, subset the selected subset group before generating slingshot input
    if (subset.flag == TRUE & is.null(subset.group) == FALSE) {
      
      # Setting the ident for subsetting
      Idents(object = seurat.x) <- "Group_Stage"
      
      # Subsetting the object and changing output label to group
      seurat.x = subset(seurat.x, idents = subset.group)
      Output.label = paste0(Project_name, "_", subset.group)
      
    }  else if (subset.flag == FALSE & is.null(subset.group) == TRUE) {
      
      # No subset or changing of output label
      Output.label = Project_name
      
    }
    
    # Check data with UMAP
    umap.plot = DimPlot(seurat.x, group.by = select.idents)
    ggsave2(plot = umap.plot, filename = paste0(Output.dir,Output.label,"_celltype_labeled_UMAP.pdf"), dpi = 700)
    
    umap.plot = DimPlot(seurat.x, group.by = "seurat_clusters")
    ggsave2(plot = umap.plot, filename = paste0(Output.dir,Output.label,"_seurat_clusters_UMAP.pdf"), dpi = 700)
    
    # Extract the objects as separate matrices for input in slingshot
    umap.dimred = seurat.x@reductions$umap@cell.embeddings
    seurat.clustering = seurat.x$seurat_clusters
    ####### CHANGE CELLTYPE HERE ####
    seurat.celltype = seurat.x$Epithelium_labelled
    RNA.counts <- as.matrix(seurat.x@assays$RNA@counts[seurat.x@assays$RNA@var.features, ])
    
    test.count = as.matrix(seurat.x@assays$RNA@counts)
    
    # Svae the matrices in the output directory
    saveRDS(object = umap.dimred, file = paste0(Output.dir, Output.label, "_umap_dimred.rds"))
    saveRDS(object = seurat.clustering, file = paste0(Output.dir, Output.label, "_seurat_clusters.rds"))
    saveRDS(object = seurat.celltype, file = paste0(Output.dir, Output.label, "_celltype_labels.rds"))
    saveRDS(object = RNA.counts, file = paste0(Output.dir, Output.label, "_RNA_counts.rds"))
    print("Matrices generated and saved")
    
    # Run default Slingshot lineage identification
    print("Running Slingshot for clusters")
    slingshot.clusters <- slingshot(Embeddings(seurat.x, "umap"), clusterLabels = seurat.x$seurat_clusters)
    ####### CHANGE CELLTYPE HERE ####
    print("Running Slingshot for celltypes")
    slingshot.celltype <- slingshot(Embeddings(seurat.x, "umap"), clusterLabels = seurat.x$Epithelium_labelled)
    
    # Save the slingshot object in the output directory
    saveRDS(object = slingshot.clusters, file = paste0(Output.dir, Output.label, "_slingshot_clusters.rds"))
    saveRDS(object = slingshot.celltype, file = paste0(Output.dir, Output.label, "_slingshot_celltype.rds"))
    
    # Run only lineages
    print("Running getLineages for clusters")
    lineages.clusters <- getLineages(data = umap.dimred, clusterLabels = seurat.clustering)
    print("Running getLineages for celltypes")
    lineages.celltype <- getLineages(data = umap.dimred, clusterLabels = seurat.celltype)
    
    # Save the lineage object in the output directory
    saveRDS(object = lineages.clusters, file = paste0(Output.dir, Output.label, "_lineages_clusters.rds"))
    saveRDS(object = lineages.celltype, file = paste0(Output.dir, Output.label, "_lineages_celltype.rds"))
    
  }
  
  # Run apply on each subset group
  sapply(slingshot.groups, function(i) generate_slingshot_input(seurat.x, subset.group = i, subset.flag = TRUE))
  
} else if (Server.run == FALSE) {
  
  # Setting up list of files to load
  list.slingshot.input = list.files(path = Output.dir, pattern = "*.rds", full.names = FALSE)
  list.slingshot.input = list.slingshot.input[grepl("*.rds$", list.slingshot.input)]
  list.slingshot.input = list.slingshot.input[grepl(Project_name, list.slingshot.input)]
   
  # Loading the slingshot input .rds
  slingshot.groups.input <- lapply(list.slingshot.input,function(i){
    readRDS(file = paste0(Output.dir, i))
  })
  
  names(slingshot.groups.input) = list.slingshot.input
  
  # Load the separate matrices for input in slingshot
  #umap.dimred = readRDS(file = paste0(Output.dir, Project_name, "_umap_dimred.rds"))
  #seurat.clustering = readRDS(file = paste0(Output.dir, Project_name, "_seurat_clusters.rds"))
  #seurat.celltype = readRDS(file = paste0(Output.dir, Project_name, "_celltype_labels.rds"))
  #RNA.counts <- readRDS(file = paste0(Output.dir, Project_name, "_RNA_counts.rds"))
  
  # Load the generated slingshot objects
  #slingshot.clusters = readRDS(file = paste0(Output.dir, Project_name, "_slingshot_clusters.rds"))
  #slingshot.celltype = readRDS(file = paste0(Output.dir, Project_name, "_slingshot_celltype.rds"))
  
  # Define a color pallete to use
  pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
  pal.celltype = c("#DD8D61", "#D3020D", "#B95FC0", "#99A3FF", "#FAA0A1", "#FF7F00", "#65C2A5")
  pal.group = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
  
  # Run default Slingshot lineage identification
  #seurat.lineages <- getLineages(data = umap.dimred, clusterLabels = seurat.clustering)
  #summary(seurat.lineages)
  
  # Run Slingshot on the seurat objecy
  #sds <- slingshot(Embeddings(seurat.x, "umap"), clusterLabels = seurat.x$seurat_clusters, 
  #                 start.clus = 13, stretch = 1)
  
  umap.dimred = slingshot.groups.input$Endo_All_Epithelium_umap_dimred.rds
  seurat.clustering = slingshot.groups.input$Endo_All_Epithelium_seurat_clusters.rds
  seurat.celltype = slingshot.groups.input$Endo_All_Epithelium_celltype_labels.rds
  Control.line = slingshot.groups.input$Endo_All_Epithelium_Control_slingshot_clusters.rds
  PCOS.line = slingshot.groups.input$Endo_All_Epithelium_PCOS_W0_slingshot_clusters.rds
  
  Metformin.line = slingshot.groups.input$Endo_All_Epithelium_PCOS_W16_Met_slingshot_clusters.rds
  Lifestyle.line = slingshot.groups.input$Endo_All_Epithelium_PCOS_W16_LS_slingshot_clusters.rds
  
  # Plot the lineages for cluster
  par(mfrow = c(1, 2))
  plot(umap.dimred[, 1:2], col = pal[seurat.clustering], cex = 0.5, pch = 16)
  for (i in levels(seurat.clustering)) {
    text(mean(umap.dimred[seurat.clustering == i, 1]), 
         mean(umap.dimred[seurat.clustering == i, 2]), labels = i, font = 2)}
  plot(umap.dimred[, 1:2], col = pal[seurat.clustering], cex = 0.5, pch = 16)
  #lines(SlingshotDataSet(seurat.lineages), lwd=2, type = 'lineages', col = 'black')
  lines(SlingshotDataSet(Control.line), lwd=2, type = 'lineages', col = pal.group[1])
  lines(SlingshotDataSet(PCOS.line), lwd=2, type = 'lineages', col = pal.group[2])
  lines(SlingshotDataSet(Metformin.line), lwd=2, type = 'lineages', col = pal.group[3])
  lines(SlingshotDataSet(Lifestyle.line), lwd=2, type = 'lineages', col = pal.group[4])
  
  # Plot the lineages for cluster
  par(mfrow = c(1, 2))
  plot(umap.dimred[, 1:2], col = pal.celltype[seurat.celltype], cex = 0.5, pch = 16)
  for (i in levels(seurat.celltype)) {
    text(mean(umap.dimred[seurat.celltype == i, 1]), 
         mean(umap.dimred[seurat.celltype == i, 2]), labels = i, font = 2)}
  plot(umap.dimred[, 1:2], col = "lightgrey", cex = 0.5, pch = 16)
  #lines(SlingshotDataSet(seurat.lineages), lwd=2, type = 'lineages', col = 'black')
  lines(SlingshotDataSet(Control.line), lwd=2, type = 'lineages', col = pal.group[1])
  lines(SlingshotDataSet(PCOS.line), lwd=2, type = 'lineages', col = pal.group[2])
  #lines(SlingshotDataSet(Metformin.line), lwd=2, type = 'lineages', col = pal.group[3])
  #lines(SlingshotDataSet(Lifestyle.line), lwd=2, type = 'lineages', col = pal.group[4])
  
  
  
  # Plot the lineages for cluster
  par(mfrow = c(1, 2))
  plot(umap.dimred[, 1:2], col = pal[seurat.clustering], cex = 0.5, pch = 16)
  for (i in levels(seurat.clustering)) {
    text(mean(umap.dimred[seurat.clustering == i, 1]), 
         mean(umap.dimred[seurat.clustering == i, 2]), labels = i, font = 2)}
  plot(umap.dimred[, 1:2], col = pal[seurat.clustering], cex = 0.5, pch = 16)
  #lines(SlingshotDataSet(seurat.lineages), lwd=2, type = 'lineages', col = 'black')
  lines(SlingshotDataSet(Control.line), lwd=2, type = 'lineages', col = 'black')
  lines(SlingshotDataSet(PCOS.line), lwd=2, type = 'lineages', col = 'red')
  
  lines(SlingshotDataSet(Control.line), lwd=2, type = 'curves', col = 'black')
  lines(SlingshotDataSet(PCOS.line), lwd=2, type = 'curves', col = 'red')
  
  # Run tradeSeq on the slingshot object
  Control.counts = slingshot.groups.input$Endo_All_Epithelium_Control_RNA_counts.rds
  PCOS.counts = slingshot.groups.input$Endo_All_Epithelium_PCOS_W0_RNA_counts.rds
  All.counts = slingshot.groups.input$Endo_All_Epithelium_RNA_counts.rds
  dim(Control.counts)
  dim(PCOS.counts)
  dim(All.counts)
  Control.GLM = fitGAM(Control.line)
  
}
