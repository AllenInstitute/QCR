# Integration of seurat objects 
# The script takes -> run dir,
#                  -> labelled Seurat objects dir,
#                  -> final integrated object name as input 
# and returns integrated objects

# load libraries 
library(data.table)
library(Rdpack)
library(Seurat)
args = commandArgs(trailingOnly = TRUE)


# functions
# Function for integration
Integrated_object <- function(seurat_object_list, final_object){
  
  object_list  <- lapply(seurat_object_list, readRDS)
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = object_list)
  
  # Can uncomment if need more space or run takes too long
  #object_list <- lapply(X = object_list, FUN = function(x) {
  #  x <- DietSeurat(x, counts = F, data = T, scale.data = T)
  #})
  
  object_list_anchors <- FindIntegrationAnchors(object.list = object_list, anchor.features = features)

  objects_combined <- IntegrateData(anchorset = object_list_anchors)

  #DefaultAssay(objects_combined) <- "integrated"
  objects_combined@assays$RNA <- NULL
  
  # PCA and UMAP
  objects_combined <- ScaleData(objects_combined, verbose = FALSE)
  objects_combined <- RunPCA(objects_combined, npcs = 30, verbose = FALSE)
  objects_combined <- RunUMAP(objects_combined, reduction = "pca", dims = 1:30)

  saveRDS(objects_combined, file = final_object)

}

# inputs
run_dir <- args[1]
setwd(run_dir)
objs_dir <- args[2]
seurat_object_list <- Sys.glob(paste0(args[2],"*labelled.RDS", collapse = ""))

#output name 
final_object <- args[3]

# Actual run
Integrated_object(seurat_object_list,final_object)
