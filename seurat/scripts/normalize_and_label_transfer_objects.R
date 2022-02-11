# Normalization and initial label transfer of Seurat objects
# The script takes -> run dir, 
#                  -> Seurat object as RDS or Rdata and 
#                  -> reference object as RDS or Rdata input 
# and returns a normalized subclass label transferred Seurat object as output

# load libraries
library(data.table)
library(Rdpack)
library(Seurat)
library(base)
args = commandArgs(trailingOnly = TRUE)

# functions
# function to load RData dir to a variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# function for label transfer
norm_label_transfer <- function(query_object, reference_object, final_object){
  
  # Normalize
  query_object <- NormalizeData(query_object)
  
  # label transfer
  transfer_anchors <- FindTransferAnchors(reference = reference_object, query = query_object,  dims = 1:30)
  predictions <- TransferData(anchorset = transfer_anchors, refdata = reference_object$confirmed_subclass, dims = 1:30)
  query_object <- AddMetaData(query_object, metadata = predictions)
  
  # find variable features
  query_object <- FindVariableFeatures(query_object, selection.method = "vst", nfeatures = 3000)
  
  
  # run PCA and UMAP
  query_object <- ScaleData(query_object)
  query_object <- RunPCA(query_object, verbose = FALSE)
  query_object <- RunUMAP(query_object, dims = 1:30, verbose = FALSE)
  
  # save the file
  saveRDS(query_object, file = paste0("./data/",final_object, collapse = ""))
}


# inputs
run_dir <- args[1]
setwd(run_dir)
query_file <- args[2]
reference_file <- args[3]

# output name
final_object <- paste0(strsplit(basename(query_file), split = ".", fixed = TRUE)[[1]][1],"_labelled.RDS", collapse = "")

# conditions
# file extension
if (endsWith(query_file, ".rds") | endsWith(query_file, ".RDS")){
  query_object <- readRDS(query_file)
}else{
  query_object <- loadRData(query_file)
}

# Reading the RDS or Rdata file for reference_object
if (endsWith(reference_file, ".rds") | endsWith(reference_file, ".RDS")){
  reference_object <- readRDS(reference_file)
}else{
  reference_object <- loadRData(reference_file)
}

# column name checks
if("confirmed_subclass" %in% colnames(reference_object@meta.data)){
  print("Reference object has the required columns")
}else{
  stop("Reference object is missing confirmed_subclass")
}

# Actual run
norm_label_transfer(query_object, reference_object, final_object)
