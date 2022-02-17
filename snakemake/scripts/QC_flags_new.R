#Quality control flag script 
# The script takes -> run dir
#                  -> query object
#                  -> ref object
#                  -> final object 
#                  -> ref subclass
#                  -> ref class
#                  -> sample id
# and gives a seurat object of the sample with QCs columns added 

# columns required in ref and Seurat object that can be edited in config file
# reference
## confirmed_subclass
## class2
# query
## sample_id
## donor_name
# glia_cutoff
# neuron_cutoff

# Loading all Libraries

library(data.table)
library(Rdpack)
library(Seurat)
library(tidyr)
library(scales)
library(reticulate)
suppressMessages(library(dplyr))
library(tibble)
library(feather)
args = commandArgs(trailingOnly = TRUE)

#functions

# Mode function 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# QCs columns
QCs <- function(query_object, reference_object, final_object, ref_subclass, ref_class, sample_id, donor_name,exclude, glia_cutoff, neuron_cutoff, resolution_value, feather_file){
  
    # clustering
    query_object <- ScaleData(query_object)
    query_object <- RunPCA(query_object, verbose = FALSE)
    query_object <- RunUMAP(query_object, dims = 1:30, verbose = FALSE)
    
    query_object <- FindNeighbors(query_object, dims = 1:30)
    query_object <- FindClusters(query_object, resolution = as.numeric(resolution_value))
  
  
  
  ## Adding class to object ##
  ref_df <- unique(reference_object@meta.data[, c(ref_subclass, ref_class)])
  query_df <- unique(query_object@meta.data[,c(sample_id,"predicted.id")])
  
  query_df$class <- ref_df[,ref_class][match(query_df$predicted.id, ref_df[,ref_subclass])]
  
  query_df_final <- as.data.frame(query_df[,c("class")])
  colnames(query_df_final) <- "class"
  rownames(query_df_final) <- query_df[,sample_id]
  
  query_object <- AddMetaData(query_object, query_df_final)
  
  
  ## QC2 ##
  meta_cell_class <- as.data.frame(table(query_object@meta.data[,c("seurat_clusters","class")]))
  meta_cell_class.wide <- pivot_wider(meta_cell_class, names_from = "class", values_from = "Freq")
  meta_cell_class.wide <- as.data.frame(meta_cell_class.wide)
  
  rownames(meta_cell_class.wide) <- meta_cell_class.wide[,1]
  meta_cell_class.wide <- meta_cell_class.wide[-1]
  meta_cell_class.wide.percent <- as.data.frame(t(apply(meta_cell_class.wide, 1, function(x) round((x/sum(x)) *100,2))))
  meta_cell_class.wide.percent$QC2_metacell_class_flag <-apply(meta_cell_class.wide.percent, 1 ,function(x) any(x >=75))
  
  to_add <- meta_cell_class.wide.percent[c("QC2_metacell_class_flag")]
  to_add$seurat_clusters <- rownames(to_add)
  
  query_df <- unique(query_object@meta.data[,c(sample_id,"seurat_clusters")])
  query_df_final <- query_df %>% 
    left_join(to_add)
  query_df_final_rownames <- query_df_final[,sample_id]
  
  query_df_final <- as.data.frame(query_df_final[,c("QC2_metacell_class_flag")])
  rownames(query_df_final) <- query_df_final_rownames
  colnames(query_df_final) <- "QC2_metacell_class_flag"
  
  query_object <- AddMetaData(query_object, query_df_final)
  
  
  ## QC3 ##
  
  if(length(unique(query_object@meta.data[,donor_name])) == 1){
    print("one donor object")
    query_df_final <- as.data.frame(query_object@meta.data[,sample_id])
    rownames(query_df_final) <- query_df_final[,1]
    query_df_final$QC3_metacell_donor_flag <- TRUE
    query_df_final <- query_df_final[-1]
    
    
    query_object <- AddMetaData(query_object, query_df_final)
    
    
  } else{
    
    # donor check
    
    meta_cell_donor <- as.data.frame(table(query_object@meta.data[, c("seurat_clusters",donor_name)]))
    meta_cell_donor.wide <- pivot_wider(meta_cell_donor, names_from = donor_name, values_from = "Freq")
    meta_cell_donor.wide <- as.data.frame(meta_cell_donor.wide)
    
    rownames(meta_cell_donor.wide) <- meta_cell_donor.wide[,1]
    meta_cell_donor.wide <- meta_cell_donor.wide[-1]
    meta_cell_donor.wide.percent <- as.data.frame(t(apply(meta_cell_donor.wide, 1, function(x) round((x/sum(x)) *100,2))))
    meta_cell_donor.wide.percent$QC3_metacell_donor_flag <-apply(meta_cell_donor.wide.percent, 1 ,function(x) ifelse(any(x > 95), FALSE, TRUE))
    
    to_add <- meta_cell_donor.wide.percent[c("QC3_metacell_donor_flag")]
    to_add$seurat_clusters <- rownames(to_add)
    
    query_df <- unique(query_object@meta.data[,c(sample_id,"seurat_clusters")])
    query_df_final <- query_df %>% 
      left_join(to_add)
    query_df_final_rownames <- query_df_final[,sample_id]
    
    query_df_final <- as.data.frame(query_df_final[,c("QC3_metacell_donor_flag")])
    rownames(query_df_final) <- query_df_final_rownames
    colnames(query_df_final) <- "QC3_metacell_donor_flag"
    
    
    query_object <- AddMetaData(query_object, query_df_final)
    
    
  }
  
  ## QC4 ##
  
  class_seuart_df <- query_object@meta.data[,c("seurat_clusters","class")]
  
  class_seuart_df <- as.data.frame(class_seuart_df %>% group_by(seurat_clusters) %>% mutate(class_final=Mode(class)))
  class_seuart_df_final <- unique(class_seuart_df[,c("seurat_clusters","class_final")])
  
  ## gene umi
  gene_umi_df <- query_object@meta.data[,c("seurat_clusters","nFeature_RNA","nCount_RNA")]
  umi.median.df <- aggregate(nCount_RNA  ~ seurat_clusters, gene_umi_df, median)
  gene.median.df <- aggregate(nFeature_RNA  ~ seurat_clusters, gene_umi_df, median)
  umi_gene_median_df <- umi.median.df %>% left_join(gene.median.df)
  umi_gene_median_df$avg <- (umi_gene_median_df$nCount_RNA + umi_gene_median_df$nFeature_RNA)/2
  
  # joining class and seurat_cluster
  
  umi_gene_median_df_final <- umi_gene_median_df %>% 
    left_join(class_seuart_df_final)
  
  #subset
  umi_gene_median_df_final_glia <- as_tibble(umi_gene_median_df_final[umi_gene_median_df_final$class_final == "glia",])
  umi_gene_median_df_final_neuron <- as_tibble(umi_gene_median_df_final[umi_gene_median_df_final$class_final %in% c("exc","inh"),])
  
  # putative bad meta cell
  top_glia <- as.character(top_frac(umi_gene_median_df_final_glia, as.numeric(glia_cutoff), c(umi_gene_median_df_final_glia$avg))$seurat_clusters)
  top_neuron <- as.character(top_frac(umi_gene_median_df_final_neuron, as.numeric(neuron_cutoff), c(umi_gene_median_df_final_neuron$avg))$seurat_clusters)
  
  bottom_glia <- as.character(top_frac(umi_gene_median_df_final_glia, -(as.numeric(glia_cutoff)), c(umi_gene_median_df_final_glia$avg))$seurat_clusters)
  bottom_neuron <- as.character(top_frac(umi_gene_median_df_final_neuron, -(as.numeric(neuron_cutoff)), c(umi_gene_median_df_final_neuron$avg))$seurat_clusters)
  
  # assigning true false
  gene_umi_df$QC4_metacell_gene_umi_high_flag <- ifelse(gene_umi_df$seurat_clusters %in% c(top_glia, top_neuron), FALSE, TRUE)
  gene_umi_df$QC4_metacell_gene_umi_low_flag <- ifelse(gene_umi_df$seurat_clusters %in%  c(bottom_glia, bottom_neuron), FALSE, TRUE)
  gene_umi_df <- gene_umi_df[c("QC4_metacell_gene_umi_high_flag","QC4_metacell_gene_umi_low_flag")]
  
  
  # add to object
  query_object <- AddMetaData(query_object, gene_umi_df)
  
  
  ## QC5 ##
  query_object_sub <- query_object@meta.data["prediction.score.max"]
  query_object_sub$QC5_subclass_prediction_flag <- ifelse(query_object_sub$prediction.score.max > 0.75 , TRUE, FALSE)
  query_object_final <- query_object_sub["QC5_subclass_prediction_flag"]
    
  # add to object
  query_object <- AddMetaData(query_object, query_object_final)
    
  ## QC6 ##
  query_object_sub <- query_object@meta.data[,c(exclude,"QC2_metacell_class_flag", "QC3_metacell_donor_flag","QC4_metacell_gene_umi_low_flag","QC4_metacell_gene_umi_high_flag","QC5_subclass_prediction_flag")]
  query_object_sub[,exclude] <- ifelse(query_object_sub[,exclude] == "No", TRUE, FALSE)
  
  query_object_sub <- as.data.frame(query_object_sub)
  #TRUE 0, FALSE 1
  query_object_sub2 <- apply(query_object_sub, 2, function(x) ifelse(x == TRUE, 0, 1))
  query_object_sub$QC6_total_flag <- apply(query_object_sub2, 1, sum)
  
  query_df_final <- query_object_sub[c("QC6_total_flag")]
  
  # add to object
  query_object <- AddMetaData(query_object, query_df_final)
 
  # pulling newly added columns to df
  query_obj_meta_cols <- query_object@meta.data[,c(sample_id, exclude, "QC2_metacell_class_flag", "QC3_metacell_donor_flag","QC4_metacell_gene_umi_low_flag","QC4_metacell_gene_umi_high_flag","QC5_subclass_prediction_flag", "QC6_total_flag")]   

  #saveRDS(query_object, file = final_object)
  seurat2h5ad(query_object, "RNA", "data", ".", final_object)
  write_feather(query_obj_meta_cols, feather_file)
}


## TOdo:  Validate the user inputs to the right type
# input 
run_dir <- args[1]
setwd(run_dir)
source("./scripts/seurat_h5ad.R")
query_object <- readRDS(args[2])
reference_file <- args[3]


# output
final_object <- args[4]

# col names
ref_subclass <- args[5]
ref_class <- args[6]
sample_id <- args[7]
donor_name <- args[8]
exclude <- args[9]

# parameters
glia_cutoff <- args[10]
neuron_cutoff <- args[11]
resolution_value <- args[12]

# additional output
feather_file <- args[13]
  
# conditions
# Reading the RDS or Rdata file for reference_object
if (endsWith(basename(reference_file), ".rds") | endsWith(basename(reference_file), ".RDS")){
  reference_object <- readRDS(reference_file)
}else{
  reference_object <- loadRData(reference_file)
}


# Actual run
QCs(query_object, reference_object, final_object, ref_subclass, ref_class, sample_id, donor_name, exclude, glia_cutoff, neuron_cutoff, resolution_value, feather_file)

