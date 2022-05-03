#' Checks each cluster for homogeneity of class level annotations. (Prev: QC2)
#'
#' @param rnaseq.data Seurat object
#' @param class.col Column of suerat object that contains `class` annotations
#' @param class.cutoff Clusters with less than class.cutoff homogeneity are flagged
#'
#' @return cluster.class.flag
#' @export
QC_cluster_class = function(rnaseq.data, class.col, class.cutoff=0.75){

  ## Tally the count of each class for each cluster
  cluster.class.count = rnaseq.data@meta.data %>% 
      group_by_at(c("seurat_clusters", `class.col`)) %>% 
      summarise(Num = n()) %>%
      spread(`class.col`, Num, fill=0) %>% 
      as.data.frame(stringsAsFactors=F)

  ## If a cluster has < 75% of one class then we flag it.
  cluster.class.dist = t(apply(cluster.class.count[,-1], 1, function(x) x/sum(x)))
  cluster.class.count$flag = apply(cluster.class.dist, 1 ,function(x) any(x >= class.cutoff))

  ## Now map TRUE/FALSE back to each cell from cluster annotation
  cluster.class.flag = cluster.class.count$flag[match(rnaseq.data@meta.data$seurat_clusters, cluster.class.count$seurat_clusters)]
  return(cluster.class.flag)
}

#' Checks for min/max genes and umi in each cell (Prev: QC4)
#'
#' Flags the top and bottom clusters for both neuron and non-neurons to be removed based on median umi and expressed genes. 
#' Careful with this flag as cells within gene/umi tolerances could be flagged. 
#'
#' @param rnaseq.data Seurat object
#' @param class.col Column of suerat object that contains `class` annotations
#' @param non.neuron.class.id Name of non-neuronal cells within the `class` annotations
#' @param non.neuron.cutoff Proportion of top/bottom clusters to flag for non-neuron clusters
#' @param neuron.class.id Name of neuronal cells within the `class` annotations
#' @param neuron.cutoff Proportion of top/bottom clusters to flag for neuron clusters
#'
#' @return gene.umi.flag
#' @export
QC_gene_umi = function(rnaseq.data, class.col, non.neuron.class.id="Non-neuron", non.neuron.cutoff=0.05, neuron.class.id="Neuron", neuron.cutoff=0.05){

  ## gene umi
  umi.gene.check = rnaseq.data@meta.data %>%
                  group_by(seurat_clusters) %>%
                  summarise(class = Mode(.data[[class.col]]), 
                            nCount_RNA = median(nCount_RNA), 
                            nFeature_RNA = median(nFeature_RNA), 
                            avg=(nCount_RNA + nFeature_RNA)/2) %>%
                  as.data.frame()

  ## Flag top clusters
  umi.gene.check.nn = umi.gene.check %>% 
                            filter(class == non.neuron.class.id) %>%
                            slice_max(order_by = avg, prop=non.neuron.cutoff) %>%
                            pull(seurat_clusters)

  umi.gene.check.neuron = umi.gene.check %>% 
                            filter(class == neuron.class.id) %>% 
                            slice_max(order_by = avg, prop=neuron.cutoff) %>%
                            pull(seurat_clusters)
  
  umi.gene.high.flag <- ifelse(rnaseq.data$seurat_clusters %in% c(umi.gene.check.nn, umi.gene.check.neuron), FALSE, TRUE)

  ## Flag bottom clusters
  umi.gene.check.nn = umi.gene.check %>% 
                            filter(class == non.neuron.class.id) %>%
                            slice_min(order_by = avg, prop=non.neuron.cutoff) %>%
                            pull(seurat_clusters)

  umi.gene.check.neuron = umi.gene.check %>% 
                            filter(class == neuron.class.id) %>% 
                            slice_min(order_by = avg, prop=neuron.cutoff) %>%
                            pull(seurat_clusters)

  umi.gene.low.flag <- ifelse(rnaseq.data$seurat_clusters %in% c(umi.gene.check.nn, umi.gene.check.neuron), FALSE, TRUE)

  ## Merge results
  gene.umi.flag = data.frame(umi.gene.high.flag = umi.gene.high.flag, umi.gene.low.flag = umi.gene.low.flag)

  ##
  return(gene.umi.flag)
}

#' Check for confidence in cell annotation
#'
#' @param rnaseq.data Seurat object
#' @param label.confidence Confidence value per cell
#' @param label.confidence.threshold Threshold to determine low confidence and acceptable annotations
#'
#' @return label.conf.flag
#' @export
QC_label_conf = function(rnaseq.data, label.confidence, label.confidence.threshold){
  ## Tally the count of each class for each cluster
  if(!is.null(rnaseq.data@meta.data$`label.confidence`)){
    label.conf.flag = ifelse(rnaseq.data@meta.data$`label.confidence` > label.confidence.threshold , TRUE, FALSE)
  }else{
    label.conf.flag = rep(NA, ncol(rnaseq.data))
  }
  return(label.conf.flag)
}

