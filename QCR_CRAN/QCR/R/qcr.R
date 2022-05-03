#' QCR
#'
#' This function computes various QC flags for a seurat object
#'
#' @param rnaseq.data Seurat object
#' @param class.col Column of suerat object that contains `class` annotations
#' @param class.cutoff Clusters with less than class.cutoff homogeneity are flagged
#' @param neuron.class.id Name of neuronal cells within the `class` annotations
#' @param neuron.cutoff Proportion of top/bottom clusters to flag for neuron clusters
#' @param non.neuron.class.id Name of non-neuronal cells within the `class` annotations
#' @param non.neuron.cutoff Proportion of top/bottom clusters to flag for non-neuron clusters
#' @param class.col.score Confidence value per cell for class annotation 
#' @param subclass.col.score Confidence value per cell for subclass annotation 
#' @param label.confidence.threshold Threshold to determine low confidence and acceptable annotations
#' @param resolution.value Seurat clustering resolution parameter if `seurat_clusters` is not defined in rnaseq.data
#'
#' @return data.frame with QCR flags
#' @export
QCR <- function(rnaseq.data, 
                class.col, class.cutoff,
                neuron.class.id, neuron.cutoff,
                non.neuron.class.id, non.neuron.cutoff,
                class.col.score=NULL, subclass.col.score=NULL,
                label.confidence.threshold=0.7,
                resolution.value=10){
  
  ## Clustering
  if(is.null(rnaseq.data$seurat_clusters)){
    rnaseq.data = seurat_cluster(rnaseq.data, 
                                    resolution.value)
  }

  ## Define QC flag object
  qcr.flags = build.qcr.object(rnaseq.data)
  
  ## QC2: Class in clusters ##
  qcr.flags$qcr.cluster.class = QC_cluster_class(rnaseq.data, class.col, class.cutoff)
  
  ## QC4: Top and bottom clusters for neuron and non-neurons based on avg of umi & gene counts
  gene.umi.flags = QC_gene_umi(rnaseq.data, 
                              class.col, 
                              non.neuron.class.id, 
                              non.neuron.cutoff, 
                              neuron.class.id, 
                              neuron.cutoff)
  qcr.flags$qcr.gene.umi.high = gene.umi.flags$umi.gene.high.flag
  qcr.flags$qcr.gene.umi.low = gene.umi.flags$umi.gene.low.flag
  
  ## QC5: Flag cells with low confidence in labeling
  qcr.flags$qcr.class.label.conf = QC_label_conf(rnaseq.data, class.col.score, label.confidence.threshold)
  qcr.flags$qcr.subclass.label.conf = QC_label_conf(rnaseq.data, subclass.col.score, label.confidence.threshold)
    
  ## QC6: Relabel previous annotation if available
  if(!is.null(rnaseq.data$exclude)){ qcr.flags$qcr.exclude = ifelse(rnaseq.data$exclude == "No", TRUE, FALSE) }else{ qcr.flags$qcr.exclude = rep(NA, ncol(rnaseq.data)) }
  if(!is.null(rnaseq.data$exclude2)){ qcr.flags$qcr.exclude2 = ifelse(rnaseq.data$exclude2 == "No", TRUE, FALSE) }else{ qcr.flags$qcr.exclude2 = rep(NA, ncol(rnaseq.data)) }

  ## For each cell determine the number of assigned flags
  qcr.flags$qcr.total = rowSums(!qcr.flags, na.rm=T)

  ## Info for downstream user
  print("QCR Done! TRUE: Keep cell; FALSE: Throw away cell")
  return(qcr.flags)
}