#' Define the object to hold QCR flags
#'
#' @param rnaseq.data Seurat object
#'
#' @return QCR flag data.frame
#' @keywords internal
build.qcr.object = function(rnaseq.data) {
    ## Define QC flag object
    qcr.flag.names = c("qcr.cluster.class",
                      "qcr.gene.umi.high",
                      "qcr.gene.umi.low",
                      "qcr.class.label.conf",
                      "qcr.subclass.label.conf",
                      "qcr.exclude",
                      "qcr.exclude2",
                      "qcr.total")
    qcr.flags = data.frame(matrix(nrow=ncol(rnaseq.data), ncol=length(qcr.flag.names), dimnames=list(colnames(rnaseq.data), qcr.flag.names)))
    return(qcr.flags)
}

#' Seurat clustering
#'
#' This function takes a seurat object and performs dim. reduction and clustering.
#'
#' @param rnaseq.data Seurat object
#' @param resolution.value Seurat clustering resolution parameter
#'
#' @return Seurat object with PCA, UMAP, clustering
#' @keywords internal
seurat_cluster = function(rnaseq.data, resolution.value){

  ## Perform normalization
  rnaseq.data <- ScaleData(rnaseq.data)
  rnaseq.data <- RunPCA(rnaseq.data, verbose = FALSE)
  rnaseq.data <- RunUMAP(rnaseq.data, dims = 1:30, verbose = FALSE)

  ## Cluster
  rnaseq.data <- FindNeighbors(rnaseq.data, dims = 1:30)
  rnaseq.data <- FindClusters(rnaseq.data, resolution = resolution.value)
  return(rnaseq.data)
}

#' Compute mode for vector
#'
#' @param x vector of characters 
#'
#' @return Most frequencly occurance unique element
#' @keywords internal
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
