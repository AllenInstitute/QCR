#' Creates a summarized data.frame aggregated by seurat_cluster
#'
#' @param query.df A data.frame 
#'
#' @return data.frame
#'
#' @keywords internal
data_summary = function(query.df){
  query.df = query.df %>% group_by(seurat_clusters)

  query.df.numeric = query.df[,unlist(lapply(query.df, is.numeric))]
  query.df.numeric$seurat_clusters = as.integer(query.df$seurat_clusters)

  summary.df  = aggregate(query.df.numeric, by=list(query.df.numeric$seurat_clusters), median)
  summary.df$nFeature_RNA = log10(summary.df$nFeature_RNA)
  summary.df$nCount_RNA   = log10(summary.df$nCount_RNA)
  return(summary.df)
}