#' Visualization code for scatter plot
#'
#' @param input_file A Seurat object with nFeature_RNA, nCount_RNA and seurat_clusters
#'
#' @return ggplot2 object
#'
#' @keywords internal
scatter.plot = function(input_file){  
  ##
  ggplot(input_file, aes(x=nFeature_RNA, y=nCount_RNA, label=seurat_clusters)) + 
          geom_point() +  
          theme_bw() + 
          NoLegend() + 
          ggtitle("Drag to select")
}

#' Visualization code for umap plot
#'
#' @param input_file A Seurat object with UMAP dim. reduction and seurat_clusters
#' @param selected The cluster ids which have been selected via brush tool on scatter plot
#' @param colour.by Seurat meta.data field to colour the umap points
#' @param range.x x-axis limits for zooming
#' @param range.y y-axis limits for zooming
#'
#' @return ggplot2 object
#'
#' @keywords internal
umap = function(input_file, selected=NULL, colour.by="cluster", range.x=NULL, range.y=NULL){

    ##
    umap.embeddings = input_file@reductions$umap@cell.embeddings
    umap.df = data.frame(x = umap.embeddings[,1],
                         y = umap.embeddings[,2],
                         cluster = input_file@meta.data$seurat_clusters)
    umap.plot = umap.df %>% 
                  ggplot(aes(x = x, y = y, colour=.data[[colour.by]], fill=.data[[colour.by]])) + 
                      geom_point(size=0.2) +
                      ylab("UMAP-2") +
                      xlab("UMAP-1") +
                      theme_bw() +
                      theme(axis.text = element_text(color = "#000000"),
                              legend.position = "right",
                              panel.background = element_rect(fill = "#FFFFFF")) +
                      guides(shape = guide_legend(override.aes = list(size = 5)),
                              fill  = guide_legend(override.aes = list(size = 5))) + NoLegend() + ggtitle("Drag then double click to zoom")

    ##
    if(!is.null(selected) & (length(selected) > 0)){ 
      umap.selected = umap.df %>% 
                        filter(cluster %in% selected)
      umap.plot = umap.plot +  
                    geom_point(data = umap.selected, aes(x=x, y=y, fill=.data[[colour.by]]), shape=21, size = 3, stroke=0.75, colour="black") 
    }

    ##
    if(!is.null(range.x)){ 
      umap.plot = umap.plot +  
                    coord_cartesian(xlim = range.x, ylim = range.y, expand = FALSE)
    }
    umap.plot
}