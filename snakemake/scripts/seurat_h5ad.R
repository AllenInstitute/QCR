#' Take a Seurat object and make it better
#'
#' This function will take a Seurat object and convert it to anndata format then save as an .h5ad file
#' for much more efficient loading and cirrocumulus.
#'
#' @info Nelson Johansen, 1/19/2022
#'
#' @param seurat.obj A standard Seurat object 
#' @param assay The assay to pull from Seurat object (typically RNA, but sometimes Integrated)
#' @param slot The gene expression data slot to use
#' @param save.dir Directory where the h5ad file will be saved
#' @param file.name Name of the h5ad file
#'
#' @return A matrix of the infile
#' @export
seurat2h5ad = function(seurat.obj, assay="RNA", slot="data", save.dir=".", file.name="annSeurat"){

    ## Reqs.
    library(Seurat)
    library(reticulate)
    library(Matrix)
    anndata = reticulate::import("anndata")

    ## Gather gene exp. data
    expr = GetAssayData(object = seurat.obj, assay = assay, slot=slot)
    
    ## Gather meta.data for obs, filter out issue points (rare)
    obs = seurat.obj[[]]
    for(name in names(obs)) {
        if(class(obs[[name]])=='data.frame') {
            obs[[name]] = NULL
        }
    }
    obs$ident = Idents(object = seurat.obj)

    ## Get gene names for var
    var = as.data.frame(rownames(expr)); colnames(var) = "gene"; rownames(var) = var$gene

    ## Record every dimensionality reduction in the Seurat assay, PCA, UMAP etc.
    obsm = list()
    for (dr in Seurat:::FilterObjects(object = seurat.obj, classes.keep = "DimReduc")) {
        print(paste0("Adding: ", dr))
        obsm[[dr]] = Seurat::Embeddings(object = seurat.obj[[dr]])
    }

    ## Create our anndata object
    adata = anndata$AnnData(X = t(expr), obs = obs, var=var, obsm=obsm)

    ## Write out the anndata object in h5ad format
    adata$write_h5ad(file.path(save.dir, paste0(file.name)))
}
