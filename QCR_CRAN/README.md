# QCR (R package)

## Example:
```
## Standard Seurat workflow
rnaseq.data = CreateSeuratObject(counts = mat, meta.data = meta.data)
rnaseq.data = NormalizeData(rnaseq.data, normalization.method = "LogNormalize", scale.factor = 1e6)
rnaseq.data = FindVariableFeatures(rnaseq.data, selection.method = "vst", nfeatures = 2000)
rnaseq.data = ScaleData(rnaseq.data, features = rownames(rnaseq.data))
rnaseq.data = RunPCA(rnaseq.data, features = VariableFeatures(object = rnaseq.data))
rnaseq.data = FindNeighbors(rnaseq.data, dims = 1:30)
rnaseq.data = FindClusters(rnaseq.data, resolution = 0.5)
rnaseq.data = RunUMAP(rnaseq.data, dims = 1:30)
## Add QCR flags!
rnaseq.data = QCR(rnaseq.data=rnaseq.data,                         ## A Seurat object
                  class.col="class",                               ## Column name for class level annotation
                  class.cutoff=0.75,                               ## Percent of class homogeneity within each cluster to determine low-quality
                  neuron.class.id="Neuron",                        ## Within `class.col` the value for neurons
                  neuron.cutoff=0.05,                              ## Percent of neuron clusters to remove based on gene & umi counts
                  non.neuron.class.id="NonNeuron",                 ## Within `class.col` the value for non-neurons
                  non.neuron.cutoff=0.05,                          ## Percent of non.neuron clusters to remove based on gene & umi counts
                  class.col.score="class.prediction.score",        ## Column name for class annotation likelihood
                  subclass.col.score="subclass.prediction.score",  ## Column name for subclass annotation likelihood
                  label.confidence.threshold=0.7,                  ## Likelihood value to determine low-confidence cell type annotation
                  resolution.value=10)                             ## For Seurat:FindClusters if not already run.
```
