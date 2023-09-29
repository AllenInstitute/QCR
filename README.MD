# QCR

## Quality Control for scRNA-seq data

QCR is a workflow for adding several Quality Control (QC) flags to scRNA-seq cells or nuclei which can be used for downstream processing. It requires normalized and label transferred cells as a Seurat object for input.

## QCR R package:

Below is the QCR R package, a follow up function to Seurat functions to add QC flags to scRNA cells 

### Example:
```
library(QCR)
library(Seurat)

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

## QCR snakemake pipeline:

Snakemake workflow which take the normalized and label transferred seurat object as input and adds several automated QC flags to the metadata of the Seurat object at cell and metacell (group of cells using Louvain clustering in Seurat) level.

### Config file:
```
cd snakemake/config_script/  

Usage: sh <config_file_maker.sh> <path to normalized and labelled seurat data (*.RDS)> >  config.yaml

Example: sh snakemake/config_script/config_yaml.sh /your-data-location/ > config.yaml
```
### Snakemake run:

```
1. Fill in config file values 
Note: example config file with default values for all keys except samples, workdir and reference can be found under "snakemake/config_script/config_default.yaml"
2. snakemake -np --use-conda --cores 30 (dry run)
3. snakemake -p  --use-conda --cores 30 (actual run)
```

## License
The license for this package is available on Github at: https://github.com/AllenInstitute/QCR/blob/master/LICENSE

## Level of Support
We are planning on occasionally updating this repo with no fixed schedule, but likely several times per year. Community involvement is encouraged through both issues and pull requests.



