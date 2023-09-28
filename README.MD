# QCR

## Quality Control for scRNA-seq data

QCR is a workflow for adding several Quality Control (QC) flags to scRNA-seq cells or nuclei which can be used for downstream processing. It consists of three parts, pre-processing, automated workflow and interactive workflow.

![plot](./QCR_workflow.png)

## Pre-processing:

Performs normalization of seurat object/query object and transfers subclass label from a reference object to seurat object cells. This is a required step as some flags are based on the cell label transfer.


```
cd seurat/

Usage: Rscript <normalize_and_label_transfer_script> <work-dir location> <seurat-object in RDS/Rdata location> <reference-obj in RDS/Rdatalocation>

Example: Rscript ./seurat/scripts/normalize_and_label_transfer_objects.R /allen/programs/celltypes/workgroups/hct/SEA-AD/RNAseq/QCR/seurat/ ./seurat/data/H2033038.RData /allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/210426_QCR/data/reference_subset.RDS
```

## Automated QC:

Snakemake workflow which take the normalized and label transferred seurat object as input and adds several automated QC flags to the metadata of the Seurat object at cell and metacell (group of cells using Louvain clustering in Seurat) level.

### Config file:
```
cd ../snakemake/config_script/  

Usage: sh <config_file_maker.sh> <path to normalized and labelled seurat data (*.RDS)> >  ../config.yaml

Example: sh ./snakemake/config_script/config_yaml.sh /allen/programs/celltypes/workgroups/hct/SEA-AD/RNAseq/QCR/seurat/data/ > ../config.yaml
```
### Snakemake run:

```
1. cd ../
2. Fill in config file values 
Note: example config file with default values for all keys except samples, workdir and reference can be found under ./snakemake/config_script/config_default.yaml
3. snakemake -np --use-conda --cores 30 (dry run)
4. snakemake -p  --use-conda --cores 30 (actual run)
```

## Interactive QC:

Shiny app for visualizing the UMAPs and scatter plots from Snakemake workflow, and also add user defined flags to cells in the Seurat object by brushing

```
COMING SOON ...
```

## QCR R package:

A follow up function to Seurat functions to add QC flags to scRNA cells when analyzing data with Seurat

### Example:
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
## License
The license for this package is available on Github at: https://github.com/AllenInstitute/QCR_HVS/blob/master/LICENSE

## Level of Support
We are planning on occasional updating this repo with no fixed schedule, but likely several times per year. Community involvement is encouraged through both issues and pull requests.



