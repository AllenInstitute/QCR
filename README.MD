# QCR

## Quality Control for scRNA-seq data

QCR is a workflow for adding several quality flags to scRNA-seq cells or nuceli which can later be used for downstream processing. It consists of three parts, pre-processing, automated workflow and interactive workflow.

## Pre-processing:

Input for this step includes a seurat object and reference data set from which cell labels are transferred to the different cells. This is a required step as some flags are based on the cell label transfer.


```
Usage: Rscript <normalize_and_label_transfer_script> <work-dir location> <seurat-object in RDS/Rdata location> <reference-obj in RDS/Rdatalocation>

Example: Rscript ./seurat/scripts/normalize_and_label_transfer_objects.R /allen/programs/celltypes/workgroups/hct/SEA-AD/RNAseq/QCR/seurat/ ./seurat/data/H2033038.RData /allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/210426_QCR/data/reference_subset.RDS
```

## Automated workflow:

This happens using Snakemake which take the Seurat object as input and adds several automated QC flags to the metadata of the Seurat object at cell and metacell (group of cells using Louvain clustering in Seurat) level.

### Config file:
```
Usage: sh <config_file_maker.sh> <path to *.RDS normalized and labelled seurat data location>

Example: sh ./snakemake/config_script/config_yaml.sh /allen/programs/celltypes/workgroups/hct/SEA-AD/RNAseq/QCR/seurat/data/ > config.yaml
```
### Snakemake run:

```
1. Fill in config file values 
Note: example config file with default values for all keys except samples, workdir and reference can be found ./snakemake/config_script/config_default.yaml
2. snakemake -np --use-conda --cores 30 (dry run)
3. snakemake -p  --use-conda --cores 30 (actual run)
```

## Interactive workflow:

This is realized using Shiny for visualizing the UMAPs and scatter plots from Snakemake workflow, and also add user defined flags to cells in the Seurat object by brushing

```
COMING SOON ...
```

## R WRAPPER:

The R wrapper is a follow up function to Seurat functions to add QC flags to cells when analyzing data with Seurat

```
COMING SOON ...
```
## License
The license for this package is available on Github at: https://github.com/AllenInstitute/QCR_HVS/blob/master/LICENSE

## Level of Support
We are planning on occasional updating this repo with no fixed schedule, but likely several times per year. Community involvement is encouraged through both issues and pull requests.



