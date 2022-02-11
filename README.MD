# QCR

## Quality Control for scRNA-seq data

QCR is a workflow for adding several quality flags to scRNA-seq cells or nuceli which can later be used for downstream processing. It consists of three parts, pre-processing, automated workflow and interactive workflow.

## Pre-processing:

Input for this step includes a seurat object and reference data set from which cell labels are transferred to the different cells. This is a required step as some flags are based on thecell label transfer. It also includes an optional step of integration where two or more samples are integrated using Seurat function for joint analyses.


Example:
```
Usage: Rscript <normalize_and_label_transfer_script> <work-dir location> <seurat-object in RDS/Rdata> <reference-obj in RDS/Rdata>
Rscript ./seurat/scripts/normalize_and_label_transfer_objects.R /allen/programs/celltypes/workgroups/hct/SEA-AD/RNAseq/QCR/seurat/ ./seurat/data/H2033038.RData /allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/210426_QCR/data/reference_subset.RDS
```

## Automated workflow:

This happens using Snakemake which take the Seurat object as input and adds several automated QC flags to the metadata of the Seurat object at cell and metacell (group of cells using Louvain clustering in Seurat) level.

Example:
```

Usage: sh <config_file_maker.sh> <path to *.RDS seurat data>
sh ./snakemake/config_script/config_yaml.sh /allen/programs/celltypes/workgroups/hct/SEA-AD/RNAseq/QCR/seurat/data/ > config.yaml
Fill in config file values
snakemake -np --use-conda --cores 30 (dry run)
snakemake -p  --use-conda --cores 30 (actual run)

```

## Interactive workflow:

This is realized using Shiny for visualizing the UMAPs and scatter plots from Snakemake workflow, and also add user defined flags to cells in the Seurat object by brushing

```
COMING SOON ....

```
## License
The license for this package is available on Github at: https://github.com/AllenInstitute/QCR_HVS/blob/master/LICENSE

## Level of Support
We are planning on occasional updating this repo with no fixed schedule, but likely several times per year. Community involvement is encouraged through both issues and pull requests.



