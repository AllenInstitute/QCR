# Shiny QC
An interactive tool for final steps of QC post pipeline

# Usage

```{R}
library(shinyQC)
shinyQC("/location/of/seurat/RDS/objects")
```

# Connecting to a Server from a local machine
If running on the cluster then open a seperate terminal and tunnel in:

```{bash}
ssh -t AIBS.id@cluster.name -L 5000:localhost:5000
```

Then you should be able to open a browser and interact with the server at:

```
http://127.0.0.1:5000
```
