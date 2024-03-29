% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flags.R
\name{QC_gene_umi}
\alias{QC_gene_umi}
\title{Checks for min/max genes and umi in each cell (Prev: QC4)}
\usage{
QC_gene_umi(
  rnaseq.data,
  class.col,
  non.neuron.class.id = "Non-neuron",
  non.neuron.cutoff = 0.05,
  neuron.class.id = "Neuron",
  neuron.cutoff = 0.05
)
}
\arguments{
\item{rnaseq.data}{Seurat object}

\item{class.col}{Column of suerat object that contains \code{class} annotations}

\item{non.neuron.class.id}{Name of non-neuronal cells within the \code{class} annotations}

\item{non.neuron.cutoff}{Proportion of top/bottom clusters to flag for non-neuron clusters}

\item{neuron.class.id}{Name of neuronal cells within the \code{class} annotations}

\item{neuron.cutoff}{Proportion of top/bottom clusters to flag for neuron clusters}
}
\value{
gene.umi.flag
}
\description{
Flags the top and bottom clusters for both neuron and non-neurons to be removed based on median umi and expressed genes.
Careful with this flag as cells within gene/umi tolerances could be flagged.
}
