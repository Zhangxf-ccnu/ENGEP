% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/engep.R
\name{gene_dataliat}
\alias{gene_dataliat}
\title{gene_dataliat}
\usage{
gene_dataliat(spa_counts, ref_list, pre_genes, n0)
}
\arguments{
\item{spa_counts}{gene expression matrix of spatial dataset (gene by cell).}

\item{ref_list}{each element in the list is a gene expression matrix (gene by cell) of
an original single reference dataset.}

\item{pre_genes}{an array contains names of genes to be predicted, if is NULL, ENGEP will
predict the intersection of unique genes of each references. If you let pre_genes = NULL,
we suggest you to use reference datasets with high variable genes.}

\item{n0}{the number of cells of sub-reference dataset, default is 8000}
}
\value{
a list contains equal-sized reference sub-datasets with common genes and predicted genes,
spatial dataset with common genes, and unique genes in references.
}
\description{
Function to randomly partition large reference into small reference.
equal-sized sub-datasets.
}
