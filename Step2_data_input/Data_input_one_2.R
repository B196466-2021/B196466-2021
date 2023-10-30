#! /usr/bin/Rscript

# Load required packages
library(Seurat)

# Get input arguments
args <- commandArgs(trailingOnly=TRUE)
data_dir <- args[1] # Directory containing the 10x data
type <- args[2] # Type of sample

# Read 10x data
data1 <- Read10X(data.dir = data_dir)

# Create Seurat object and set sample type
result1 <- CreateSeuratObject(counts = data1, project = type)
result1$type <- type

# Normalize data and find variable features
compare <- result1
rm(data1,result1)
gc()
compare <- NormalizeData(compare)
compare <- FindVariableFeatures(compare, selection.method = "vst", nfeatures = 2000)