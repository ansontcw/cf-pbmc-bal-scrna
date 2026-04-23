#!/usr/bin/env Rscript

# This script performs HLCA v2 annotation in R v.4.2.2

# load packages ----------------------------------------------------------------
suppressPackageStartupMessages({
library(Seurat)
library(Azimuth)
library(qs)
})
set.seed(1990)

# parse arguments --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: script.R <input_qs_file> <output_qs_file>")
}
input_file <- args[1]
output_file <- args[2]

if(!file.exists(output_file)) {
  # read object ------------------------------------------------------------------
  seu <- qread(input_file, nthreads=16)
  
  # run Azimuth ------------------------------------------------------------------
  seu <- RunAzimuth(seu,
                    "/group/canc2/anson/reference/annotation/HLCA_v2",
                    umap.name = 'azimuth.umap',
                    verbose = TRUE,
                    assay = "SCT",
                    k.weight = 50,
                    n.trees = 20,
                    mapping.score.k = 100)
  
  # write output - ---------------------------------------------------------------
  # save object
  qsave(seu, output_file, nthreads=16)
}

print("HLCA v2 annotation Done!")
