#!/usr/bin/env Rscript

# This script performs integration and clustering analysis

# Load packages
suppressPackageStartupMessages({
library(harmony)
library(Seurat)
library(here)
library(qs)
library(ggplot2)
library(clustree)
library(cowplot)
library(forcats)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pals)
library(AnnotationHub)
library(ensembldb)
library(msigdbr)
library(Homo.sapiens)
library(scuttle)
library(Cepo)
library(patchwork)
library(speckle)
library(scCustomize)
library(speckle)
})
set.seed(1990)
options(future.globals.maxSize = 500000000000, future.seed=TRUE)
date <- Sys.Date()
celltype <- "bal"

if(!dir.exists(here("data","plots",celltype))) {
  dir.create(here("data","plots",celltype), recursive = TRUE)
}

if(!dir.exists(here("data","SCEs","merge"))) {
  dir.create(here("data","SCEs","merge"), recursive = TRUE)
}

if(!dir.exists(here("data","SCEs","integration"))) {
  dir.create(here("data","SCEs","integration"), recursive = TRUE)
}

# SCTransform, merge, and PCA --------------------------------------------------
out <- here("data","SCEs","merge",paste0(celltype,".SCT.obj.list.qs"))
s <- here("data","SCEs","merge",paste0(celltype,"preprocessed.merged.SEU.qs"))

if(!file.exists(out)) {
  seu <- qread(here("data","SCEs","preprocessed",paste0(celltype,".preprocessed.SEU.qs")),nthreads=16)
  DefaultAssay(seu) <- "RNA"
  seu <- DietSeurat(seu, assays = "RNA", dimreducs = NULL, graphs = NULL)

  # split object by experiment
  obj.list <- SplitObject(seu, split.by="batchID")

  # normalization
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- SCTransform(x, vst.flavor = "v2", verbose=TRUE)
  })

  # save object list
  qsave(obj.list, file = out, nthreads=16)

  # find most variable features across samples to integrate
  var.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures=3000)

  # merge normalized samples
  seu <- merge(obj.list[[1]], y=obj.list[2:length(obj.list)], merge.data=TRUE)
  DefaultAssay(seu) <- "SCT"

  # manually set variable features of merged Seurat object
  VariableFeatures(seu) <- var.features

  # PCA
  seu <- RunPCA(seu, verbose=TRUE)

  # save merged object
  qsave(seu, file=s, nthreads=16)

  obj.list <- NULL

} else {if (!file.exists(s)) {
  # read object list
  obj.list <- qread(out, nthreads=16)

  # find most variable features across samples to integrate
  var.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures=3000)

  # merge normalized samples
  seu <- merge(obj.list[[1]], y=obj.list[2:length(obj.list)], merge.data=TRUE)
  DefaultAssay(seu) <- "SCT"

  # manually set variable features of merged Seurat object
  VariableFeatures(seu) <- var.features

  # PCA
  seu <- RunPCA(seu, verbose=TRUE)

  qsave(seu, file=s, nthreads=16)

} else {
  seu <- qread(s, nthreads=16)
}}

# Harmony integration ----------------------------------------------------------

# run Harmony
print("Running Harmony integration...")
DefaultAssay(seu) <- "SCT"
set.seed(1990)
seu <- RunHarmony(seu,
                  group.by.vars=c("batchID"),
                  reduction="pca",
                  assay.use="SCT",
                  reduction.save="harmony",
                  plot_convergence=TRUE,
                  .options = harmony_options(max.iter.cluster=50),
                  #early_stop=FALSE,
                  max_iter=100)

# save integrated object
out <- here("data","SCEs","integration",paste0(celltype,".preprocessed.merged.integrated.SEU.qs"))
if(!file.exists(out)) {
  qsave(seu, file = out, nthreads=16)
}

# Leiden clustering ------------------------------------------------------------
if(!dir.exists(here("data","SCEs","clustering"))) {
  dir.create(here("data","SCEs","clustering"), recursive = TRUE)
}
out <- here("data","SCEs","clustering",paste0(celltype,".preprocessed.merged.integrated.clustered.SEU.qs"))

if(!file.exists(out)) {
  DefaultAssay(seu) <- "SCT"
  seu <- FindNeighbors(seu, reduction="harmony", dims=1:50)
  
  plan("multicore", workers=20)
  print("Running Leiden clustering...")
  seu <- FindClusters(seu, algorithm = 4,
                      method = "igraph",
                      resolution = seq(0.2,3,0.2))
  
  # run UMAP
  seu <- RunUMAP(seu, assay="SCT", reduction = "harmony", dims=1:50)

  # # PrepSCTFindMarker
  # plan("multicore", workers=17)
  # seu <- PrepSCTFindMarkers(seu)

  # save clustered object
  qsave(seu, file = out, nthreads=16)
  
  # run clustree
  ggsave(plot=clustree(seu, prefix="SCT_snn_res."), device="pdf",
         width = 10,
         height = 30,
         file=here("data","plots",paste0(celltype,".clustree.",date,".pdf")))
}

# Annotation - Wong
# 3 Reference-based annotation using Wong 2025 dataset
## 3.1 Prepare Wong reference data
out <- "/group/canc2/anson/reference/annotation/Wong_2025/wong_2025.processed.SEU.qs"
seu.ref <- qread(out,nthreads=24)
DefaultAssay(seu.ref) <- "SCT"

## 3.2 Find transfer anchors
# load query object
seu.query <- seu
DefaultAssay(seu.query) <- "SCT"

# Find anchors
if (!file.exists("/group/canc2/anson/reference/annotation/Wong_2025/anchors_Wong_250429.qs")) {
  anchors <- FindTransferAnchors(reference=seu.ref, 
                                 query=seu.query,
                                 normalization.method="SCT",
                                 reference.reduction="pca",
                                 dims=1:50,
                                 query.assay="SCT",
                                 reference.assay="SCT")
  qsave(anchors, file="/group/canc2/anson/reference/annotation/Wong_2025/anchors_Wong_250429.qs", nthreads=16)
} else {
  anchors <- qread("/group/canc2/anson/reference/annotation/Wong_2025/anchors_Wong_250429.qs", nthreads=16)
}

# MapQuery
seu.query <- MapQuery(anchorset = anchors, 
                      reference = seu.ref, 
                      query = seu.query,
                      refdata=list(celltype="celltype",
                                   sampleID="sampleID"), 
                      reference.reduction = "pca", 
                      reduction.model = "umap")

### Save object
out <- here("data","SCEs","annotation","bal.preprocessed.merged.integrated.clustered.Wong.SEU.qs")
if (!file.exists(out)) {qsave(seu.query, file=out, nthreads=16)}

print("Wong annotation done. Please proceed to HLCA azimuth annotation.")
