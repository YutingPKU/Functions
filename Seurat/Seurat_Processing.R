##################################################
## Project: functions for seurat obj processing
## Script purpose: functions for seurat obj processing
## Date: 2022-09-21
## Author: Yuting Liu
##################################################
library(Seurat)
library(tidyverse)


#' Pre-process for seurat obj by standard log2 method
#'
#' @param combined A seurat obj intended to be pre-processed
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
#' @param ndims Which dimensions to use as input features, used only if features is NULL
#' @param res Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
#'
#' @return seurat obj with normalized, clustering and dimension reduction
#' 

Seurat_PreProcess <- function(combined, vars.to.regress, ndims, res){
  combined %<>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
    ScaleData(vars.to.regress = vars.to.regress) %>%
    RunPCA(features = VariableFeatures(object = .)) %>%
    FindNeighbors(dims = 1:ndims) %>%
    FindClusters(resolution = res) %>%
    RunUMAP(dims = 1:ndims) 
  
  return(combined)
}
