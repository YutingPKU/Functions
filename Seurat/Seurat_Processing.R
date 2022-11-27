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
#' @param normal A boolean value indicating whther to normalize
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
#' @param ndims Which dimensions to use as input features, used only if features is NULL
#' @param res Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
#'
#' @return seurat obj with normalized, clustering and dimension reduction
#' 

Seurat_PreProcess <- function(combined, normal = F, vars.to.regress, ndims, res){
  if(normal){
    combined <- NormalizeData(combined)
  }
  combined %<>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
    ScaleData(vars.to.regress = vars.to.regress) %>%
    RunPCA(features = VariableFeatures(object = .)) %>%
    FindNeighbors(dims = 1:ndims) %>%
    FindClusters(resolution = res) %>%
    RunUMAP(dims = 1:ndims) 
  
  return(combined)
}



#' Add supplementary information for FindAllMarkers results. 
#' Add exp.pct and exp.avg per cluster, exp.pct/exp.avg difference in target and background cluster
#'
#' @param obj A seurat obj intended to be added meta
#' @param mk Results from FindAllMarkers
#' @param celltype Vectors for all clusters intended to compare
#' 
#' @return dataframe with supplementary information
#' 

AddMeta_Mk <- function(obj, mk, celltype){
  
  obj <- ScaleData(obj, features = mk$gene)
  
  # per cluster
  df.ls <- as.character(unique(mk$cluster)) %>% map(~{
    mg <- mk$gene[which(mk$cluster == .x)]
    gp <- DotPlot(obj, features = mg)
    df <- gp$data
    
    # calculate pct and avg.exp per cluster
    freq.df <- dcast(df, features.plot ~ id, value.var='pct.exp')
    freq.df[,-1] <- freq.df[,-1]/100
    exp.df <- dcast(df, features.plot ~ id, value.var='avg.exp')
    
    # set the target and max-non-target value
    bg.cluster <- setdiff(celltype, .x)
    freq.df$target <- freq.df[,.x]
    freq.df$bg <- apply(freq.df[, bg.cluster], 1, max)
    freq.df$diff <- freq.df$target - freq.df$bg
    
    exp.df$target <- exp.df[,.x]
    exp.df$bg <- apply(exp.df[, bg.cluster], 1, max)
    exp.df$diff <- exp.df$target - exp.df$bg
    
    # merge freq and exp
    colnames(freq.df)[-1] <- paste0('pct.', colnames(freq.df)[-1])
    colnames(exp.df)[-1] <- paste0('exp.', colnames(exp.df)[-1])
    
    
    df <- cbind(freq.df, exp.df[, -1])
    df <- df[match(mg, df$features.plot), ]
    
    return(df)
    
  })
  
  df.all <- do.call('rbind', df.ls)
  
  mk <- cbind(mk, df.all)
  
  return(mk)
}



