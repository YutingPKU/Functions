##################################################
## Project: functions for ArchR obj finding peak-2-gene links and visualization
## Script purpose: functions for ArchR obj finding peak-2-gene links and visualization
## Date: 2022-11-15
## Author: Yuting Liu
##################################################
library(ArchR)
library(tidyverse)
library(parallel)


#' get peak2gene links given genes of interesting, given as geneSymbol form
#' @param geneSymbol character specificying the gene of interesing
#' @return peak2gene links in GRanges

getPeak2GeneLinks_givenGene  <- function (ArchRProj = NULL, geneSymbol = NULL,
                                          corCutOff = 0.45, FDRCutOff = 1e-04, 
                                          varCutOffATAC = 0.25, varCutOffRNA = 0.25, resolution = 1, 
                                          returnLoops = TRUE) 
{
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  ArchR:::.validInput(input = FDRCutOff, name = "FDRCutOff", valid = "numeric")
  ArchR:::.validInput(input = varCutOffATAC, name = "varCutOffATAC", 
              valid = "numeric")
  ArchR:::.validInput(input = varCutOffRNA, name = "varCutOffRNA", 
              valid = "numeric")
  ArchR:::.validInput(input = resolution, name = "resolution", valid = c("integer", 
                                                                 "null"))
  ArchR:::.validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  if (is.null(ArchRProj@peakSet)) {
    return(NULL)
  }
  if (is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)) {
    return(NULL)
  }
  else {
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= 
                       FDRCutOff), , drop = FALSE]
    if (!is.null(varCutOffATAC)) {
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC), ]
    }
    if (!is.null(varCutOffRNA)) {
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA), ]
    }
    
    if(!is.null(geneSymbol)){
      rna.gr <- metadata(p2g)$geneSet
      loci <- grep(paste(paste0('^', geneSymbol, '$'), collapse = "|"), rna.gr$name)
      rna.gr.ofi <- rna.gr[loci,]
      p2g <- p2g[which(p2g$idxRNA %in% loci), ]
    }
    
    if(nrow(p2g)  ==0){
      return(NULL)
    }
    if (returnLoops) {
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      if (!is.null(resolution)) {
        summitTiles <- floor(start(peakSummits)/resolution) * 
          resolution + floor(resolution/2)
        geneTiles <- floor(start(geneStarts)/resolution) * 
          resolution + floor(resolution/2)
      }
      else {
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneStarts)
      }
      loops <- ArchR:::.constructGR(seqnames = seqnames(peakSummits[p2g$idxATAC]), 
                            start = summitTiles[p2g$idxATAC], end = geneTiles[p2g$idxRNA])
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      loops <- loops[order(mcols(loops)$value, decreasing = TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      loops <- SimpleList(Peak2GeneLinks = loops)
      return(loops)
    }
    else {
      return(p2g)
    }
  }
}


#' get peak2gene links given peaks of interesting, peak id in the form of seATAC idxATAC
#'
#' @param peakid numeric vector given the idxATAC for peaks of interesting
#' @return peak2gene links in GRanges
#' 
getPeak2GeneLinks_subsetPeak <- function (ArchRProj = NULL, peakid = NULL, corCutOff = 0.45, FDRCutOff = 1e-04, 
                                          varCutOffATAC = 0.25, varCutOffRNA = 0.25, resolution = 1, 
                                          returnLoops = TRUE) 
{
  
  if (is.null(ArchRProj@peakSet)) {
    return(NULL)
  }
  if (is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)) {
    return(NULL)
  }
  else {
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= 
                       FDRCutOff), , drop = FALSE]
    if (!is.null(varCutOffATAC)) {
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC), ]
    }
    if (!is.null(varCutOffRNA)) {
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA), ]
    }
    # subset peaks of interesting
    if(!is.null(peakid)) {
      p2g <- p2g[which(p2g$idxATAC %in% peakid), ]
    }
    
    if (returnLoops) {
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      if (!is.null(resolution)) {
        summitTiles <- floor(start(peakSummits)/resolution) * 
          resolution + floor(resolution/2)
        geneTiles <- floor(start(geneStarts)/resolution) * 
          resolution + floor(resolution/2)
      }
      else {
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneStarts)
      }
      loops <- ArchR:::.constructGR(seqnames = seqnames(peakSummits[p2g$idxATAC]), 
                                    start = summitTiles[p2g$idxATAC], end = geneTiles[p2g$idxRNA])
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      loops <- loops[order(mcols(loops)$value, decreasing = TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      loops <- SimpleList(Peak2GeneLinks = loops)
      return(loops)
    }
    else {
      return(p2g)
    }
  }
}


# this change is trying to plot peak2gene heatmap given peaks of interesting and given cell type order
#' @param peakid numeric vector given the idxATAC for peaks of interesting
#' @param celltype.level character vector given the level for groupBy values
#' @return complexheatmap showing the ATAC and correspoding RNA matrix for peak2gene links
#' 

plotPeak2GeneHeatmap_subsetPeak <- function (ArchRProj = NULL, peakid = NULL, celltype.level = NULL,
                                             corCutOff = 0.45, FDRCutOff = 1e-04, 
                                             varCutOffATAC = 0.25, varCutOffRNA = 0.25, k = 25, nPlot = 25000, 
                                             limitsATAC = c(-2, 2), limitsRNA = c(-2, 2), groupBy = "Clusters",
                                             palGroup = NULL, palATAC = paletteContinuous("solarExtra"), 
                                             palRNA = paletteContinuous("blueYellow"), verbose = TRUE, 
                                             returnMatrices = FALSE, seed = 1, logFile = createLogFile("plotPeak2GeneHeatmap")) 
{
  
  #.logThis(mget(names(formals()), sys.frame(sys.nframe())),  "peak2GeneHeatmap Input-Parameters", logFile = logFile)
  if (is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)) {
    stop("No Peak2GeneLinks Found! Try addPeak2GeneLinks!")
  }
  
  ## Section: get peak-to-gene link
  ##################################################
  
  ccd <- getCellColData(ArchRProj, select = groupBy)
  p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
  p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= 
                     FDRCutOff), , drop = FALSE]
  if (!is.null(varCutOffATAC)) {
    p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC), ]
  }
  if (!is.null(varCutOffRNA)) {
    p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA), ]
  }
  
  # subset peaks of interesting
  if(!is.null(peakid)) {
    p2g <- p2g[which(p2g$idxATAC %in% peakid), ]
  }
  
  if (nrow(p2g) == 0) {
    stop("No peak2genelinks found with your cutoffs!")
  }
  if (!file.exists(metadata(p2g)$seATAC)) {
    stop("seATAC does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  if (!file.exists(metadata(p2g)$seRNA)) {
    stop("seRNA does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  
  ## Section: get atac and rna matrix
  ##################################################
  mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
  mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
  p2g$peak <- paste0(rowRanges(mATAC))
  p2g$gene <- rowData(mRNA)$name
  gc()
  mATAC <- assay(mATAC)
  mRNA <- assay(mRNA)
  
  
  ## Section: get kn groups
  ##################################################
  #.logDiffTime(main = "Determining KNN Groups!", t1 = tstart, verbose = verbose, logFile = logFile)
  KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, 
                "list")
  KNNGroups <- lapply(seq_along(KNNList), function(x) {
    KNNx <- KNNList[[x]]
    names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
  }) %>% unlist
  cD <- DataFrame(row.names = paste0("K_", seq_len(ncol(mATAC))), 
                  groupBy = KNNGroups)
  pal <- paletteDiscrete(values = gtools::mixedsort(unique(ccd[, 
                                                               1])))
  if (!is.null(palGroup)) {
    pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% 
                                                                        names(pal)]
  }
  colorMap <- list(groupBy = pal)
  attr(colorMap[[1]], "discrete") <- TRUE
  mATAC <- ArchR:::.rowZscores(mATAC)
  mRNA <- ArchR:::.rowZscores(mRNA)
  rownames(mATAC) <- NULL
  rownames(mRNA) <- NULL
  colnames(mATAC) <- paste0("K_", seq_len(ncol(mATAC)))
  colnames(mRNA) <- paste0("K_", seq_len(ncol(mRNA)))
  rownames(mATAC) <- paste0("P2G_", seq_len(nrow(mATAC)))
  rownames(mRNA) <- paste0("P2G_", seq_len(nrow(mRNA)))
  rownames(p2g) <- paste0("P2G_", seq_len(nrow(p2g)))
  
  ## Section: cluster p2g 
  ##################################################
  #.logDiffTime(main = "Ordering Peak2Gene Links!", t1 = tstart,  verbose = verbose, logFile = logFile)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  k1 <- kmeans(mATAC, k)
  if (nrow(mATAC) > nPlot) {
    nPK <- nPlot * table(k1$cluster)/length(k1$cluster)
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x) {
      idx <- sample(splitK[[x]], floor(nPK[x]))
      k <- rep(x, length(idx))
      DataFrame(k = k, idx = idx)
    }) %>% Reduce("rbind", .)
  }
  else {
    kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
  }
  bS <- ArchR:::.binarySort(t(ArchR:::.groupMeans(t(mATAC[kDF[, 2], ]), kDF[, 1])), clusterCols = TRUE, cutOff = 1)
  rowOrder <- rownames(bS[[1]])
  colOrder <- colnames(bS[[1]])
  
  # order column by given celltype
  if(!is.null(celltype.level)){
    tmp <- cD
    tmp$kcluster = rownames(tmp)
    tmp <- tmp[order(match(tmp$groupBy, celltype.level)),]
    tmp$groupBy <- factor(tmp$groupBy, levels = celltype.level)
    kcluster.order.ls <- split(tmp$kcluster, f = tmp$groupBy)
    kcluster.order.ls %<>% map(~{
      vec <- colOrder[which(colOrder %in% .x)]
      return(vec)
    })
    colOrder.new <- unlist(kcluster.order.ls)
    colOrder <- colOrder.new
  }
  
  
  kDF[, 3] <- as.integer(mapLabels(paste0(kDF[, 1]), newLabels = paste0(seq_along(rowOrder)), 
                                   oldLabels = rowOrder))
  if (returnMatrices) {
    out <- SimpleList(ATAC = SimpleList(matrix = mATAC[kDF[, 
                                                           2], colOrder], kmeansId = kDF[, 3], colData = cD[colOrder, 
                                                                                                            , drop = FALSE]), RNA = SimpleList(matrix = mRNA[kDF[, 
                                                                                                                                                                 2], colOrder], kmeansId = kDF[, 3], colData = cD[colOrder, 
                                                                                                                                                                                                                  , drop = FALSE]), Peak2GeneLinks = p2g[kDF[, 2], 
                                                                                                                                                                                                                  ])
    return(out)
  }
  
  ## Section: heatmap
  ##################################################
  #.logDiffTime(main = "Constructing ATAC Heatmap!", t1 = tstart, verbose = verbose, logFile = logFile)
  htATAC <- ArchR:::.ArchRHeatmap(mat = mATAC[kDF[, 2], colOrder], 
                                  scale = FALSE, limits = limitsATAC, color = palATAC, 
                                  colData = cD[colOrder, , drop = FALSE], colorMap = colorMap, 
                                  clusterCols = FALSE, clusterRows = FALSE, split = kDF[, 
                                                                                        3], labelRows = FALSE, labelCols = FALSE, draw = FALSE, 
                                  name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks"))
  
  
  #.logDiffTime(main = "Constructing RNA Heatmap!", t1 = tstart, verbose = verbose, logFile = logFile)
  htRNA <- ArchR:::.ArchRHeatmap(mat = mRNA[kDF[, 2], colOrder], scale = FALSE, 
                                 limits = limitsRNA, color = palRNA, colData = cD[colOrder, 
                                                                                  , drop = FALSE], colorMap = colorMap, clusterCols = FALSE, 
                                 clusterRows = FALSE, split = kDF[, 3], labelRows = FALSE, 
                                 labelCols = FALSE, draw = FALSE, name = paste0("RNA Z-Scores\n", 
                                                                                nrow(mRNA), " P2GLinks"))
  
  htATAC + htRNA
}



# this modify is trying to find more p-2-g links for cell types with limited cell numbers
# to this end, we try to balance the cell type percentages in seed cells 
#' @param  seed.cell  character vector indicating the seed cells for finding knn clusters, if NULL, random sample knnIteration cells
#' @return add seATAC and seRNA object in the ArchRProj output directory
#' 
addPeak2GeneLinks_givenCellSeed <- function (ArchRProj = NULL, reducedDims = "IterativeLSI", prefix = '', 
                                             useMatrix = "GeneIntegrationMatrix", 
                                             dimsToUse = 1:30, scaleDims = NULL, corCutOff = 0.75, cellsToUse = NULL, 
                                             seed.cell = NULL,
                                             k = 100, knnIteration = 500, overlapCutoff = 0.8, maxDist = 250000, 
                                             scaleTo = 10^4, log2Norm = TRUE, predictionCutoff = 0.4, 
                                             addEmpiricalPval = FALSE, seed = 1, threads = max(floor(getArchRThreads()/2), 
                                                                                               1), verbose = TRUE, logFile = createLogFile("addPeak2GeneLinks")) 
{
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
           "addPeak2GeneLinks Input-Parameters", logFile = logFile)
  .logDiffTime(main = "Getting Available Matrices", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  
  
  AvailableMatrices <- getAvailableMatrices(ArchRProj)
  if ("PeakMatrix" %ni% AvailableMatrices) {
    stop("PeakMatrix not in AvailableMatrices")
  }
  if (useMatrix %ni% AvailableMatrices) {
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }
  ArrowFiles <- getArrowFiles(ArchRProj)
  tstart <- Sys.time()
  dfAll <- .safelapply(seq_along(ArrowFiles), function(x) {
    cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], 
                                                    paste0(useMatrix, "/Info/CellNames")))
    pSx <- tryCatch({
      h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    }, error = function(e) {
      if (getArchRVerbose()) 
        message("No predictionScore found. Continuing without predictionScore!")
      rep(9999999, length(cNx))
    })
    DataFrame(cellNames = cNx, predictionScore = pSx)
  }, threads = threads) %>% Reduce("rbind", .)
  .logDiffTime(sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
                       sum(dfAll[, 2] < predictionCutoff), nrow(dfAll), round(sum(dfAll[, 
                                                                                        2] < predictionCutoff)/nrow(dfAll), 3)), t1 = tstart, 
               verbose = verbose, logFile = logFile)
  keep <- sum(dfAll[, 2] >= predictionCutoff)/nrow(dfAll)
  dfAll <- dfAll[which(dfAll[, 2] > predictionCutoff), ]
  set.seed(seed)
  
  peakSet <- getPeakSet(ArchRProj)
  .logThis(peakSet, "peakSet", logFile = logFile)
  geneSet <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, 
                                                 width = 1), name = geneSet$name, idx = geneSet$idx)
  .logThis(geneStart, "geneStart", logFile = logFile)
  
  
  # get knn clusters
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, 
                       corCutOff = corCutOff, dimsToUse = dimsToUse)
  if (!is.null(cellsToUse)) {
    rD <- rD[cellsToUse, , drop = FALSE]
  }
  # seed cells for generating knn clusters
  if(is.null(seed.cell)){
    idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= 
                    knnIteration)
  }else{
    idx <- seed.cell
  }
  
  
  .logDiffTime(main = "Computing KNN", t1 = tstart, verbose = verbose, 
               logFile = logFile)
  knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx, ], k = k)
  .logDiffTime(main = "Identifying Non-Overlapping KNN pairs", 
               t1 = tstart, verbose = verbose, logFile = logFile)
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * 
                                                         k))
  knnObj <- knnObj[keepKnn == 0, ]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), 
               t1 = tstart, verbose = verbose, logFile = logFile)
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x) {
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)
  
  
  # get pseudobulk rna/atac matrix
  .logDiffTime(main = "Getting Group RNA Matrix", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  groupMatRNA <- ArchR:::.getGroupMatrix(ArrowFiles = getArrowFiles(ArchRProj), 
                                         featureDF = geneDF, groupList = knnObj, useMatrix = useMatrix, 
                                         threads = threads, verbose = FALSE)
  rawMatRNA <- groupMatRNA
  .logThis(groupMatRNA, "groupMatRNA", logFile = logFile)
  .logDiffTime(main = "Getting Group ATAC Matrix", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  groupMatATAC <- ArchR:::.getGroupMatrix(ArrowFiles = getArrowFiles(ArchRProj), 
                                          featureDF = peakDF, groupList = knnObj, useMatrix = "PeakMatrix", 
                                          threads = threads, verbose = FALSE)
  rawMatATAC <- groupMatATAC
  .logThis(groupMatATAC, "groupMatATAC", logFile = logFile)
  
  # normalize pseudobulk rna/atac matrix
  .logDiffTime(main = "Normalizing Group Matrices", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  groupMatRNA <- t(t(groupMatRNA)/colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC)/colSums(groupMatATAC)) * 
    scaleTo
  if (log2Norm) {
    groupMatRNA <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)
  }
  names(geneStart) <- NULL
  
  # generate SE object for rna/atac matrix
  seRNA <- SummarizedExperiment(assays = SimpleList(RNA = groupMatRNA, 
                                                    RawRNA = rawMatRNA), rowRanges = geneStart)
  metadata(seRNA)$KNNList <- knnObj
  .logThis(seRNA, "seRNA", logFile = logFile)
  names(peakSet) <- NULL
  seATAC <- SummarizedExperiment(assays = SimpleList(ATAC = groupMatATAC, 
                                                     RawATAC = rawMatATAC), rowRanges = peakSet)
  metadata(seATAC)$KNNList <- knnObj
  .logThis(seATAC, "seATAC", logFile = logFile)
  rm(groupMatRNA, groupMatATAC)
  gc()
  
  
  # find peak-gene pair
  .logDiffTime(main = "Finding Peak Gene Pairings", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  o <- DataFrame(findOverlaps(.suppressAll(resize(seRNA, 2 * 
                                                    maxDist + 1, "center")), resize(rowRanges(seATAC), 1, 
                                                                                    "center"), ignore.strand = TRUE))
  o$distance <- distance(rowRanges(seRNA)[o[, 1]], rowRanges(seATAC)[o[, 
                                                                       2]])
  colnames(o) <- c("B", "A", "distance")
  if (addEmpiricalPval) {
    .logDiffTime(main = "Computing Background Correlations", 
                 t1 = tstart, verbose = verbose, logFile = logFile)
    nullCor <- ArchR:::.getNullCorrelations(seATAC, seRNA, o, 1000)
  }
  .logDiffTime(main = "Computing Correlations", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), 
                                     assay(seATAC), assay(seRNA))
  o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation/sqrt((pmax(1 - o$Correlation^2, 
                                       1e-17, na.rm = TRUE))/(ncol(seATAC) - 2)))
  o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", 
               "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", 
                     "VarQATAC", "VarQRNA")
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart
  if (addEmpiricalPval) {
    out$EmpPval <- 2 * pnorm(-abs(((out$Correlation - mean(nullCor[[2]]))/sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), 
             showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", 
                      paste0(prefix,  "seATAC-Group-KNN.rds"))
  .safeSaveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", 
                      paste0(prefix, "seRNA-Group-KNN.rds"))
  .safeSaveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA
  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out
  .logDiffTime(main = "Completed Peak2Gene Correlations!", 
               t1 = tstart, verbose = verbose, logFile = logFile)
  .endLogging(logFile = logFile)
  ArchRProj
}


# this function tyring to get peak2gene links based on correlation between given seATAC and seRNA matrix
#' @param prefix character specifying the prefix for seRNA/seATAC SCE object


addPeak2GeneLinks_givenATACRNA <- function(ArchRProj=NULL, maxDist=250000, addEmpiricalPval = FALSE, prefix = ''){
  outATAC <- paste0(getOutputDirectory(ArchRProj), '/Peak2GeneLinks/', prefix, 'seATAC-Group-KNN.rds')
  outRNA <- paste0(getOutputDirectory(ArchRProj), '/Peak2GeneLinks/', prefix, 'seRNA-Group-KNN.rds')
  
  
  seATAC <- readRDS(outATAC)
  seRNA <- readRDS(outRNA)
  
  peakSet <- rowRanges(seATAC)
  geneStart <- rowRanges(seRNA)
  
  o <- DataFrame(findOverlaps(ArchR:::.suppressAll(resize(seRNA, 2 * 
                                                            maxDist + 1, "center")), resize(rowRanges(seATAC), 1, 
                                                                                            "center"), ignore.strand = TRUE))
  o$distance <- distance(rowRanges(seRNA)[o[, 1]], rowRanges(seATAC)[o[, 
                                                                       2]])
  colnames(o) <- c("B", "A", "distance")
  if (addEmpiricalPval) {
    nullCor <- ArchR:::.getNullCorrelations(seATAC, seRNA, o, 1000)
  }
  
  o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), 
                                     assay(seATAC), assay(seRNA))
  o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation/sqrt((pmax(1 - o$Correlation^2, 
                                       1e-17, na.rm = TRUE))/(ncol(seATAC) - 2)))
  o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", 
               "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", 
                     "VarQATAC", "VarQRNA")
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart
  if (addEmpiricalPval) {
    out$EmpPval <- 2 * pnorm(-abs(((out$Correlation - mean(nullCor[[2]]))/sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }
  
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA
  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out
  return(ArchRProj)
  
  
}



