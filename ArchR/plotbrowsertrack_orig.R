##################################################
## Project: functions for ArchR plotBrowserTrack
## Script purpose: modify gene tracks and loop tracks
## Date: 2022-11-15
## Author: Yuting Liu
##################################################
library(ArchR)
library(tidyverse)
library(parallel)
library(ggplot2)


# modify gene tracks, only show givn genes
plotGeneTrack_givenGene <- function(ArchRProj = NULL, geneSymbol = NULL, upstream=20000, downstream=20000) {
  
  humanFormat <- function(n) {
    stopifnot(is.numeric(n), n > 0)
    
    size <- log10(n)
    
    if(size < 3) {
      as.character(n)
    } else if(size < 6) {
      paste0(round(n / 1e3), "K")
    } else if(size < 9) {
      paste0(round(n / 1e6), "M")
    } else {
      paste0(round(n / 1e9), "G")
    }
  }
  
  
  gene <- geneSymbol
  stopifnot(gene %in% ArchRProj@geneAnnotation$genes$symbol)
  
  # get regions
  # gene body region
  gr <- ArchRProj@geneAnnotation$genes
  gr.sub <- gr[match(gene, gr$symbol), ]
  
  # plot region
  region <- GenomicRanges::resize(gr.sub, 1, "start")
  region <- extendGR(region, upstream = upstream, downstream = downstream)
  
  
  # extract gene position
  chr <- seqnames(gr.sub) %>% as.vector()
  startPos <- start(gr.sub) %>% as.vector()
  endPos <- end(gr.sub) %>% as.vector()
  strand <- strand(gr.sub) %>% as.vector()
  
  # extract plot psotion
  start <- start(region) %>% as.vector()
  end <- end(region) %>% as.vector()
  
  # extract segments data
  segmentData <- data.table(
    start = seq(startPos, endPos, length.out = 10)[2:9]
  )
  
  if(strand == "-") {
    segmentData[, end := start - 1][]
  }else{
    segmentData[, end := start + 1][]
  }
  
  sizeBar <- log10(endPos - startPos) %>% floor() %>% 10^.
  
  
  # extract exon data
  exonMeta <- ArchRProj@geneAnnotation$exons %>% as.data.table()
  exonData <- exonMeta[symbol == gene]
  
  
  g2 <- ggplot(exonData) +
    geom_linerange(x = 0, ymin = startPos, ymax = endPos, size = 0.5) +
    geom_linerange(aes(x = 0, ymin = start, ymax = end), size = 3) +
    geom_segment(
      data = segmentData, aes(y = start, yend = end), x = 0, xend = 0,
      arrow = arrow(length = unit(.4, "lines")), size = 0.3) +
    # annotate("segment", x = 1, xend = 1, y = end - sizeBar, yend = end, size = 1.5) +
    #  annotate("text", x = 1.2, y = end, label = humanFormat(sizeBar), hjust = 1, size = 3) +
    annotate('text', x = 0, y = end-sizeBar , label = gene, size = 2) +
    scale_y_continuous(expand = c(0, 0), limits = c(start, end)) +
    scale_x_continuous(expand = expansion(c(.1, .1))) +
    labs(x = "", y = "") +
    coord_flip() +
    theme(
      aspect.ratio = 1/10,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
      #plot.margin = margin()
    ) + theme(plot.margin = unit(c(0, 0,0,0), "cm")) 
  
  
  baseSize = 7
  borderWidth = 0.4
  facetbaseSize = 7
  g3 <- g2+ theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, 
                        baseRectSize = borderWidth) + 
    theme(axis.title.x = element_blank(),  
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) + 
    theme(legend.text = element_text(size = baseSize), strip.text.y = element_text(size = facetbaseSize, 
                                                                                   angle = 0)) + 
    guides(fill = guide_legend(override.aes = list(colour = NA, 
                                                   shape = "c", size = 3)), color = "none") + theme(legend.position = "bottom") + 
    theme(legend.title = element_text(size = 5), legend.text = element_text(size = 7), 
          legend.key.size = unit(0.75, "line"), legend.background = element_rect(color = NA), 
          strip.background = element_blank())
  
  
  
  return(g3)
}



# this modified function return plot lists, so it's more esaier to replace any tracks like gene track 
#' @param ylim.bar a numveric vectors of length two, specificying the bar height range in absolute value form, such as c(0,2.3)
#' @param ylim a numveric vectors of length two, specificying the signal range in quantile form, such as c(0.01, 0.99)


#' @return plotlists object

plotBrowserTrack_returnPlotList <-  function(ArchRProj = NULL, region = NULL, groupBy = "Clusters", 
                                             useGroups = NULL, plotSummary = c("bulkTrack", "featureTrack", 
                                                                               "loopTrack", "geneTrack"), 
                                             sizes = c(10, 1.5, 3, 4), 
                                             features = getPeakSet(ArchRProj), loops = getCoAccessibility(ArchRProj), 
                                             geneSymbol = NULL, useMatrix = NULL, log2Norm = TRUE, upstream = 50000, 
                                             downstream = 50000, tileSize = 250, minCells = 25, normMethod = "ReadsInTSS", 
                                             threads = getArchRThreads(), ylim = NULL, ylim.bar = NULL, pal = NULL, baseSize = 7, 
                                             scTileSize = 0.5, scCellsMax = 100, borderWidth = 0.4, tickWidth = 0.4, 
                                             facetbaseSize = 7, geneAnnotation = getGeneAnnotation(ArchRProj), 
                                             title = "", verbose = TRUE, logFile = createLogFile("plotBrowserTrack")) 
{
  
  geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  
  if (is.null(region)) {
    if (!is.null(geneSymbol)) {
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% 
                               tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), 
                                   tolower(geneSymbol)))]
      print(region)
      region <- GenomicRanges::resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  
  
  if (is.null(geneSymbol)) {
    useMatrix <- NULL
  }
  if (!is.null(useMatrix)) {
    featureMat <- ArchR:::.getMatrixValues(ArchRProj = ArchRProj, 
                                           matrixName = useMatrix, name = mcols(region)$symbol)
    if (log2Norm) {
      featureMat <- log2(featureMat + 1)
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, 
                                       drop = FALSE)[rownames(featureMat), 1]
  }
  ggList <- lapply(seq_along(region), function(x) {
    plotList <- list()
    if ("bulktrack" %in% tolower(plotSummary)) {
      
      alter <- function(){
        plotList$bulktrack <- ArchR:::.bulkTracks(ArchRProj = ArchRProj, 
                                                  region = region[x], tileSize = tileSize, groupBy = groupBy, 
                                                  threads = threads, minCells = minCells, pal = pal, 
                                                  ylim = ylim, baseSize = baseSize, borderWidth = borderWidth, 
                                                  tickWidth = tickWidth, facetbaseSize = facetbaseSize, 
                                                  normMethod = normMethod, geneAnnotation = geneAnnotation, 
                                                  title = title, useGroups = useGroups, tstart = tstart, 
                                                  logFile = logFile) + theme(plot.margin = unit(c(0.35, 
                                                                                                  0.75, 0.35, 0.75), "cm"))
      }
      plotList$bulktrack <- bulkTracks_givenBarRange(ArchRProj = ArchRProj, 
                                           region = region[x], tileSize = tileSize, groupBy = groupBy, 
                                           threads = threads, minCells = minCells, pal = pal, 
                                           ylim = ylim, ylim.bar = ylim.bar, baseSize = baseSize, borderWidth = borderWidth, 
                                           tickWidth = tickWidth, facetbaseSize = facetbaseSize, 
                                           normMethod = normMethod, geneAnnotation = geneAnnotation, 
                                           title = title, useGroups = useGroups, tstart = tstart, 
                                           logFile = logFile) + theme(plot.margin = unit(c(0.35, 
                                                                                           0.75, 0.35, 0.75), "cm"))
    }
    if ("sctrack" %in% tolower(plotSummary)) {
      
      plotList$sctrack <- ArchR:::.scTracks(ArchRProj = ArchRProj, 
                                            region = region[x], tileSize = tileSize, groupBy = groupBy, 
                                            threads = threads, minCells = 5, maxCells = scCellsMax, 
                                            pal = pal, baseSize = baseSize, borderWidth = borderWidth, 
                                            tickWidth = tickWidth, scTileSize = scTileSize, 
                                            facetbaseSize = facetbaseSize, geneAnnotation = geneAnnotation, 
                                            title = title, useGroups = useGroups, tstart = tstart, 
                                            logFile = logFile) + theme(plot.margin = unit(c(0.35, 
                                                                                            0.75, 0.35, 0.75), "cm"))
    }
    if ("featuretrack" %in% tolower(plotSummary)) {
      if (!is.null(features)) {
        
        plotList$featuretrack <- ArchR:::.featureTracks(features = features, 
                                                        region = region[x], facetbaseSize = facetbaseSize, 
                                                        hideX = TRUE, title = "Peaks", logFile = logFile) + 
          theme(plot.margin = unit(c(0.1, 0.75, 0.1, 
                                     0.75), "cm"))
      }
    }
    if ("looptrack" %in% tolower(plotSummary)) {
      if (!is.null(loops)) {
        
        plotList$looptrack <- ArchR:::.loopTracks(loops = loops, 
                                                  region = region[x], facetbaseSize = facetbaseSize, 
                                                  hideX = TRUE, hideY = TRUE, title = "Loops", 
                                                  logFile = logFile) + theme(plot.margin = unit(c(0.1, 
                                                                                                  0.75, 0.1, 0.75), "cm"))
      }
    }
    if ("genetrack" %in% tolower(plotSummary)) {
      
      plotList$genetrack <- ArchR:::.geneTracks(geneAnnotation = geneAnnotation, 
                                                region = region[x], facetbaseSize = facetbaseSize, 
                                                title = "Genes", logFile = logFile) + theme(plot.margin = unit(c(0.1, 
                                                                                                                 0.75, 0.1, 0.75), "cm"))
    }
    plotSummary <- tolower(plotSummary)
    names(sizes) <- plotSummary
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]
    sizes <- sizes[tolower(names(plotList))]
    if (!is.null(useMatrix)) {
      suppressWarnings(ArchR:::.combinedFeaturePlot(plotList = plotList, 
                                                    log2Norm = log2Norm, featureMat = featureMat, 
                                                    feature = region[x]$symbol[[1]], useMatrix = useMatrix, 
                                                    pal = pal, sizes = sizes, baseSize = baseSize, 
                                                    facetbaseSize = facetbaseSize, borderWidth = borderWidth, 
                                                    tickWidth = tickWidth))
    }
    else {
      
      return(plotList)
    }
  })
  if (!is.null(mcols(region)$symbol)) {
    names(ggList) <- mcols(region)$symbol
  }
  else {
    if (length(ggList) == 1) {
      ggList <- ggList[[1]]
    }
  }
  ggList
}



# this modified function will plot the bulk signal tracks given y range
#' @param ylim.bar a numveric vectors of length two, specificying the bar height range in absolute value form, such as c(0,2.3)
#' @param ylim a numveric vectors of length two, specificying the signal range in quantile form, such as c(0.01, 0.99)


bulkTracks_givenBarRange <- function (ArchRProj = NULL, region = NULL, tileSize = 100, minCells = 25, 
                                      groupBy = "Clusters", useGroups = NULL, normMethod = "ReadsInTSS", 
                                      threads = 1, ylim = NULL, ylim.bar = NULL, baseSize = 7, borderWidth = 0.4, 
                                      tickWidth = 0.4, facetbaseSize = 7, geneAnnotation = getGeneAnnotation(ArchRProj), 
                                      title = "", pal = NULL, tstart = NULL, verbose = FALSE, logFile = NULL) 
{

  if (is.null(tstart)) {
    tstart <- Sys.time()
  }
  df <- ArchR:::.groupRegionSumArrows(ArchRProj = ArchRProj, groupBy = groupBy, 
                                      normMethod = normMethod, useGroups = useGroups, minCells = minCells, 
                                      region = region, tileSize = tileSize, threads = threads, 
                                      verbose = verbose, logFile = logFile)
  ArchR:::.logThis(split(df, df[, 3]), ".bulkTracks df", logFile = logFile)
  if (!is.null(ylim) ) {
    ylim <- quantile(df$y, ylim)
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  } 
  else {
    ylim <- c(0, quantile(df$y, probs = c(0.999)))
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }
  
  if(is.null(ylim.bar)){
    ylim.bar <- ylim
  }
  
  uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  if (!is.null(useGroups)) {
    uniqueGroups <- unique(useGroups)
  }
  df$group <- factor(df$group, levels = uniqueGroups)
  title <- paste0(as.character(seqnames(region)), ":", start(region) - 
                    1, "-", end(region), " ", title)
  allGroups <- gtools::mixedsort(unique(getCellColData(ArchRProj = ArchRProj, 
                                                       select = groupBy, drop = TRUE)))
  if (is.null(pal)) {
    pal <- suppressWarnings(paletteDiscrete(values = allGroups))
  }
  p <-
    ggplot(df, aes_string("x", "y", color = "group", fill = "group")) +
    geom_area(stat = "identity") + facet_wrap(facets = ~ group,
                                              strip.position = "right",
                                              ncol = 1) + ylab(
                                                sprintf(
                                                  "Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)",
                                                  round(min(ylim.bar), 2),
                                                  round(max(ylim.bar), 2),
                                                  normMethod
                                                )
                                              ) +
    scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(start(region), end(region)),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = ylim.bar,
                       expand = c(0, 0)) + theme_ArchR(
                         baseSize = baseSize,
                         baseRectSize = borderWidth,
                         baseLineSize = tickWidth,
                         legendPosition = "right",
                         axisTickCm = 0.1
                       ) + theme(
                         panel.spacing = unit(0,
                                              "lines"),
                         axis.title.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(),
                         strip.text = element_text(
                           size = facetbaseSize,
                           color = "black",
                           margin = margin(0, 0.35, 0, 0.35,
                                           "cm")
                         ),
                         strip.text.y = element_text(angle = 0),
                         strip.background = element_rect(color = "black")
                       ) + guides(fill = "none",
                                  colour = "none") + ggtitle(title)
  p
}


