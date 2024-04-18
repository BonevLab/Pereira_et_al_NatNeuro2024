library(scales)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(Gviz)
#library(ggforce)
#library(Sushi)
#library(ggrastr)
#library(paletteer)

#' Optimized discrete color palette generation
#'
#' This function assesses the number of inputs and returns a discrete color palette that is tailored to provide the most
#' possible color contrast from the designated color set.
#'
#' @param values A character vector containing the sample names that will be used. Each entry in this character vector will be
#' given a unique color from the designated palette set.
#' @param set The name of a color palette provided in the `ArchRPalettes` list object.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteDiscrete <- function(
  values = NULL,
  set = "stallion",  
  reverse = FALSE
){
  
  .validInput(input = set, name = "set", valid = c("character"))
  .validInput(input = values, name = "values", valid = c("character", "factor"))
  .validInput(input = reverse, name = "reverse", valid = c("boolean"))
  
  values <- gtools::mixedsort(values)
  n <- length(unique(values))
  pal <- ArchRPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))] #mixed sort gets 1,2,3,4..10,11,12
  
  if(n > length(palOrdered)){
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  }else{
    palOut <- palOrdered[seq_len(n)]
  }
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  names(palOut) <- unique(values)
  
  return(palOut)
  
}

#' Continuous Color Palette
#'
#' @param set The name of a color palette provided in the `ArchRPalettes` list object.
#' @param n The number of unique colors to generate as part of this continuous color palette.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteContinuous <- function(
  set = "solarExtra", 
  n = 256, 
  reverse = FALSE
){
  
  .validInput(input = set, name = "set", valid = c("character"))
  .validInput(input = n, name = "n", valid = c("integer"))
  .validInput(input = reverse, name = "reverse", valid = c("boolean"))
  
  pal <- ArchRPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  return(palOut)
  
}


extractFeatures <- function(object,features='seurat_clusters',cells_toInclude='all',cells_toExclude='none',min.cutoff=NA,max.cutoff=NA){
  plot.data <- FetchData(object,vars=c('UMAP_1','UMAP_2',features))
  plot.data$labels <- Idents(object)
  plot.data$reps <- object$orig.ident
  if(cells_toInclude!='all'){
    plot.data <- plot.data[grepl(paste(cells_toInclude,collapse="|"),plot.data$labels),]
  }
  if(cells_toInclude!='none'){
    plot.data <- plot.data[!grepl(paste(cells_toExclude,collapse="|"),plot.data$labels),]
  }
  colnames(plot.data)[3] <- 'feature'
  #  feature <- plot.data[,3]
  if(!(features=='seurat_clusters'|features=='orig.ident')){
    min.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = min(plot.data[, 'feature']),
          no = cutoff
        ))
      },
      cutoff = min.cutoff,
      feature = 'feature'
    )
    max.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = max(plot.data[, feature]),
          no = cutoff
        ))
      },
      cutoff = max.cutoff,
      feature = 'feature'
    )
    plot.data[, 3] <- as.vector(sapply(
      X = 3,
      FUN = function(index) {
        plot.data.feature <- as.vector(x = plot.data[, index])
        min.use <- SetQuantile(cutoff = min.cutoff[index - 2], plot.data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[index - 2], plot.data.feature)
        plot.data.feature[plot.data.feature < min.use] <- min.use
        plot.data.feature[plot.data.feature > max.use] <- max.use
        return(plot.data.feature)
      }
    ))
  }
  return(plot.data)
}

groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    #  data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

seqplots_heatmap <- function(tracks,peaks,window,bin,ignore_strand=T,idx=NULL,show_n=T,clusters=NULL,which_index=1,cols,zlims=NULL,raster_quality=5,show_leg=F,type = 'mf'){
  peaks_df <- read.table(peaks,header=F)
  res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = type,add_heatmap=T,xmin=window,xmax=window,bin = bin,ignore_strand = ignore_strand)
  hm_names <- names(res$data[[1]])
  hm_l <- list()
  if(!is.null(clusters)){
    if (is.null(idx)){
      idx <- order(clusters,rowSums(res$data[[1]][[which_index]]$heatmap),decreasing=T)   #Reorder all based on signal in first hm
    }
    row_split=clusters[idx]
  } else {
    idx <- order(rowSums(res$data[[1]][[which_index]]$heatmap),decreasing=T)
    row_split <- NULL
  }
  if(!is.null(zlims)){
    if(is.vector(zlims)){
      zlims <- matrix(rep(zlims,length(hm_names)),ncol=2,byrow = T)
    }
  }
  for (i in 1:length(hm_names)) {                 
    res_hm <- res$data[[1]][[i]]$heatmap
    res_hm <- res_hm[idx,]
    if(i==1){
      la <- rowAnnotation(foo = anno_block(labels=ifelse(show_n,paste0('n=',table(clusters)),NULL),gp = gpar(fill = rev(colorpalette('brbg',max(as.numeric(clusters)))))))
    } else {
      la=NULL
    }
    hm <- Heatmap(res_hm,name = hm_names[[i]],cluster_rows = F,cluster_columns=F,show_row_names = F,show_column_names = F,row_title = NULL,left_annotation = la,cluster_row_slices = F,show_row_dend = F,row_split = row_split,
                  col =  colorRamp2(seq(zlims[i,1],zlims[i,2],length.out=length(cols)), cols),use_raster = TRUE, raster_quality = raster_quality,show_heatmap_legend = show_leg,
                  border=T,heatmap_legend_param=list(at=seq(zlims[i,1],zlims[i,2],length.out=3),direction = "horizontal",title='CPM',xjust=0.5,yjust=0.5))
    hm_l[[hm_names[[i]]]] <- hm
  }
  return(list(hm_l=hm_l,idx=idx))
}




plotMisha <- function(main_f,object,targetGene,outDir='figures/',gene_track=readRDS(paste0('/home/hpc/bonev/annotations/',genome,'/',genome,'_gviz.RDS')),out_f=NULL,upstream=1e5,downstream=1e5,plotClusters=levels(object),chipTracksToExtract = gtrack.ls('scATAC.E14'),conditions='',targetEnh=NULL,plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=TRUE, arcs=FALSE, rna=FALSE, chip=TRUE,axis=FALSE,scATAC=TRUE, genes=TRUE,ideogram=TRUE),cluster_cols,sample_cells=NULL,arcIntervals=NULL,...){
  suppressWarnings(suppressMessages(source(paste0(main_f,'config.R'))))
  targetCell = cells[1]
  targetRegion = targetGene
  targetEnh <- NULL
  scoreTrackToExtract=conditions
  domainsToPlot = paste0("hic.",targetCell,".ins_250_domains_expanded")
  rnaTracksToExtract = ''
  maxDist <- max(upstream,downstream)
  chipRes = 10
  zoomRegion <- c(50-upstream/maxDist*50,50+downstream/maxDist*50)
  pointCEX = 0.4
  plotScale=FALSE
  cex.axis=1.2
  chipYlim <- ''
  rnaYlim=c(0)
  window_scale=1.1
  img_factor=1
  #plotOrder <- list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=TRUE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=TRUE, genes=TRUE,ideogram=TRUE)
  plotRatios <- list(unitHeight=120, scores=2.5/window_scale, VP=1.5,UMI4C=2.5, loops=2.2, rna=0.6, chip=0.6,meth=0.6, domains=0.15, genes=0.7, arcs=0.7, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15)
  # if(maxDist>=1e5){plotRatios$genes <- plotRatios$genes*1.5} 
  vTypeScore <- 'max'
  vTypeChip <- 'avg'
  zoomInterval <- zoomRegion
  rnaNames <- c('')
  #  plotClusters <- c('NSC','IPC','PN1','PN2')
  chipNames <- ''
  chipColors <- ''
  methNames <- ''
  methColors <- ''
  rnaColors <- NA
  gene_color <- 'brown'
  figure_mode <- TRUE
  geneNames=TRUE
  radius=25e3
  if (!is.null(sample_cells)){
    scCellsToPlot <- sample(colnames(object),sample_cells)
  } else {
    scCellsToPlot <- NULL
  }
  scATAC_window <- 100
  plot.dgram=FALSE
  binSize <- 2e3
  loopThr <- -101
  arcThr <- 59
  widthCorrFactor <- 15
  arcCol <- c(4,3,2)
  arcType <- c(3,4,2)
  getArcScores <- TRUE
  viewpointOnly <- FALSE
  tssOnly <- FALSE
  exactRegion <- TRUE
  ann2D <- gintervals.ls('domains_expanded')[1]
  annIntervals <- gintervals.ls('domains_expanded')[1]
  leftMargin <- 5
  rightMargin <- 3.35
  suppressWarnings(suppressMessages(source('/home/hpc/bonev/projects/hic/dfg/scripts/temp_functions.R')))
  geneTableDir <- paste0(main_f,'data/extractedGenes/')
  dir.create(geneTableDir,showWarnings = F,recursive = T)
  if(exists('chipTracksToExtract')){
    chipYlim <- matrix(,length(chipTracksToExtract),2)
  }
  #########################
  if (length(unlist(strsplit(targetRegion,split = ',')))==1){
    targetGene=targetRegion
  } else if (length(unlist(strsplit(targetRegion,split = ',')))==3){
    targetCoordinates=targetRegion
  } else if (length(unlist(strsplit(targetRegion,split = ',')))==6){
    squareCoordinates=targetRegion
  } else {
    stop('Unknown interval format. Exiting ...')
  }
  ######################################
  if(exists('targetCoordinates') & exactRegion == TRUE){maxDist <- (as.numeric(as.character(unlist(strsplit(targetCoordinates,','))[3]))-as.numeric(as.character(unlist(strsplit(targetCoordinates,','))[2])))/2}
  
  ###########################
  ##################################################################################################################################################################
  ##################################################################################################################################################################
  if (exists("scoreTrackToExtract")){
    if (exists('targetGene')){
      outName <- paste0(targetGene,'_',maxDist,'_',paste0(scoreTrackToExtract,collapse='and'),'_',paste0(zoomRegion,collapse='-'))
    } else if (exists('targetCoordinates')) {
      outName <- paste0(paste0(targetCoordinates,collapse='_'),'_',maxDist,'_',scoreTrackToExtract,'_',paste0(zoomRegion,collapse='-'))
    } else { 
      outName <- paste0(paste0(squareCoordinates,collapse='_'),'_',maxDist,'_',scoreTrackToExtract) 
    }
  } else {
    outName <- paste0(ifelse(exists("targetGene"),targetGene,targetCoordinates),'_',maxDist)
  }
  #message(scoreTrack)
  if (!is.null(out_f)){outName <- out_f}
  ###########################
  availableTables <- list.files(geneTableDir, full.names=F)
  plotMar <- c(leftMargin,rightMargin)
  ###########################
  ### LOAD INTERVALS
  tssCoordinates <- gintervals.load(tss_f)
  tssCoordinates <- tssCoordinates[!duplicated(tssCoordinates[,5]),]  #Fix duplicated names
  rownames(tssCoordinates) <- tssCoordinates$geneName
  geneCoordinates <- gintervals.load(genes_f)
  geneCoordinates <- geneCoordinates[!duplicated(geneCoordinates[,5]),]  #Fix duplicated names
  rownames(geneCoordinates) <- geneCoordinates$geneName
  
  annIntervals <- gintervals.load(annIntervals)
  if(!is.null(arcIntervals)){
    if(!is.data.frame(arcIntervals)){   #Assume is granges
      arcIntervals <- data.frame(chrom1=seqnames(arcIntervals),start1=as.numeric(start(arcIntervals)),end1=as.numeric(end(arcIntervals )),chrom2=arcIntervals$gene_chr,start2=as.numeric(arcIntervals$gene_start),end2=arcIntervals$gene_start+1,Correlation=arcIntervals$Correlation,gene_name=as.character(arcIntervals$gene_name))
    }
    arcIntervals <- arcIntervals[arcIntervals$gene_name==targetGene,]
    arcColor_by='Correlation'
    arcColors <- colorRampPalette(c('darkblue','blue','grey80','red','darkred'))
  }
  #### Evaluate additional arguments #####
  list2env(list(...), envir = environment())
  ##################
  ann2D <- gintervals.load(ann2D)
  ann2D <- ann2D[1,]
  ann2D$cluster=1
  if (!is.null(targetEnh)){
    annIntervals <- rbind(tssCoordinates[targetGene,1:3],targetEnh)
    annIntervals$cluster <- targetGene
    annIntervals$cluster[2:nrow(annIntervals)] <- 'Enhancer'
    annIntervals <- intervals.expand(annIntervals,expansion = 1000)
  } else {
    #plotOrder$anno <- FALSE
  }
  arguments <- list(...)
  paste(arguments)
  ##################################################################################################################################################################
  #radius <- 5e3
  if (!is.null(arcIntervals)){
    loopDetails <- arcIntervals
    # radius <- round(abs(distance_f$dist)/20,-3)
    # loopDetails <- gintervals.2d(distance_f[,1],distance_f[,2],distance_f[,3],distance_f[,8],distance_f[,9],distance_f[,10])
    # if(abs(distance_f$dist)>=1e6){
    #   pointCEX = 0.1
    # } else if(abs(distance_f$dist)>=1e5){
    #   pointCEX = 0.25 
    # } else {
    #   pointCEX = 0.5 
    # }
  } else {
    loopDetails <- gintervals.2d(1,1,2,1,1,2)
  }
  
  genePairs <- c()
  ##################################################################################################################################################################
  
  ##################################################################################################################################################################
  if(!exists('binSize')){binSize=5e+3}
  ##################################################################################################################################################################
  # next
  
  plotSecondary <- 0
  if (exists('squareCoordinates'))
  {
    plotSecondary <- 1
    conditions <- ''
    arcsToDraw <- ''
    zoomRegion <- '0|100'
    targetGene <- paste0(squareCoordinates,collapse='|')
    
    chr <- as.character(unlist(strsplit(squareCoordinates, ","))[1])
    start1 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[2]))
    stop1 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[3]))
    start2 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[5]))
    stop2 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[6]))
    ##########################################
    fullInterval <- gintervals.2d(chr,start1,stop1,chr,start2,stop2)
    currentCoordinates <- gintervals(chr,start1,stop2)
    currentIntrv <- gintervals(chr,start1,stop2)
  }else if (exists('targetCoordinates'))
  {
    posData <- unlist(strsplit(targetCoordinates, ","))
    mid <- floor(as.numeric(as.character(posData[2]))+(as.numeric(as.character(posData[3]))-as.numeric(as.character(posData[2])))/2)
    currentCoordinates <- gintervals(as.character(posData[1]),mid,mid+1)
    #targetGene <- targetCoordinates
  }else{
    currentCoordinates <- tssCoordinates[which(tssCoordinates$geneName == targetGene),]
    centerTarget=TRUE
    if (nrow(currentCoordinates) == 0){stop(paste0('\n-------------------------\n',targetGene,' NOT FOUND IN DATABASE...\nTRY INTERVAL INPUT... (not available yet, but soon...)'),'\n-------------------------\n');}
  }
  
  ###########################
  # next
  
  plotSecondary <- 0
  
  ###########################
  ### COMPUTE 2D INTERVALS TO EXTRACT
  if(plotSecondary == 0)
  {
    iter2d <- makeIter2D_TSS(currentCoordinates,maxDist,binSize)
    currentIntrv <- gintervals(unique(iter2d$chrom1), min(iter2d$start2), max(iter2d$end2))
    fullInterval <- gintervals.2d(currentIntrv$chrom,min(iter2d$start2),max(iter2d$end2),currentIntrv$chrom,min(iter2d$start2),max(iter2d$end2))
  }
  
  ###########################
  plotLength <- currentIntrv$end-currentIntrv$start
  center <- currentIntrv$start+plotLength/2
  
  ###########################
  
  if(exists('zoomRegion'))
  {
    plotRegionStart <- currentIntrv$start+plotLength*(as.numeric(as.character(zoomInterval[1]))/100)
    plotRegionEnd <- currentIntrv$start+plotLength*(as.numeric(as.character(zoomInterval[2]))/100)
    currentIntrv <- gintervals(currentIntrv$chrom, plotRegionStart, plotRegionEnd)
  }
  
  ############ Extract HiC KNN scores for the full region
  if (exists('scoreTrackToExtract'))
  {
    if (length(scoreTrackToExtract)>=1)
    {
      checkFile <- paste(currentCoordinates$chrom,currentCoordinates$start,currentCoordinates$end,maxDist,paste0(targetGene,collapse='_'),paste0(scoreTrackToExtract,collapse='and'),sep='_')
      if (!(checkFile %in% availableTables))
      { 
        print('EXTRACT TABLE...')
        startTime <- proc.time()
        fullInterval2 <- intervals.2d.expand(fullInterval,(fullInterval$end1-fullInterval$start1)*2,(fullInterval$end2-fullInterval$start2)*2) 
        extractedScores <- list()
        for (scoreTrack in scoreTrackToExtract){
          extractedScores[[scoreTrack]] <- gextract(scoreTrack, fullInterval2,band=c(-2e8,-1000))
        }  
        save(extractedScores, file=paste0(geneTableDir,'/',checkFile))
        stopTime <- proc.time() - startTime
        print(paste(checkFile, 'DONE... in ', stopTime['elapsed'], 'seconds'))
      }else{
        print('USING MEMORY DATA...')			
        extractedScores <- get(load(paste0(geneTableDir,'/',checkFile)))
      }
    } else {plotOrder[['scores']]=FALSE}
  }else{scoreTrackToExtract=''; plotOrder[['scores']]=FALSE}
  
  ############ Extract HiC data for the marginals
  if (exists('conditions'))
  {
    if (nchar(conditions[1]) > 1)
    {
      #   scoreTracks <- gtrack.ls(conditions)
      scoreTracks <- conditions
      vTracks <- paste0("v_",scoreTracks)
      for(set in scoreTracks){gvtrack.create(paste0("v_",set),set, vTypeScore)}
      extReads <- gextract(vTracks, iter2d, iterator=iter2d)
    }else{conditions=''; plotOrder[['VP']]=FALSE; plotOrder[['loops']]=FALSE}
  }else{conditions=''; plotOrder[['VP']]=FALSE; plotOrder[['loops']]=FALSE}
  ########################################################
  ############ Extract HiC data for the UMI4C
  if (exists('umi4c_tracks'))
  {
    if (nchar(umi4c_tracks[1]) > 1)
    {
      library("umi4cPackage")
      p4cLoadConfFiles(conf_dir = umi4c_conf_dir)
      umi4c_list <- list()
      for(set in 1:length(umi4c_tracks)){
        bait_coord <- as.numeric(gtrack.attr.get(umi4c_tracks[set][1], "Bait_coord"))
        umi4c_list[[set]] <- p4cNewProfile(umi4c_tracks[set], scope_5=abs(bait_coord-currentIntrv$start), scope_3=abs(bait_coord-currentIntrv$end))
      }
      if(!exists('umi4c_ylim')){umi4c_ylim=NULL}
      if(!exists('umi4c_main')){umi4c_main=NULL}
    }else{umi4c_tracks=''; plotOrder[['UMI4C']]=FALSE; plotOrder[['loops']]=FALSE}
  }else{umi4c_tracks=''; plotOrder[['UMI4C']]=FALSE; plotOrder[['loops']]=FALSE}
  ########################################################
  
  ############ Extract HiC data for the arcs
  
  ########################################################
  
  ############ Extract ChIPseq data for the full region
  if (exists('chipTracksToExtract'))
  {
    if (nchar(chipTracksToExtract[1]) > 1)
    {
      #tracks <- as.vector(unlist(sapply(chipTracksToExtract,function(x){gtrack.ls(x)},simplify=T)))
      tracks <- chipTracksToExtract
      #tracks <- gtrack.ls(chipTracksToExtract)
      chipTracks <- paste0("v_",tracks)
      for(set in tracks){gvtrack.create(paste0("v_",set),set, vTypeChip)}
      if (is.null(chipRes)){
        chipRes <- 10		
        if(plotLength > 1e+6 & plotLength < 5e+6){chipRes <- 100
        } else if(plotLength > 5e+6){chipRes <- 1000}
      }
      chipData <- gextract(chipTracks, currentIntrv, iterator=chipRes)
      chipTracksToPlot <- unlist(strsplit(chipTracksToExtract, ","))
      if(plotSecondary ==1)
      {
        chipData_UP <- subset(chipData, chipData$start > fullInterval$start1 & chipData$end < fullInterval$end1)
        chipData_DOWN <- subset(chipData, chipData$start > fullInterval$start2 & chipData$end < fullInterval$end2)
      }
      
    }else{chipTracksToExtract=''; plotOrder[['chip']]=FALSE}
  }else{chipTracksToExtract=''; plotOrder[['chip']]=FALSE}
  ########################################################
  ############ Extract MEth data for the full region
  if (exists('methTracksToExtract'))
  {
    #methTracksToExtract <- methTracksToExtract
    if (nchar(methTracksToExtract[1]) > 1)
    {
      #methTracks <- as.vector(unlist(sapply(methTracksToExtract,function(x){gtrack.ls(x)},simplify=T)))
      methTracks <- methTracksToExtract
      methList <- list()
      for (i in 1:length(methTracks)){
        temp_meth <- gextract(methTracks, currentIntrv, iterator=methTracks[i])
        methList[[i]] <- unique(temp_meth) 
      }
      methData <- do.call("rbind", methList)
      methTracksToPlot <- unlist(strsplit(methTracksToExtract, ","))
      if(plotSecondary ==1)
      {
        methData_UP <- subset(methData, methData$start > fullInterval$start1 & methData$end < fullInterval$end1)
        methData_DOWN <- subset(methData, methData$start > fullInterval$start2 & methData$end < fullInterval$end2)
      }
      
    }else{methTracksToExtract=''; plotOrder[['meth']]=FALSE}
  }else{methTracksToExtract=''; plotOrder[['meth']]=FALSE}
  ########################################################
  ############ Extract RNAseq data for the full region
  if (exists('rnaTracksToExtract'))
  {
    if (nchar(rnaTracksToExtract[1]) >1)
    {
      rTracks <- as.vector(unlist(sapply(rnaTracksToExtract,function(x){gtrack.ls(x)},simplify=T)))
      #rTracks <- gtrack.ls(rnaTracksToExtract)
      rnaTracks <- paste0("v_",rTracks)
      for(set in rTracks){gvtrack.create(paste0("v_",set),set, vTypeChip)}
      ####################
      # chipRes <- 10
      if(plotLength > 1e+6 & plotLength < 5e+6){chipRes <- 100
      } else if(plotLength > 5e+6){chipRes <- 1000}
      ####################
      rnaData <- gextract(rnaTracks, currentIntrv, iterator=chipRes)
      rnaData[is.na(rnaData)] <- 0
      
      if(plotSecondary ==1)
      {
        rnaData_UP <- subset(rnaData, rnaData$start > fullInterval$start1 & rnaData$end < fullInterval$end1)
        rnaData_DOWN <- subset(rnaData, rnaData$start > fullInterval$start2 & rnaData$end < fullInterval$end2)
      }
      
    }else{rnaTracksToExtract=''; plotOrder[['rna']]=FALSE}
  }else{rnaTracksToExtract=''; plotOrder[['rna']]=FALSE}
  ########################################################
  # next
  if(!exists('domainsToPlot')){domainsToPlot=''; plotOrder[['domains']]=FALSE}
  
  ##################################################################################################################################################################
  
  ############ PLOT DATA
  # plotGenes <- subset(geneCoordinates, geneCoordinates$chrom == as.character(currentIntrv$chrom) & geneCoordinates$start > currentIntrv$start & geneCoordinates$end < currentIntrv$end)
  plotGenes <- subset(geneCoordinates, geneCoordinates$chrom == as.character(currentIntrv$chrom) & geneCoordinates$end > currentIntrv$start & geneCoordinates$start < currentIntrv$end)
  plotAnn <- subset(annIntervals, annIntervals$chrom == as.character(currentIntrv$chrom) & annIntervals$end > currentIntrv$start & annIntervals$start < currentIntrv$end)
  
  plot2Dann <- subset(ann2D, ann2D$chrom == as.character(currentIntrv$chrom) & ann2D$end > currentIntrv$start & ann2D$start < currentIntrv$end)
  
  hicNames <- rev(paste0('v_',unlist(strsplit(conditions,"\\|"))))
  plotLim <- c(currentIntrv$start,currentIntrv$end)
  
  ############ Set output image parameters
  ######################################################
  scoreSets <- length(conditions)
  chipSets <- length(chipTracksToExtract)
  methSets <- length(methTracksToExtract)
  rnaSets <- length(rnaTracksToExtract)
  domainSets <- length(domainsToPlot)
  umi4cSets <- ifelse(plot.dgram,length(umi4c_tracks)*2,length(umi4c_tracks))
  ######################################################
  if(plotOrder[['scores']] == TRUE){plotOrder[['scores']] <- length(scoreTrackToExtract)}
  if(plotOrder[['UMI4C']] == TRUE){plotOrder[['UMI4C']] <- umi4cSets}
  if(plotOrder[['chip']] == TRUE){
    plotOrder[['chip']] <- chipSets
    atac_ylim <- c(0,ifelse(ncol(chipData[,grep('ATAC',colnames(chipData),ignore.case = F)])==0,2,max(chipData[,grep('ATAC',colnames(chipData),ignore.case = F)],na.rm=T))*1.1)
  }
  if(plotOrder[['meth']] == TRUE){plotOrder[['meth']] <- methSets}
  if(plotOrder[['rna']] == TRUE){
    plotOrder[['rna']] <- rnaSets
    rna_ylim <- c(-(max(rnaData[,grep('_Rev',colnames(rnaData))],na.rm=T)*1.1),max(rnaData[,grep('_For',colnames(rnaData))],na.rm=T)*1.1)
  }
  if(plotOrder[['domains']] == TRUE){plotOrder[['domains']] <- domainSets}
  
  ######################################################
  plotSets <- sum(unlist(plotOrder))
  imgGrid <- matrix(1:plotSets, nrow=plotSets)
  ######################################################
  plotHeights <- c()
  for (set in names(plotOrder))
  {
    currentSet <- plotOrder[[set]]
    currentRatio <- plotRatios[[set]]
    if (currentSet == TRUE | currentSet == 1){plotHeights <- append(plotHeights, currentRatio)
    }else if(currentSet > 1){plotHeights <- append(plotHeights, rep(currentRatio,currentSet))}
  }
  
  ############ START THE PLOT
  ########################################################
  imgHeight <- sum(plotHeights)*plotRatios[['unitHeight']]
  imgWidth <- img_factor*10*plotRatios[['unitHeight']]
  #pngName <- paste0(outDir,outName,'.png')
  #if(plotOrder[['scores']] == FALSE)
  #{
  if (!plotOrder[['VP']]){
    pdfName <- paste0(outDir,outName,'.pdf')
  } else {
    pdfName <- paste0(outDir,outName,'_',targetCell,'.pdf')
  }
  imgHeight <- imgHeight/96
  imgWidth <- imgWidth/96
  pdf(file=pdfName,width=imgWidth, height=imgHeight)
  #}else{png(file=pngName,width=imgWidth,height=imgHeight)}
  #png(file=pngName,width=imgWidth, height=imgHeight)
  
  #########################################################
  grid_layout <- grid.layout(nrow=nrow(imgGrid), heights=plotHeights, respect=F)
  plot.new()
  pushViewport(viewport(layout=grid_layout))
  
  i=1
  for (set in names(plotOrder))
  {
    
    if (plotOrder[[set]] == FALSE){next
    } else {
      # print(set)
      if(set == 'genes') {
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0.1,plotMar[1]-0.3,0.1,plotMar[2])))
      } else if(set == 'ideogram'){
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0,plotMar[1]+1,0,plotMar[2]+1)))
      } else if(set == 'scATAC'){
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0,plotMar[1]-0.3,0,plotMar[2]-1.3)))
      } else if(set != 'scores'){
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1))
      }
      par(fig=gridFIG())
      par(new=TRUE)
      # if(set == 'scores'){
      #   for (scoreSet_counter in 1:length(extractedScores)){
      #     if(scoreSet_counter>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
      #     plotSingleKNN(extractedScores[[scoreSet_counter]], plotAnn = FALSE, plotMar, plotLim, plot2Dann, loopThr, loopDetails, radius,pointCEX=pointCEX,plotScale=plotScale,cex.axis=cex.axis,window_scale=window_scale)
      #     if(length(extractedScores)>1&scoreSet_counter!=length(extractedScores)){par(new=TRUE);upViewport();i=i+1}
      #   }
      if(set == 'scores'){
        for (scoreSet_counter in 1:length(extractedScores)){
          pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0,plotMar[1],0.05,plotMar[2])))
          p <- plotSingleKNN_ggplot(extractedScores[[scoreSet_counter]], plotAnn = FALSE, plotMar, plotLim, plot2Dann,plotThr = loopThr,loopDetails =  loopDetails,radius =  radius,pointCEX=pointCEX,plotScale=plotScale,cex.axis=cex.axis,window_scale=window_scale)
          print(p, vp=viewport(layout.pos.row = i, layout.pos.col = 1))
          if(scoreSet_counter!=length(extractedScores)){
            popViewport();upViewport();i=i+1
          }
        }
      }else if(set == 'UMI4C'){
        if(plot.dgram){
          for (scoreSet_counter in 1:length(umi4c_tracks)){
            if(scoreSet_counter>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}         #FALSE
            plot(umi4c_list[[scoreSet_counter]],min_win_cov=min_umi4c_covs[scoreSet_counter],plot.trend=T,plot.dgram=F,add_xcoord=F,plotMar=c(0,plotMar[1],0.5,plotMar[2]),plotMar2=c(0,plotMar[1],0.5,plotMar[2]),plotOma=c(0,0,0,0),plot.colorbar=plot.colorbar,plot.misha=T,ylim=umi4c_ylim,main=umi4c_main[scoreSet_counter],cex.axis=cex.axis)
            par(new=TRUE);upViewport();i=i+1      #i=2
            pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)
            plot(umi4c_list[[scoreSet_counter]],min_win_cov=min_umi4c_covs[scoreSet_counter],plot.trend=F,plot.dgram=T,add_xcoord=F,plotMar=c(0,plotMar[1],0.5,plotMar[2]),plotMar2=c(0.5,plotMar[1],0,plotMar[2]),plotOma=c(0,0,0,0),plot.colorbar=plot.colorbar,plot.misha=T)
            if(length(umi4c_tracks)>1&scoreSet_counter!=length(umi4c_tracks)){par(new=TRUE);upViewport();i=i+1}          
          }
        }else{
          for (scoreSet_counter in 1:length(umi4c_tracks)){
            if(scoreSet_counter>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
            plot(umi4c_list[[scoreSet_counter]],min_win_cov=min_umi4c_covs[scoreSet_counter],plot.trend=T,plot.dgram=F,add_xcoord=F,plotMar=c(0.5,plotMar[1],0.5,plotMar[2]),plotMar2=c(0.5,plotMar[1],0,plotMar[2]),plotOma=c(0,0,0,0),plot.colorbar=plot.colorbar,plot.misha=T,ylim=umi4c_ylim,main=umi4c_main[scoreSet_counter],cex.axis=cex.axis)
            if(length(umi4c_tracks)>1&scoreSet_counter!=length(umi4c_tracks)){par(new=TRUE);upViewport();i=i+1}
          }
        }
      }else if(set == 'genes'){plot_genes(genome,gene_track=gene_track[[as.character(currentIntrv$chrom)]], currentIntrv, gene_stacking="squish",gene_color=gene_color,gene_size=0.7,fontsize=20,targetGene=plotGenes[plotGenes$geneName%in%targetGene,],collapseTranscripts='meta',geneNames=geneNames)
      }else if(set == 'rna'){
        for (rnaSet in 1:rnaSets){
          if (rnaSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          combineRNA(rnaData, rnaTracksToExtract[rnaSet], plotMar, yLim=rna_ylim,cex.axis=cex.axis,setNames=rnaNames[rnaSet],rnaColors=rnaColors[rnaSet],figure_mode=figure_mode)
          if(rnaSets>1&rnaSet!=rnaSets){par(new=TRUE);upViewport();i=i+1}
        } 
      }else if(set == 'chip'){
        for (chipSet in 1:chipSets){
          if (chipSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          if(grepl('ATAC',chipTracksToPlot[chipSet],ignore.case = F)){
            chipYlim[chipSet,] <- t(as.matrix(atac_ylim))
          } 
          plotChIP(chipData, chipTracksToPlot[chipSet], plotMar, yLim=chipYlim,cex.axis=cex.axis,setNames=chipNames[chipSet],chipColors=chipColors[chipSet],figure_mode=figure_mode,chipIndex=chipSet)
          if(chipSets>1&chipSet!=chipSets){par(new=TRUE);upViewport();i=i+1}
        } 
      }else if(set == 'meth'){
        for (methSet in 1:methSets){
          if (methSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          plotMeth(methData, methTracksToPlot[methSet], plotMar, yLim=c(0,100),cex.axis=cex.axis,setNames=methNames[methSet],methColors=methColors[methSet],figure_mode=figure_mode,methIndex=methSet,currentIntrv=currentIntrv)
          if(methSets>1&methSet!=methSets){par(new=TRUE);upViewport();i=i+1}
        }  
      }else if(set == 'VP'){plotViewpoint(hicNames, plotLim, targetGene, extReads, plotMar)
      }else if(set == 'loops'){plotLoops(hicNames, plotLim, targetGene, extReads, plotMar)
      }else if(set == 'domains'){loadDomains(domainsToPlot, plotLim, currentIntrv, plotMar)
      }else if(set == 'anno'){plotSimpleBED(plotAnn,plotMar = plotMar,plotLim = plotLim)
      }else if(set == 'arcs'){
        arcIntervals <- subset(arcIntervals, arcIntervals$start1 > plotLim[1] & arcIntervals$end2 < plotLim[2])
        plotArcs(arcIntervals,plotMar=plotMar,color_by=arcColor_by,cols=arcColors,plotLim=plotLim,flip=F,plotScale=F)
      }else if(set == 'axis'){plotAxis(plotLim, plotMar,currentIntrv$chrom,cex.axis=cex.axis,figure_mode=figure_mode)
      }else if(set == 'ideogram'){plotIdeogram(genome,plotLim,currentIntrv$chrom)
      }else if(set == 'scATAC'){
        p <- scCoveragePlot(object=object,region=currentIntrv,window = scATAC_window,annotation=NULL,idents = plotClusters,cells=scCellsToPlot,fragment.path=fragment.path,plotMar=plotMar,cluster_cols = cluster_cols)
        print(p, vp=viewport(layout.pos.row = i, layout.pos.col = 1))
        # 	}else if(set == 'arcs'){plotArcs(arcReads, arcsToDraw, plotLim, widthCorrection=25, loopThr=50)
      }else{plot(0,xlim=plotLim,xlab='', ylab='', yaxt='n', xaxt='n', frame=T, main=set)}
      par(new=TRUE)
      if(set == 'genes'|set=='ideogram'|set=='scATAC'|set=='scores') {popViewport()}
      upViewport()
      i=i+1
    }  
  }
  
  dev.off()
}

plotSingleKNN_ggplot <- function(set1, plotAnn=F, plotMar=c(3.5,3), plotLim, plot2Dann=F, plotThr=-101, loopDetails=FALSE, radius=25e3,pointCEX=0.05,plotScale=TRUE,cex.axis=2,window_scale=1){
  if(is.null(set1)){return()}
  if(window_scale==1){
    set1 <- subset(set1, set1$start1 > plotLim[1] & set1$start2 < plotLim[2])
  } else {
    set1 <- subset(set1, (set1$start1 + set1$start2)/2 > plotLim[1] & (set1$start1 + set1$start2)/2 < plotLim[2])
  }
  set1 <- subset(set1, set1$start1 < set1$start2)
  set1 <- subset(set1, set1[,7] > plotThr)
  ########################################################
  ########################################################
  winSize <- plotLim[2]-plotLim[1]
  yLim <- c(0,winSize)
  # 	yLim <- c(0, max((set1$start2 - set1$start1)))
  ########################################################
  set1 <- set1[order(set1[,7]),]
  plot_label <- colnames(set1)[7]
  plot_label <- gsub('.score_k200','',plot_label)
  plot_label <- gsub('.score_k100','',plot_label)
  plot_label <- gsub('hic.','',plot_label)
  plot_label <- gsub('E14_','',plot_label)
  plot_label <- gsub('ncx_','',plot_label)
  plot_label <- gsub('_IUE24h','',plot_label)
  plot_label <- gsub('NGN2','Neurog2',plot_label)
  plot_label <- ''
  if(grepl('3DRAM',plot_label)|grepl('Methyl',plot_label)|grepl('Bonev',plot_label)){
    plot_label <- ''
  }
  plot_label <- ''
  colnames(set1)[7] <- 'score'
  plot_points <- data.frame(x=(set1$start1+set1$start2)/2,y=(set1$start2-set1$start1)/2,score=set1$score)
  p <- ggplot(plot_points,aes(x=x, y=y, color=score)) + theme(legend.position="none") + geom_point_rast(size=pointCEX,raster.dpi = 72) +
    coord_cartesian(xlim=plotLim,ylim=yLim/window_scale) 
  p <- p + scale_color_gradientn(colors=col.scores[1:201],limits=c(-100,100),breaks=c(-100,0,100),name='') +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_cowplot() + theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=1)) + annotate(geom='text',label=plot_label,x=plotLim[1]+0.95*(plotLim[2]-plotLim[1]),y=0.9*yLim[2]/window_scale,fontface =2,size=10)
  loopDetails <- intervals.2d.centers(loopDetails)
  plot_anno <- data.frame(x=(loopDetails$start1+loopDetails$start2)/2,y=(loopDetails$start2-loopDetails$start1)/2)
  plot_anno$radius <- radius
  p <- p + ggforce::geom_circle(aes(x0=x,y0=y,r=radius),data=plot_anno,color='black',inherit.aes=F,linetype=2) 
  return(p)
}


centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}



getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}

groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gs <- lapply(unique(groups), function(x){
    if (sparse){
      matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
    }else{
      matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gs) <- unique(groups)
  return(gs)
}


ggFootprint <- function(
  seFoot = NULL,
  name = NULL,
  pal = NULL,
  smoothWindow = NULL,
  flank = NULL,
  flankNorm = NULL,
  baseSize = 6,
  normMethod = NULL,
  which_clusters=NULL,
  aggregate_groups=FALSE
){
  
  errorList <- list()
  
  #Get Footprint Info
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="footprint"),], name)
  biasMat <- getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="bias"),], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
  biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]
  
  errorList$footMat <- footMat
  errorList$biasMat <- biasMat
  errorList$footDF <- footDF
  errorList$biasDF <- biasDF
  
  #Smooth Foot and Bias Mat because of sparsity
  if(!is.null(smoothWindow)&!(aggregate_groups)){
    footMat <- apply(footMat, 2, function(x) centerRollMean(x, smoothWindow))
    biasMat <- apply(biasMat, 2, function(x) centerRollMean(x, smoothWindow))
  }
  
  if(aggregate_groups){
    footMat <- cbind(rowSums(footMat[,grep('Rep1',colnames(footMat))]),rowSums(footMat[,grep('Rep2',colnames(footMat))]))
    colnames(footMat) <- c('Rep1','Rep2')
    biasMat <- cbind(rowSums(biasMat[,grep('Rep1',colnames(biasMat))]),rowSums(biasMat[,grep('Rep2',colnames(biasMat))]))
    colnames(biasMat) <- c('Rep1','Rep2')
  }
  
  
  #Normalize Foot and Bias Mat
  idx <- which(abs(footDF$x) >= flank - flankNorm)
  footMat <- t(t(footMat) / colMeans(footMat[idx, ,drop=FALSE]))
  biasMat <- t(t(biasMat) / colMeans(biasMat[idx, ,drop=FALSE]))
  
  errorList$footMatNorm <- footMat
  errorList$biasMatNorm <- footMat
  
  #Norm Foot By Bias
  if(tolower(normMethod) == "none"){
    title <- ""
  }else if(tolower(normMethod) == "subtract"){
    title <- "Tn5 Bias Subtracted\n"
    footMat <- footMat - biasMat
  }else if(tolower(normMethod) == "divide"){
    title <- "Tn5 Bias Divided\n"
    footMat <- footMat / biasMat
  }else{
    stop("normMethod not recognized!")
  }
  
  #Get Mean and SD for each Assay
  if (aggregate_groups){
    footMatMean <- groupMeans(footMat, colnames(footMat))
    footMatSd <- groupSds(footMat, colnames(footMat))
    biasMatMean <- groupMeans(biasMat, colnames(biasMat))
    biasMatSd <- groupSds(biasMat, colnames(biasMat))
    smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) centerRollMean(x, 11)))
  } else {
    footMatMean <- groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
    footMatSd <- groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
    biasMatMean <- groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
    biasMatSd <- groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
    smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) centerRollMean(x, 11)))
  }
  

  
  errorList$footMatMean <- footMatMean
  errorList$footMatSd <- footMatSd
  errorList$biasMatMean <- biasMatMean
  errorList$biasMatSd <- biasMatSd
  errorList$smoothFoot <- smoothFoot
  
  #Create Plot Data Frames
  plotIdx <- seq_len(nrow(footMatMean)) #sort(unique(c(1, seq(1, nrow(footMatMean), smoothWindow), nrow(footMatMean))))
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x){
    data.frame(
      x = footDF$x, 
      mean = footMatMean[,x], 
      sd = footMatSd[,x], 
      group = colnames(footMatMean)[x]
    )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))
  
  plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x){
    data.frame(
      x = biasDF$x, 
      mean = biasMatMean[,x], 
      sd = biasMatSd[,x], 
      group = colnames(biasMatMean)[x]
    )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))
  
  errorList$plotFootDF <- plotFootDF
  errorList$plotBiasDF <- plotBiasDF
  
  if(is.null(pal)){
    pal <- paletteDiscrete(values=gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
  }
  
  plotMax <- plotFootDF[order(plotFootDF$mean,decreasing=TRUE),]
  plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 50, ] #<= flank - flankNorm,]
  plotMax <- plotMax[!duplicated(plotMax$group),]
  plotMax <- plotMax[seq_len(ceiling(nrow(plotMax) / 4)), ]
  plotMax$x <- 25
  if(!aggregate_groups){plotFootDF$group <- factor(plotFootDF$group,levels=which_clusters)}
  ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) + 
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
    geom_line() + 
    scale_color_manual(name='',values = pal) + 
    scale_fill_manual(name='',values = pal) + 
    xlab("Distance to motif center (bp)") +
    coord_cartesian(
      expand = FALSE, 
      ylim = c(quantile(plotFootDF$mean, 0.0001), 1.15*quantile(smoothFoot, 0.999)), 
      xlim = c(-flank,flank)
    ) +
    guides(fill = FALSE) + 
    guides(color = FALSE) + ylab(paste0(title,"Normalized Insertions")) 
  
  if(!aggregate_groups){plotBiasDF$group <- factor(plotBiasDF$group,levels=which_clusters)}
  ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) + 
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
    geom_line() +
    scale_color_manual(name='',values = pal) + 
    scale_fill_manual(name='',values = pal) + 
    xlab("Distance to center (bp)") +
    coord_cartesian(
      expand = FALSE, 
      ylim = c(quantile(plotBiasDF$mean, 0.0001), 1.05*quantile(plotBiasDF$mean, 0.999)), 
      xlim = c(-flank,flank)
    ) + ylab("Tn5-Bias Normalized Insertions") + 
    theme(legend.position = "bottom", legend.box.background = element_rect(color = NA)) 
  return(list(ggFoot=ggFoot,ggBias=ggBias))
}
  


palette.breaks = function(n, colors, breaks){
  colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
  
  for(i in 2:(length(colors)) ){
    colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
  }
  colspec = c( colspec,
               colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
  )
  colspec
}

plotSingleKNN_ggplot <- function(set1, plotAnn=F, plotMar=c(3.5,3), plotLim, plot2Dann=F, plotThr=-101, loopDetails=FALSE, radius=25e3,pointCEX=0.05,plotScale=TRUE,cex.axis=2,window_scale=1){
  if(is.null(set1)){return()}
  if(window_scale==1){
    set1 <- subset(set1, set1$start1 > plotLim[1] & set1$start2 < plotLim[2])
  } else {
    set1 <- subset(set1, (set1$start1 + set1$start2)/2 > plotLim[1] & (set1$start1 + set1$start2)/2 < plotLim[2])
  }
  set1 <- subset(set1, set1$start1 < set1$start2)
  set1 <- subset(set1, set1[,7] > plotThr)
  ########################################################
  ########################################################
  winSize <- plotLim[2]-plotLim[1]
  yLim <- c(0,winSize)
  # 	yLim <- c(0, max((set1$start2 - set1$start1)))
  ########################################################
  set1 <- set1[order(set1[,7]),]
  plot_label <- colnames(set1)[7]
  plot_label <- gsub('.score_k200','',plot_label)
  plot_label <- gsub('.score_k100','',plot_label)
  plot_label <- gsub('hic.','',plot_label)
  plot_label <- gsub('E14_','',plot_label)
  plot_label <- gsub('ncx_','',plot_label)
  plot_label <- gsub('_IUE24h','',plot_label)
  plot_label <- gsub('NGN2','Neurog2',plot_label)
  plot_label <- ''
  colnames(set1)[7] <- 'score'
  plot_points <- data.frame(x=(set1$start1+set1$start2)/2,y=(set1$start2-set1$start1)/2,score=set1$score)
  p <- ggplot(plot_points,aes(x=x, y=y, color=score)) + theme(legend.position="none") + geom_point_rast(size=pointCEX,raster.dpi = 72) +
    coord_cartesian(xlim=plotLim,ylim=yLim/window_scale) 
  p <- p + scale_color_gradientn(colors=col.scores[1:201],limits=c(-100,100),breaks=c(-100,0,100),name='') +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_cowplot() + theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=1)) + annotate(geom='text',label=plot_label,x=plotLim[1]+0.95*(plotLim[2]-plotLim[1]),y=0.9*yLim[2]/window_scale,fontface =2,size=10)
  loopDetails <- intervals.2d.centers(loopDetails)
  plot_anno <- data.frame(x=(loopDetails$start1+loopDetails$start2)/2,y=(loopDetails$start2-loopDetails$start1)/2)
  plot_anno$radius <- radius
  #p <- p + ggforce::geom_circle(aes(x0=x,y0=y,r=radius),data=plot_anno,color='black',inherit.aes=F,linetype=2) 
  return(p)
  }

scCoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 1,
  cluster_cols=hue_pal(length(idents)),
  cells = NULL,
  idents = NULL,
  singleCells=TRUE,
  average=FALSE,
  plotMar=c(5,3)
) {
  if (!is.null(cells)){
    cells <- cells[cells %in% colnames(object)]
  } else {
    cells <- colnames(x = object)
  }
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (!is(object = region, class2 = 'GRanges')) {
    region <- makeGRangesFromDataFrame(region)
  }
  reads <- GetReadsInRegion(
    object = object,
    region = region,
    cells = cells,
    group.by = group.by,
    fragment.path = fragment.path,
    verbose = FALSE
  )
  cells.per.group <- CellsPerGroup(
    object = object,
    group.by = group.by
  )
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  
  if (average){
    reads.per.group <- AverageCounts(
      object = object,
      group.by = group.by,
      verbose = FALSE
    )
    coverages <- suppressWarnings(CalculateCoverages(
      reads = reads,
      cells.per.group = cells.per.group,
      reads.per.group = reads.per.group,
      window = window/2,
      verbose = FALSE
    ))
    if (downsample > 1) {
      warning('Requested downsampling <0%, retaining all positions')
      downsample <- 1
    }
    
    stepsize <- 1 / downsample
    total_range <- end.pos - start.pos
    steps <- ceiling(x = (total_range / stepsize))
    retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
    downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
    ymax <- signif(x = max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)
    downsampled_coverage <- downsampled_coverage[!is.na(x = downsampled_coverage$coverage), ]
    
    p <- ggplot(data = downsampled_coverage, mapping = aes(x = position, y = coverage, fill = group)) +
      geom_bar(stat = 'identity', width=1) +
      facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      ylab(label = paste0('0 - ', as.character(x = ymax))) +
      ylim(c(0, ymax)) +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        legend.position = 'none',
        strip.text.y = element_text(angle = 0)
      )
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  }
  
  if(singleCells){
    reads <- makeGRangesFromDataFrame(reads,keep.extra.columns = T)
    sliding_region <- slidingWindows(region,window,window)[[1]]
    mat <- t(countInsertions2(sliding_region,reads,cells=cells,by = 'cell',window=window))
    mat@x[mat@x > 0] <- 1
    df.long <- melt(as.matrix(mat))
    df.long$group <- Idents(object)[match(df.long$Var1,names(Idents(object)))]
    cells_group <- paste0(names(cells.per.group))
    names(cells_group) <- names(cells.per.group)
    cells_group <- ''
    plotLabel <- as.list(cells_group)
    plot_labeller <- function(variable,value){
      return(plotLabel[value])
    }
    df_sub <- subset(df.long, value == 1)
    df_sub$color <- df_sub$group
    if (sum(!(df.long$Var1%in%df_sub$Var1))>0){
      df_neg <- as.data.frame(df.long[!(df.long$Var1%in%df_sub$Var1),]%>% group_by(Var1) %>% dplyr::slice(1))
      df_neg$color <- NA
      df_comb <- rbind(df_sub,df_neg)
    } else {df_comb <- df_sub}
    
    df_comb <- df_comb[sample(row.names(df_comb),nrow(df_comb)),] %>% arrange(group) %>% mutate(Var1=factor(Var1, levels=unique(Var1)))
    row.names(df_comb) <- 1:nrow(df_comb)
    #df_comb <- df_comb[sample(1:nrow(df_comb),nrow(df_comb)),]
    p1 <- ggplot(data = df_comb, mapping = aes(x = Var2, y = Var1, fill = color,size=4)) +
      geom_tile() + scale_x_continuous(expand = c(0,0),limits=c(start(region),end(region))) +
      scale_fill_manual(name='',values=as.character(cluster_cols),na.value = 'white')+
      facet_grid(rows=vars(group), scales="free_y",space='free_y',labeller =plot_labeller,margins=F) +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line=element_blank(),
        legend.position = 'none',
        strip.text.y = element_text(hjust=-10),
        strip.background = element_rect(
          color="black", fill=cluster_cols, size=1, linetype="solid"
        ),
        plot.margin = unit(c(0,0,0,0), "cm")
      )
    p1 <- rasterize(p1, layers='Point', dpi=300)
    g <- ggplot_gtable(ggplot_build(p1))
    stripr <- which(grepl('strip-r', g$layout$name))
    fills <- cluster_cols
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    p1 <- g
  } 
  
  if (!is.null(x = annotation)) {
    gr <- GRanges(
      seqnames = chromosome,
      IRanges(start = start.pos, end = end.pos)
    )
    hits <- DataFrame(findOverlaps(gr,transcripts(annotation),type='any'))
    hits <- transcripts(annotation)[hits$subjectHits]
    hits <- AnnotationDbi::select(org.Mm.eg.db,keys=gsub('\\..*','',hits$tx_name),keytype = "REFSEQ",columns = "SYMBOL")
    genes <- suppressMessages(expr = ggbio::autoplot(object = annotation, which=gr,names.expr = hits$SYMBOL))
    gene.plot <- genes@ggplot +
      xlim(start.pos, end.pos) +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      theme_classic()
  }
  if (singleCells) {
    if (average){
      if(!is.null(x = annotation)){
        plotlist <- list(p1,p,gene.plot)
        rel.heights=c(2,1,0.6)
      }
      plotlist <- list(p1,p)
      rel.heights=c(2,1)
    }else if(!is.null(x = annotation)){
      plotlist <- list(p1,gene.plot)
      rel.heights=c(1,0.6)
    } else {
      plotlist <- list(p1)
      rel.heights=c(1)     
    }
  } else if(!is.null(x = annotation)){
    plotlist <- list(p,gene.plot)
    rel.heights=c(1,0.6)
  } else {
    plotlist <- list(p)
    rel.heights=c(1,0.6)    
  }
  
  p <- suppressWarnings(cowplot::plot_grid(
    plotlist = plotlist,
    ncol = 1,
    axis = 'btlr',
    rel_heights = rel.heights,
    align = 'v',
    greedy = FALSE
  ))
  # p <- p + theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  message('Succesfully extracted scATAC')
  return(p)
}
