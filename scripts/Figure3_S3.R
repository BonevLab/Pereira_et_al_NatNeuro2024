library(Seurat)
library(Signac)
library(ggplot2)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
library(Hmisc)
library(SummarizedExperiment)
require(LSD)
library(ggrepel)
library(ggpointdensity)
library(grid)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(tglkmeans)
library(monaLisa)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

## Set main folder as current dir ####
require(ArchR)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")
options(future.globals.maxSize = 100 * 1024 ^ 3)


source('scripts/figures/config.R')
source('scripts/figures/plot_functions.R')
source('scripts/SC/aux_functions.R')
source('scripts/aux_functions.R')

source('/home/hpc/bonev/projects/hic/dfg/config.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/main_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/aux_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))


object <- readRDS('results/SC/seurat_IDs.RDS')
cell_markers <- read.table('data/SC/cell_markers.tsv')
rna_fc <- read.table('data/RNA/res.Ngn2-GFP.txt')

p2glinks <- readRDS('results/P2G-Links.RDS')
seRNA_all <- readRDS("results/SC/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
gtfFile <- '/home/hpc/bonev/annotations/mm10/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf'
gtf <- getGeneGTF(gtfFile)
gtf <- unique(gtf[gtf$gene_name%in%row.names(seRNA_all)])
gene.length=gtf$exonLength[match(row.names(seRNA_all),gtf$gene_name)]

uf_rpkm <- edgeR::rpkmByGroup(assay(seRNA_all)[,c(1,2,5:8)],group=factor(c('AST','AST','iN_1','iN_1','iN_2','iN_2'),levels=c('AST','iN_1','iN_2')),
                              gene.length=gene.length,log=F,prior.count=0)


Figure3A <- function(archr_obj,group.by,useMatrix="PeakMatrix",which_features,fdr.cutoff=0.1,log2FC.cutoff=0.5,tss.coords,distance=5000,k=5,features,cols,cluster_cols,anno.size,out_f,width,height){
  markersPeaks <- getMarkerFeatures(
    ArchRProj = archr_obj,
    useGroups=which_features,
    useMatrix = useMatrix, 
    groupBy = group.by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  mat <- plotMarkerHeatmap(
    seMarker = markersPeaks, 
    cutOff = paste0("FDR <= ",fdr.cutoff," & Log2FC >= ",log2FC.cutoff),
    transpose = FALSE
    ,returnMatrix = T
  )
  if(is.data.frame(tss.coords)){
    colnames(tss.coords) <- c('chrom','start','end','gene','name','strand')
    tss.coords <- makeGRangesFromDataFrame(tss.coords,keep.extra.columns = T)
  }
  obj_features <- StringToGRanges(rownames(mat), sep = c(":", "-"))
  overlaps <- distanceToNearest(obj_features,tss.coords,ignore.strand=T)
  mat <- as.data.frame(mat)
  mat$nearestGene <- tss.coords$gene[overlaps@to]
  mat$dist <- start(resize(obj_features,width = 1,'center'))[overlaps@from] - start(tss.coords)[overlaps@to]
  ### Invert sign of distance for - strand genes ####
  mat$dist[as.vector(strand(tss.coords)=='-')[overlaps@to]] <- mat$dist[as.vector(strand(tss.coords)=='-')[overlaps@to]]*(-1)
  #Split into distal and promoter-associated
  mat_d <- mat[abs(mat$dist)>distance,]
  mat_p <- mat[(500>mat$dist)&(mat$dist>(-2000)),]
  mat_d$labels <- paste0(mat_d$nearestGene,':',mat_d$dist)
  mat_p$labels <- mat_p$nearestGene
  #K-means clustering
  mat_k_d <- TGL_kmeans(mat_d[,1:length(which_features)], k=k, id_column = FALSE,reorder_func = hclust,seed=42,max_iter=500,hclust_intra_clusters = F,parallel = F)
  mat_k_p <- TGL_kmeans(mat_p[,1:length(which_features)], k=k, id_column = FALSE,reorder_func = hclust,seed=42,max_iter=100,hclust_intra_clusters = F,parallel = F)
  #Reorder cluster levels based on max column
  clust_order_d <- order(max.col(mat_k_d$centers, "first"))
  clust_order_p <- order(max.col(mat_k_p$centers, "first"))
  mat_d$cluster <- factor(mat_k_d$cluster,levels=clust_order_d)
  levels(mat_d$cluster) <- 1:k
  mat_p$cluster <- factor(mat_k_p$cluster,levels=clust_order_p)
  levels(mat_p$cluster) <- 1:k
  #Create hm
  n_clusters <- length(which_features)
  col.list <- cluster_cols
  names(col.list) <- colnames(mat)[1:n_clusters]
  ha <- columnAnnotation(foo = anno_text(colnames(mat)[1:n_clusters], location = 0.5,rot=0, just = "center",gp = gpar(fill = col.list, col = "black", border = "black", fontsize = 12),height = unit(height/50,'inch')))
  la <- rowAnnotation(foo = anno_block(gp = gpar(fill = rev(colorpalette('brbg',k)))))
  ra_d = rowAnnotation(foo = anno_mark(at = which(mat_d$labels%in%features), labels = mat_d$labels[which(mat_d$labels%in%features)],labels_gp = gpar(fontsize = anno.size)))
  ra_p = rowAnnotation(foo = anno_mark(at = which(mat_p$labels%in%features), labels = mat_p$labels[which(mat_p$labels%in%features)],labels_gp = gpar(fontsize = anno.size)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm_d <- Heatmap(as.matrix(mat_d[,1:n_clusters]),name='Distal',row_title = paste0('Distal: ',nrow(mat_d)),row_split = mat_d$cluster,show_row_names = F,cluster_row_slices = F,cluster_rows = F,left_annotation = la,show_row_dend = F,cluster_columns = F,show_column_names = F,use_raster = TRUE, raster_quality = 10, right_annotation = ra_d,col = cols,top_annotation = ha,
                  border=T,height=unit(height*0.65, "inch"),heatmap_legend_param=list(at=c(-1.2,0,1.2),direction = "horizontal",title='Row Z-scores',legend_width = unit(width/6, "inch")))
  hm_p <- Heatmap(as.matrix(mat_p[,1:n_clusters]),name='Promoter',row_title = paste0('Promoter: ',nrow(mat_p)),row_split = mat_p$cluster,show_row_names = F,cluster_row_slices = F,cluster_rows = F,left_annotation = la,show_row_dend = F,cluster_columns = F,show_column_names = F,use_raster = TRUE, raster_quality = 10, right_annotation = ra_p,col = cols,
                  border=T,height=unit(height*0.25, "inch"),show_heatmap_legend = F)
  hm <- hm_d %v% hm_p
  draw(hm,main_heatmap='Distal',merge_legends=T, ht_gap = unit(0.1, "inch"),heatmap_legend_side = "bottom", 
       annotation_legend_side = "bottom") 
  dev.off()
  return(list(all=mat,Distal=mat_d,Promoter=mat_p))
}

Figure3B <- function(res,which_features,pwm_f,fpkm_f,genome,mcparams,FDR.cutoff=1,enr.cutoff=0.25,cols){
  mat <- res[[which_features]]
  pwms <- readRDS(pwm_f)
  ##Filter pwms for expressed motifs 
  motif.names <- name(pwms)
  motif.names <- capitalize(tolower(motif.names))
  simple_motif_names <- sapply(strsplit(motif.names, "\\(|\\:|\\."), function(x) x[[1]])
  fpkm_mat <- read.table(fpkm_f,header = T,sep='\t')
  motif_idx <- simple_motif_names%in%row.names(fpkm_mat)[rowMaxs(as.matrix(fpkm_mat[,1:4]),na.rm=T)>=1]
  pwms <- pwms[motif_idx,]
  ######
  pwms <- toPWM(pwms)
  pwms <- pwms[!duplicated(names(pwms))]
  peaks <- featureToGR(gsub('\\:','-',row.names(test)),pattern = '-')
  names(peaks) <- test$labels
  peakseqs <- getSeq(genome, peaks)
  se <- calcBinnedMotifEnrR(seqs = peakseqs, bins = factor(test$cluster),pwmL = pwms,BPPARAM = mcparams)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  seSel <- se[sel1&sel2, ]
  #hcl <- hclust(dist(assay(seSel,"log2enr")), method = "ward.D")
  
  df <- assay(seSel,"log2enr")
  row.names(df) <- rowData(seSel)$motif.name
  df <- df[ order(max.col(df, "first")), ]
  col.list=rev(colorpalette('brbg',max(as.numeric(mat$cluster))))
  names(col.list) <- 1:max(as.numeric(mat$cluster))
  if(!is.null(features)){
    la = rowAnnotation(foo = anno_mark(at = which(row.names(df)%in%features), labels = row.names(df)[which(row.names(df)%in%features)],labels_gp = gpar(fontsize = anno.size)))
  } else {
    la = rowAnnotation(foo = anno_mark(at = 1:nrow(df), labels = row.names(df),labels_gp = gpar(fontsize = anno.size)))
  }
  ta <- HeatmapAnnotation(cluster = names(col.list) ,col = list(cluster=col.list),show_legend = F,show_annotation_name=F,border=T)
  hm <- Heatmap(df,name='log2enr',cluster_rows = F,show_row_names = F,right_annotation = la,cluster_columns = F,show_column_names = F,show_row_dend = F,use_raster = TRUE, raster_quality = 10, col = cols,top_annotation = ta,border=T,heatmap_legend_param=list(at=c(-1.2,0,1.2)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  draw(hm)
  dev.off()
}

Figure3C <- function(object,out_f,features,n_cols,height,width){
  motifs <- object@assays$peaks@motifs@motif.names
  motifs <- as.data.frame(do.call("rbind", motifs))
  motifs$motifs <- row.names(motifs)
  colnames(motifs)[1] <- 'names'
  plot_motifs <- motifs[match(features,motifs$names),]
  plot_motifs <- plot_motifs[!is.na(plot_motifs$motifs),]
  p1 <- MotifPlot(
    object = object,assay = 'peaks',
    motifs = as.character(plot_motifs$motifs),
    use.names=T,ncol=n_cols)
  p1 <- p1 +theme_cowplot(font_size = 12) + theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line=element_blank())
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
}

Figure3D <- function(mat,nVar,object,which_clusters,which_conditions,features,clust_method='ward.D',cols,cluster_cols,anno.size,out_f,width,height,expr=NULL){
  motif.names <- object@assays$peaks@motifs@motif.names
  motif.names <- as.data.frame(do.call("rbind", motif.names))
  motif.names$motifs <- row.names(motif.names)
  row.names(mat) <- motif.names$V1[match(row.names(mat),motif.names$motifs)]
  if (!is.null(expr)){
    tf_name <- sapply(strsplit(as.character(row.names(mat)), "\\(|\\:|\\."), function(x) x[[1]])
    tf_name <- capitalize(tolower(tf_name))
    expr <- expr[match(tf_name,row.names(expr)),]
    mat <- mat[complete.cases(expr),]
  }
  sub_mat <- mat[row.names(mat)%in%features,]
  mat <- head(mat[order(matrixStats::rowVars(mat), decreasing = TRUE),],nVar)
  mat <- rbind(mat,sub_mat[!row.names(sub_mat)%in%row.names(mat),])
  
  idents_indx1 <- Idents(object)
  idents_indx1 <- idents_indx1[match(colnames(mat),names(idents_indx1))]
 
  idents_indx2 <- object$orig.ident
  idents_indx2 <- idents_indx2[match(colnames(mat),names(idents_indx2))]
  
  intraMean1 <- groupMeans(mat, groups = idents_indx1, sparse = F, na.rm = TRUE)
  intraMean2 <- groupMeans(mat, groups = idents_indx2, sparse = F, na.rm = TRUE)
  intraMean2 <- intraMean2
  test <- TGL_kmeans(cbind(intraMean1,intraMean2),k = 5,id_column = F,seed = 42)
  rowClust <- hclust(dist(cbind(intraMean1,intraMean2),method='euclidean'),method = clust_method)
  mat <- mat[rowClust$order,]
  intraMean1 <- intraMean1[rowClust$order,]
  intraMean2 <- intraMean2[rowClust$order,]
  
  la = rowAnnotation(foo = anno_mark(at = which(row.names(intraMean1)%in%(features)), labels = row.names(intraMean1)[which(row.names(intraMean1)%in%unique(c(features)))],labels_gp = gpar(fontsize = anno.size)))
  names(cluster_cols) <- which_clusters
  names(condition_cols) <- which_conditions
  
  ha1 = HeatmapAnnotation(cluster = colnames(intraMean1) ,col = list(cluster=cluster_cols),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical')))
  ha2 = HeatmapAnnotation(condition = colnames(intraMean2) ,col = list(condition=condition_cols),show_legend = T,show_annotation_name=F,annotation_legend_param=list(condition=list(direction='vertical')))
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm1 <- Heatmap(intraMean1,cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = T,use_raster = TRUE, raster_quality = 3, left_annotation = la,col = cols,top_annotation = ha1,
                heatmap_legend_param=list(direction = "vertical",title='',legend_height = unit(height/6, "inch")))
  hm2 <- Heatmap(intraMean2,cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = T,use_raster = TRUE, raster_quality = 3, col = cols,top_annotation = ha2,
                heatmap_legend_param=list(direction = "vertical",title='',legend_height = unit(height/6, "inch")))
    draw(hm1+hm2, merge_legends=T) 
  dev.off()
  return(list(intraMean1,intraMean2))
}

Figure3E <- function(object,motifs,assay='chromvar',out_f,cols2,point.size,height,width,plot_filled=F,theme=NULL,stroke=0.1){
  plot_list2 <- list()
  for (i in seq_along(motifs)){
    DefaultAssay(object) <- assay
    feature_f <- motifs[i]
    DefaultAssay(object) <- assay
    motif.names <- object@assays$peaks@motifs@motif.names
    motif.names <- as.data.frame(do.call("rbind", motif.names))
    motif.names$motifs <- row.names(motif.names)
    if(motifs[i]%in%motif.names$V1){
      plot.data <- extractFeatures(object,features=motif.names$motifs[motif.names$V1==motifs[i]],cells_toInclude = c('all'),cells_toExclude = 'none')
    } else {
      plot.data <- extractFeatures(object,features=motifs[i],cells_toInclude = c('all'),cells_toExclude = 'none')
    }
    #cols_lower <- colorpalette(cols2,round(abs(min(plot.data[, 3],na.rm=T)))*2,rev=F)[1:round(abs(min(plot.data[, 3],na.rm=T)))]
    #cols_upper <- colorpalette(cols2,round(abs(max(plot.data[, 3],na.rm=T)))*2,rev=F)[round(abs(max(plot.data[, 3],na.rm=T))+1):(round(abs(max(plot.data[, 3],na.rm=T)))*2)]
    cols_lower <- colorRampPalette(ArchR::ArchRPalettes$solarExtra)(round(abs(min(plot.data[, 3],na.rm=T)))*2)[1:round(abs(min(plot.data[, 3],na.rm=T)))]
    cols_upper <- colorRampPalette(ArchR::ArchRPalettes$solarExtra)(round(abs(max(plot.data[, 3],na.rm=T)))*2)[round(abs(max(plot.data[, 3],na.rm=T))+1):(round(abs(max(plot.data[, 3],na.rm=T)))*2)]
    cols <- c(cols_lower,cols_upper)
    #cols <- colorRampPalette(ArchRPalettes$solarExtra)
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_color_gradientn(colors=cols,name='',breaks=c(min(plot.data$feature),0,max(plot.data$feature)),labels=c(round(min(plot.data$feature)),0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 0.975),legend.background =element_blank(),legend.key = element_blank())
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colors=cols,name='',breaks=c(min(plot.data$feature),0,max(plot.data$feature)),labels=c(round(min(plot.data$feature)),0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.25, 0.95),legend.background =element_blank(),legend.key = element_blank())
      p <- p + guides(fill = guide_colourbar(barwidth = 4, barheight = 1))
    }
    if(!is.null(theme)){
      p <- p + theme_border}
    plot_list2[[feature_f]] <- p + ggtitle(paste0(feature_f))
  }
  p2 <- Reduce("+", plot_list2)
  width=length(motifs)*width
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p2)
  dev.off()
}

Figure3G <- function(selected_peaks,scaleMax=NULL,sortBy='ATAC',seRNA_all,sePB_all,out_f,features,padj=1,minLFC=0,minSdRatio=0,zCutoff=0,which_clusters,cluster_cols,scale=T,cols1,cols2,anno.size=12,height,width){
  seRNA <- seRNA_all[row.names(seRNA_all)%in%selected_peaks$gene_name,]
  sePB <- sePB_all[match(gsub('_','-',selected_peaks$peakName),row.names(sePB_all)),]
  uf <- suppressMessages(uniqueFeatures(
    edgeR::cpm(assay(sePB),log=TRUE,prior.count=1),
    groups = colData(sePB)$Group,
    padj = padj,
    minSdRatio = minSdRatio,
    minLFC = minLFC,
    zCutoff = zCutoff,
    breakPt = "last",
    groupMin = 20,
    maxGroupSize = 2,
    clusterCols = F
  ))
  uf_rna <- suppressMessages(uniqueFeatures(
    edgeR::cpm(assay(seRNA),log=TRUE,prior.count=1),
    groups = colData(seRNA)$Group,
    padj = padj,
    minSdRatio = minSdRatio,
    minLFC = minLFC,
    zCutoff = zCutoff,
    breakPt = "last",
    groupMin = 0,
    maxGroupSize = 2,
    clusterCols = F
  ))
  mat <- uf$groupMat
  indx <- rep(TRUE,nrow(mat))
  if (length(unique(which_clusters))!=length(unique(colnames(mat)))){
    mat <- mat[,colnames(mat)%in%which_clusters]
    bin_mat <- uf$binaryMat
    if(sortBy=='ATAC'){
      indx <- rowMaxs(bin_mat[,colnames(bin_mat)%in%which_clusters])>0
      mat <- mat[indx,]
    }
  }
  if (scale) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
  }
  if (is.numeric(scaleMax)){
    mat[mat > scaleMax] <- scaleMax
    mat[mat < -scaleMax] <- -scaleMax
  }
  n_clusters <- ncol(mat)
  mat_r <- uf_rna$groupMat
  if (length(unique(which_clusters))!=length(unique(colnames(mat_r)))){
    mat_r <- mat_r[,colnames(mat_r)%in%which_clusters]
    bin_mat_r <- uf_rna$binaryMat
    if(sortBy=='RNA'){
      indx_r <- rowMaxs(bin_mat_r[,colnames(bin_mat_r)%in%which_clusters])>0
      mat_r <- mat_r[indx_r,]
    }
  }
  mat_r_orig <- mat_r
  if (scale) {
    mat_r <- sweep(mat_r - rowMeans(mat_r), 1, matrixStats::rowSds(mat_r), `/`)
  }
  if (is.numeric(scaleMax)){
    mat_r[mat_r > scaleMax] <- scaleMax
    mat_r[mat_r < -scaleMax] <- -scaleMax
  }
  n_clusters_r <- ncol(mat_r)
  linksSig_IDs <- as.data.frame(cbind(as.vector(seqnames(selected_peaks)),as.numeric(start(selected_peaks)),as.numeric(end(selected_peaks ))))
  selected_peaks$IDs <- unite(linksSig_IDs,col='peaks_ID',sep='-')[,1]
  anno <- as.data.frame(mcols(selected_peaks),stringsAsFactors=FALSE)
  row.names(anno) <- paste0(anno$gene_name,':',anno$distance)
  anno <- anno[uf$rowOrder,]
  row.names(mat) <- paste0(anno$gene_name,':',anno$distance)
  anno$rna_order <- match(anno$gene_name,rownames(uf_rna$groupMat))
  anno$atac_order <- match(anno$peakName,gsub('-','_',rownames(uf$groupMat)))
  indx <- (!is.na(anno$rna_order))&(!is.na(anno$atac_order))
  comb_mat <- mat_r[match(anno$gene_name,row.names(mat_r)),]
  anno <- anno[indx,]
  mat <- mat[indx,]
  comb_mat <- comb_mat[indx,]
  if(sortBy=='ATAC'){
    indx_order <- order(anno$atac_order,anno$rna_order)
  } else {
    indx_order <- order(anno$rna_order,anno$atac_order)
  }
  anno <- anno[indx_order,]
  mat <- mat[indx_order,]
  comb_mat <- mat_r[match(anno$gene_name,row.names(mat_r)),]
  mat <- mat[complete.cases(comb_mat),]
  comb_mat <- comb_mat[complete.cases(comb_mat),]
  anno_orig <- anno
  col.list1 <- cluster_cols[1:n_clusters]
  names(col.list1) <- colnames(mat) 
  col.list2 <- cluster_cols[1:n_clusters_r]
  names(col.list2) <- colnames(comb_mat)
  ha1 = HeatmapAnnotation(ATAC_clusters = colnames(mat[,1:n_clusters]) ,col = list(ATAC_clusters=col.list1),show_legend = T,show_annotation_name=F,annotation_legend_param=list(ATAC_clusters=list(direction='vertical')))
  ha2 = HeatmapAnnotation(RNA_clusters = colnames(comb_mat[,1:n_clusters_r]) ,col = list(RNA_clusters=col.list2),show_legend = T,show_annotation_name=F,annotation_legend_param=list(RNA_clusters=list(direction='vertical')))
  ha1 = HeatmapAnnotation(foo = anno_mark(at = 1:n_clusters, labels = colnames(mat[,1:n_clusters]),labels_gp = gpar(fontsize = 12),labels_rot = 0),ATAC_clusters=colnames(mat[,1:n_clusters]),col = list(ATAC_clusters=col.list1),show_legend = F,show_annotation_name=F,annotation_legend_param=list(ATAC_clusters=list(direction='vertical')))
  ha2 = HeatmapAnnotation(foo = anno_mark(at = 1:n_clusters_r, labels = colnames(comb_mat[,1:n_clusters_r]),labels_gp = gpar(fontsize = 12),labels_rot = 0),RNA_clusters=colnames(comb_mat[,1:n_clusters_r]),col = list(RNA_clusters=col.list2),show_legend = F,show_annotation_name=F,annotation_legend_param=list(RNA_clusters=list(direction='vertical')))
  la1 <- columnAnnotation(foo = anno_text(colnames(mat), location = 0.5,rot=0, just = "center",gp = gpar(fill = col.list1, col = "black", border = "black", fontsize = 12),height = unit(height/50,'inch')))
  la2 <- columnAnnotation(foo = anno_text(colnames(comb_mat), location = 0.5,rot=0, just = "center",gp = gpar(fill = col.list2, col = "black", border = "black", fontsize = 12),height = unit(height/50,'inch')))
  if(!is.null(features)){
    ra1 = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%features),side='left', labels = row.names(mat)[which(row.names(mat)%in%features)],labels_gp = gpar(fontsize = anno.size)))
  } else {
    ra1 = NULL
  }
  hm1 <- Heatmap(mat,name='Accessibility',column_title = 'Accessibility',cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10, left_annotation = ra1,col = cols1,top_annotation = la1,
                 heatmap_legend_param=list(direction = "horizontal",title='Accessibility Z-score',legend_width = unit(width/6, "inch")))
  hm2 <- Heatmap(comb_mat,name='Expression',column_title = 'Expression',cluster_rows = F,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10,col = cols2,top_annotation = la2,
                 heatmap_legend_param=list(direction = "horizontal",title='Expression Z-score',legend_width = unit(width/6, "inch")))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  draw(hm1+hm2, merge_legends=F, ht_gap = unit(0.2, "inch"),heatmap_legend_side='bottom',annotation_legend_side = "bottom") 
  dev.off()
  if(sortBy=='RNA'){
    res <- merge(anno,uf_rna$binaryMat,by.x='gene_name',by.y='row.names',sort=F)
  } else {
    row.names(uf$binaryMat) <- gsub('-','_',row.names(uf$binaryMat))
    res <- merge(res,uf$binaryMat,by.x='peakName',by.y='row.names',sort=F)
  }
  return(res)
}

Figure3H <- function(df,control_df,include_prom=T,file_f,object,target_TF,cols,features=NULL,chip_f=NULL,rna_fc,GO=FALSE,anno.size=12,point.size=1,leg.text.size=10,width,height,n_go){
  expr_matrix <- GetAssayData(object,slot = 'data',assay='RNA')
  indx <- sapply(strsplit(as.character(target_TF), "\\(|\\:|\\."), function(x) x[[1]])
  indx <- capitalize(tolower(indx))
  motif.names <- object@assays$peaks@motifs@motif.names
  motif.matrix <- object@assays$peaks@motifs@data
  colnames(motif.matrix) <- as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])
  if(include_prom){
    df_prom <- df
    df_prom$Correlation <- 1
    ranges(df_prom) <- IRanges(start=df$gene_start-2000,end=df$gene_start+200)
    df_prom <- df_prom[!duplicated(df_prom$gene_name)]
    all_features <- StringToGRanges(row.names(object), sep = c("-", "-"))
    prom_overlaps <- findOverlaps(df_prom,all_features)
    #   prom_overlaps <- prom_overlaps[!is.na(prom_overlaps@to)]
    df_prom <- df_prom[prom_overlaps@from]
    ranges(df_prom) <- ranges(all_features[prom_overlaps@to])
    df_prom$peakName <- paste0(seqnames(df_prom),'_',start(df_prom),'_',end(df_prom))
    df <- c(df,df_prom)
  }
  meta.features <- GetAssayData(object = object, slot = 'meta.features')
  background <- MatchRegionStats(
    meta.feature = meta.features,
    query.feature = meta.features[match(gsub('_','-',df$peakName),row.names(meta.features)),],
    n = length(df)
  )
  if (!is.null(chip_f)){
    chip <- read.table(chip_f)
    colnames(chip)[1:3] <- c('chrom','start','end')
    df_features <- StringToGRanges(df$peakName, sep = c("_", "_"))
    chip_features <- makeGRangesFromDataFrame(chip)
    overlaps <- GenomicRanges::countOverlaps(df_features,chip_features)
    sel_df <- cbind(mcols(df),overlaps)
    colnames(sel_df)[ncol(sel_df)] <- target_TF
    linkage_df <- sel_df[sel_df[,target_TF]==1,]
    if(!is.null(control_df)){
      control_df <- control_df[control_df$Correlation<0.35&(control_df$Correlation>(-0.35))]
      control_df <- control_df[control_df$gene_name%in%df$gene_name,]
      control_features <- StringToGRanges(control_df$peakName, sep = c("_", "_"))
      control_overlaps <- GenomicRanges::countOverlaps(control_features,chip_features)
      control_sel_df <- cbind(mcols(control_df),control_overlaps)
      colnames(control_sel_df)[ncol(control_sel_df)] <- target_TF
    } else {
      control_df <- background
      control_features <- StringToGRanges(background, sep = c("-", "-"))
      control_sel_df <- as.data.frame(x =GenomicRanges::countOverlaps(control_features,chip_features))
      colnames(control_sel_df)[ncol(control_sel_df)] <- target_TF
    }
  } else {
    if(!is.null(control_df)){
      message('Calculating significance using provided control peaks')
      control_df <- control_df[control_df$Correlation<0.35&(control_df$Correlation>(-0.35))]
      control_df <- control_df[control_df$gene_name%in%df$gene_name,]
      control_motif.matrix <- motif.matrix[match(gsub('_','-',control_df$peakName),row.names(motif.matrix)),]
      control_sel_df <- control_motif.matrix
      control_sel_df <- cbind(mcols(control_df),control_sel_df[,target_TF,drop=F])
    } else {
      message('Calculating significance using GC-matched peaks')
      control_sel_df <- motif.matrix[match(background,row.names(motif.matrix)),]
    }
    motif.matrix <- motif.matrix[match(gsub('_','-',df$peakName),row.names(motif.matrix)),]
    sel_df <- cbind(mcols(df),motif.matrix[,target_TF,drop=F])
    linkage_df <- sel_df[sel_df[,target_TF]==1,]
  }
  
  linkageScore <- split(linkage_df$Correlation^2, as.character(as.vector(linkage_df$gene_name))) %>% 
    lapply(., sum) %>% 
    unlist(use.names=TRUE)
  
  query.counts <- split(as.numeric(sel_df[,target_TF]), as.character(as.vector(sel_df$gene_name))) %>% 
    lapply(., sum) %>% 
    unlist(use.names=TRUE) 
  
  query.length <- split(as.numeric(sel_df[,target_TF]), as.character(as.vector(sel_df$gene_name))) %>% 
    lapply(., length) %>% 
    unlist(use.names=TRUE) 
  
  
  plot_df <- as.data.frame(cbind(linkageScore=linkageScore[match(names(query.counts),names(linkageScore))],query.counts,query.length,background.counts=sum(control_sel_df[,target_TF]),background.length=nrow(control_sel_df),Pvalue=NA,rho=NA))
  row.names(plot_df) <- names(query.counts)
  plot_df$Gene <- row.names(plot_df)
  plot_df <- plot_df[!is.na(plot_df$linkageScore),]
  plot_df <- adply(plot_df,1,function(x){
    #   dat <- data.frame(Background=c(x[1,4],x[1,5] - x[1,4]),Query=c(x[1,2],x[1,3]-x[1,2]))        ###  Hypergeometric test is more appropriate here
    #   pvalue <- -log10(fisher.test(dat,alternative='less')$p.value)                                ###  Hypergeometric test is more appropriate here
    cor <- cor.test(expr_matrix[row.names(x),],expr_matrix[indx,],method = 'pearson')$estimate
    pvalue <- -log10(phyper(
      q = x[1,2],
      m = x[1,4],
      n = x[1,5] - x[1,4],
      k = x[1,3],
      lower.tail = FALSE
    ))
    return(c(rho=as.numeric(cor),Pvalue=pvalue)) 
  },.parallel = T)
  plot_df$Pvalue[is.infinite(plot_df$Pvalue)] <- NA
  plot_df <- plot_df[order(plot_df$Pvalue,decreasing=T),]
  plot_df$FC=rna_fc$log2FoldChange[match(plot_df$Gene,rna_fc$gene_symbol)]
  #### Plot the ranked plot labeling top 10 genes ########
  p <- ggplot(plot_df,aes(x=linkageScore,y=Pvalue))  + geom_point(color='grey',size = point.size*1.1) + geom_point(aes(fill = FC),alpha=0.8, pch = I(21),size = point.size,data = plot_df[plot_df$Pvalue>=(-log10(0.01)),]) +ggtitle(label = target_TF) + ylim(c(0,max(plot_df$Pvalue)*1.1))
  p <- p + scale_fill_gradientn(colours = cols,name='',breaks=c(-5,5),labels=c('-5','+5'),limits=c(-5,5)) + xlab('linkageScore') + ylab(expression(-Log[10]~(P)))
  p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.75, 0.97),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.2), vjust = 9,size=leg.text.size))+ guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
  if (is.null(features)){
    sel_df <- plot_df[plot_df$rho>0.1&plot_df$Pvalue>1.32,]
    features <- head(sel_df[order(sel_df$linkageScore,decreasing=T),'Gene'],10)
  }
  p <- p + ggrepel::geom_text_repel(
    data = plot_df, size = anno.size,box.padding = 0.8, min.segment.length = 0,max.iter = 10000,
    aes(x=linkageScore,y=Pvalue,color=NULL,label=ifelse(Gene%in%features, as.character(Gene), "")),force=10)
  
  if (GO){
    #############################
    ##### GO ANALYSIS ###########
    #############################
    
    geneList <- plot_df$Gene[plot_df$Pvalue>(-log10(0.05))]
    gene.df <- bitr(geneList, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene         = geneList,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    ego <- simplify(ego)
    p1 <- dotplot(ego, showCategory=n_go,) 
    p1$theme <- theme_cowplot()
    p1 <- p1 + scale_color_gradient(low='red',high='blue',breaks=scales::pretty_breaks(n=3),labels=scales::label_scientific())
    p1 <- p1 + scale_size(breaks=scales::pretty_breaks(n=3)) + scale_y_discrete(labels=label_wrap(40))
    p1 <- p1 + theme(legend.position = c(0.8,0.35)) 
    
    # pathview::pathview(gene.data  = gene.df$ENTREZID,
    #          pathway = kk@result$ID[2],
    #          species    = "mmu")
    p <- p+ p1 + plot_layout(ncol=2,widths=c(2.5,2))
    #########################
  }
  
  pdf(paste0('plots/figures/',file_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
  return(list(plot_df,p))
}





FigureS3A <- function(archr_obj,group.by,which_features,features,fdr.cutoff=0.1,log2FC.cutoff=0.5,tss.coords,point.size=1,distance=5000,cols,anno.size,out_f,width,height,xlim=NULL,ylim=NULL){
  markersPeaks <- getMarkerFeatures(
    ArchRProj = archr_obj,
    useGroups=which_features[2],
    bgdGroups = which_features[1],
    useMatrix = "PeakMatrix", 
    groupBy = group.by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  mat <- data.frame(row.names=paste0(rowData(markersPeaks)$seqnames,':',rowData(markersPeaks)$start,'-',rowData(markersPeaks)$end),
                    Log2FC=assay(markersPeaks,'Log2FC')[,1],FDR=assay(markersPeaks,'FDR')[,1])
  if(is.data.frame(tss.coords)){
    colnames(tss.coords) <- c('chrom','start','end','gene','name','strand')
    tss.coords <- makeGRangesFromDataFrame(tss.coords,keep.extra.columns = T)
  }
  obj_features <- StringToGRanges(rownames(mat), sep = c(":", "-"))
  overlaps <- distanceToNearest(obj_features,tss.coords,ignore.strand=T)
  mat <- as.data.frame(mat)
  mat$nearestGene <- tss.coords$gene[overlaps@to]
  mat$dist <- start(resize(obj_features,width = 1,'center'))[overlaps@from] - start(tss.coords)[overlaps@to]
  ### Invert sign of distance for - strand genes ####
  mat$dist[as.vector(strand(tss.coords)=='-')[overlaps@to]] <- mat$dist[as.vector(strand(tss.coords)=='-')[overlaps@to]]*(-1)
  #Split into distal and promoter-associated
  mat$class <- 'Distal'
  mat$class[(500>mat$dist)&(mat$dist>(-2000))] <- 'Prom'
  mat$labels <- paste0(mat$nearestGene,':',mat$dist)
  mat$labels[mat$class=='Prom'] <- mat$nearestGene[mat$class=='Prom']
  p <- ggplot(mat,aes(x=Log2FC,y=-log10(FDR))) + geom_point(color='grey',size = point.size) + geom_point(fill=cols[1],alpha=0.8, pch = I(21),size = point.size,data = mat[mat$FDR<fdr.cutoff&mat$Log2FC<(log2FC.cutoff*(-1)),]) + geom_point(fill=cols[2],alpha=0.8, pch = I(21),size = point.size,data = mat[mat$FDR<fdr.cutoff&mat$Log2FC>log2FC.cutoff,])  
  p <- p + xlab(expression(Log[2]~Fold~Change)) + ylab(expression(-Log[10]~(P)))
  p <- p + ggrepel::geom_text_repel(
    data = mat, size = anno.size,seed = 42,
    box.padding =0.8, min.segment.length = 0,max.iter = 10000,max.overlaps = Inf,
    aes(x=Log2FC,y=-log10(FDR),color=NULL,label=ifelse(labels%in%features, as.character(labels), "")),force=10)
  p <- p + ggtitle(paste0(which_features[1],' vs ',which_features[2]),subtitle=paste0(' Up: ',nrow(mat[mat$FDR<fdr.cutoff&mat$Log2FC>log2FC.cutoff,]),' Down:',nrow(mat[mat$FDR<fdr.cutoff&mat$Log2FC<(log2FC.cutoff*(-1)),])))
  p <- ggrastr::rasterise(p)
  if(!is.null(ylim)){
    p <- p+ ylim(ylim)
  }
  if(!is.null(xlim)){
    p <- p+ xlim(xlim)
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
  return(mat)
}


FigureS3B <- function(object,features,cols,out_f,height,width){
  motif.names <- unlist(object@assays$peaks@motifs@motif.names)
  motif.features <- names(motif.names)[motif.names%in%features]
  plot.data <- extractFeatures(object,features=motif.features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$feature[is.na(plot.data$feature)] <- 0
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = as.character(cols)) + geom_boxplot(width=0.1, fill="white",outlier.shape = NA) + theme(legend.position=c('none')) + xlab('') + ylab('ChromVar Deviation') + ggtitle(features)
  #return(p)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}




### Plot Figures

require(ArchR)
archr_obj <- loadArchRProject("data/SC/ArchR/seurat_atac/",showLogo = F)

res <- Figure3A(archr_obj=archr_obj,group.by='Clusters',which_features = c('AST','iN_1','iN_2'),
                tss.coords=read.table('data/mm10_TSS.bed'),
                features=c('Gfap','Aldoc','Fos','Tead3','Dcx','Rbfox3',"Hes5",
                           "Rbfox3:148916","Dab1:31912","Fos:-43134","Sox11:43842","Sox11:-14591","Sox11:-170662","Sox11:-495844","Sox11:-589737",
                           "Slc1a3:50970","Slc1a3:-72546","Slc1a3:-105947","Pcdh9:696185","Robo2:-115687"),
                cols=colorRamp2(seq(-1.2,1.2,length.out=length(heatmap_colors)), heatmap_colors),
                cluster_cols=cluster_colors[c(1,3,4)],anno.size=12,out_f='Figure3A',width=6,height=10)

res <- Figure3B(res=res,which_features = 'Distal',genome=BSgenome.Mmusculus.UCSC.mm10,
                fpkm_f='data/RNA/FPKM_AVGcounts_annot.txt',pwm_f='results/SC/combined_pwm.RDS',
                features=c("NEUROG2(var.2)","FOS::JUN","RFX4","IRF9","TGIF2","MEIS2","TEAD3","NEUROD2","ZEB1","TCF12(var.2)","JPD2","ZBTB18","HES6"),
                FDR.cutoff=2,enr.cutoff=0.5,mcparams=BiocParallel::MulticoreParam(10L),
                cols=colorRamp2(seq(-1.2,1.2,length.out=9), rev(colorpalette('rdbu',9))),
                cluster_cols=cluster_colors[c(1,3,4)],anno.size=12,out_f='Figure3B',width=6,height=10)

Figure3C(object=object,out_f='Figure3B_2',features=c("RFX4","TEAD3","FOS::JUN","TGIF2","NEUROG2(var.2)","NEUROD2"),
         n_cols=2,height=3.5,width=6)

res <- Figure3A(archr_obj=archr_obj,group.by='Cond',which_features = c('Astro','GFP','Ngn2','PmutNgn2'),
                tss.coords=read.table('data/mm10_TSS.bed'),
                features=c('Gfap','Aldoc','Fos','Tead3','Dcx','Rbfox3',"Hes5",
                           "Rbfox3:148916","Dab1:31912","Fos:-43134","Sox11:43842","Sox11:-14591","Sox11:-170662","Sox11:-495844","Sox11:-589737",
                           "Slc1a3:50970","Slc1a3:-72546","Slc1a3:-105947","Pcdh9:696185","Robo2:-115687","Cdkn1c:59587","Rnd2:17320","Auts2:23090","Gria1:203029","Slit1:9189"),
                cols=colorRamp2(seq(-1.2,1.2,length.out=length(heatmap_colors)), heatmap_colors),
                cluster_cols=sample_colors,anno.size=12,out_f='Figure3C',width=6,height=10)

res <- Figure3B(res=res,which_features = 'Distal',genome=BSgenome.Mmusculus.UCSC.mm10,
                fpkm_f='data/RNA/FPKM_AVGcounts_annot.txt',pwm_f='results/SC/combined_pwm.RDS',
                features=c("NEUROG2(var.2)","FOS::JUN","RFX4","IRF9","TGIF2","MEIS2","TEAD3","NEUROD2","TCF12(var.2)","JPD2","ZBTB18","TCF3","NFIB","HES6","NFIC(var.2)"),
                FDR.cutoff=2,enr.cutoff=0.5,mcparams=BiocParallel::MulticoreParam(10L),
                cols=colorRamp2(seq(-1.2,1.2,length.out=9), rev(colorpalette('rdbu',9))),
                cluster_cols=sample_colors,anno.size=12,out_f='Figure3D',width=6,height=10)

Figure3C(object=object,out_f='Figure3D_2',features=c("HES6","NFIB","TGIF2","TCF12(var.2)","TCF3"),
         n_cols=2,height=3.5,width=6)


mat <- Figure3D(mat=as.matrix(object.sub@assays$chromvar@data),
                nVar=100,object=object.sub,expr = uf_rpkm[rowMaxs(uf_rpkm,na.rm=T)>=3,],which_clusters = c('AST','iN_1','iN_2'),which_conditions = c('Astro','GFP','Ngn2','PmutNgn2'),
                features=c('TEAD2','FOS::JUNB','RFX4','SOX4','NEUROG2(var.2)','LHX2','YY1','TCF3','BHLHE22','TGIF2'),
                anno.size=12,clust_method='ward.D',out_f='Figure3_chromvar',
                cols=colorRamp2(seq(-8,8,length.out=11), rev(colorpalette('rdbu',11))),
                cluster_cols=cluster_cols[c(1,3,4)],width=8,height=10)

Figure3E(object=object,assay = 'chromvar',out_f='Figure3E',motifs=c('RFX4','NEUROG2(var.2)','TCF3'),cols2='matlablike2',point.size=1,width=5,height=5,plot_filled=T,stroke=0.1,theme = theme_border)


motif_footprinting(archr_obj=archr_obj,flank=500,motif='NEUROG2(var.2)',regions=NULL,
                   normalize_by='Subtract',out_f='Figure3F_1',cols=cluster_colors[c(1,3,4)],
                   group.by='Clusters',which_clusters=c('AST','iN_1','iN_2'),plot_bias=F,anno.size=12,key.size=4,width=6,height=4.5)
motif_footprinting(archr_obj=archr_obj,flank=500,motif='NEUROG2(var.2)',regions=NULL,
                   normalize_by='Subtract',out_f='Figure3F_2',cols=sample_colors,
                   group.by='Cond',which_clusters=c('Astro','GFP','Ngn2','PmutNgn2'),plot_bias=F,anno.size=12,key.size=4,width=6,height=4.5)


seRNA_all <- readRDS("results/SC/Cluster_PseudoBulk_RNA-Summarized-Experiment.RDS")
sePB_all <- readRDS("results/SC/Cluster_PseudoBulkATAC-Summarized-Experiment.RDS")
res <- Figure3G(selected_peaks=p2glinks$posCor,sortBy = 'RNA',scaleMax=F,seRNA_all =seRNA_all,sePB_all =sePB_all,out_f='Figure3G',
                  which_clusters = c('AST','iN_1','iN_2'),
                  features=c("Rbfox3:247336","Sox9:-280433","Sox9:234828","Gfap:-154240","Sox11:277561","Sox11:-233215","Sox4:-435166","Sox4:235153","Hes6:348009","Elavl2:132949","Sema5a:-344033","Sema5a:101208","Dll1:127766",
                             "Slc1a3:-105947","Slc1a3:38106","Cdkn1c:59587","Rnd2:17320","Auts2:23090","Auts2:409859","Auts2:-358036","Gria1:203029","Gria1:-455123","Slit1:131953","Slit1:46416"),
                  cols1=colorRamp2(seq(-2,2,length.out=length(heatmap_colors)), heatmap_colors),cluster_cols=cluster_colors[c(1,3,4)],
                  cols2=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))),anno.size=12,height=10,width=8)
res$labels <- paste0(res$gene_name,':',res$distance)
write.table(res,'results/P2G_ClusterBinaryMat_posCor.tsv',quote=F,col.names=T,row.names=F,sep='\t')

DefaultAssay(object) <- 'RNA'
o <- Figure3H(df=p2glinks$posCor,control_df=NULL,include_prom = T,
               file_f='Figure3H_1',target_TF='NEUROG2(var.2)',
               cols=colorRampPalette(c('blue','white','red'))(101),
               features=c('Cplx2','Sox11','Rbfox3','Rnd2','Igfbpl1','Kirrel3','Ralgds','Auts2','Cntn2','Slit1','Elavl2','Dll1','Onecut2','Zmiz1','Rybp','Tubb3','Plxna2'),GO=TRUE,
               chip_f=NULL,n_go=10,rna_fc=rna_fc,
               width=12,height=5,point.size=3,anno.size=5,leg.text.size=10,object=object)

o2 <- Figure3H(df=p2glinks$posCor,control_df=NULL,include_prom = T,
              file_f='Figure3H_1',target_TF='RFX4',
              cols=colorRampPalette(c('blue','white','red'))(101),
              features=c('Cplx2','Sox11','Rbfox3','Rnd2','Igfbpl1','Kirrel3','Ralgds','Auts2','Cntn2','Slit1','Elavl2','Dll1','Onecut2','Zmiz1','Rybp','Tubb3','Plxna2'),GO=TRUE,
              chip_f=NULL,n_go=15,rna_fc=rna_fc,
              width=14,height=5,point.size=3,anno.size=5,leg.text.size=10,object=object)

source('/home/hpc/bonev/projects/hic/dfg/config.R')

df <- p2glinks$posCor
motif.names <- object@assays$peaks@motifs@motif.names
motif.matrix <- object@assays$peaks@motifs@data
colnames(motif.matrix) <- as.vector(motif.names[match(names(motif.names),colnames(motif.matrix))])
motif.matrix <- motif.matrix[match(gsub('_','-',df$peakName),row.names(motif.matrix)),]
anno <- as.data.frame(motif.matrix[,'NEUROG2(var.2)',drop=F])
anno <- as.data.frame(stringr::str_split(row.names(anno)[anno$`NEUROG2(var.2)`>=1], pattern ='\\.' , n = 3, simplify = TRUE))
anno <- read.table('results/beds/allATAC_enh_Ngn2motif.bed')
anno <- gintervals(anno[,1],as.numeric(anno[,2]),as.numeric(anno[,3]))
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long


plotMisha(main_f=main_f,targetGene='Cplx2',outDir='plots/figures/',out_f='Figure3I',upstream=4.5e5,downstream=1e5,
          window_scale=2.2,pointCEX=0.5, binSize=5e3,conditions=score_tracks[2:4],
          chipNames=c('AST','iN_1','iN_2'),figure_mode=TRUE,
          chipTracksToExtract=c("scATAC.repro_AST","scATAC.repro_iN_1","scATAC.repro_iN_2"),
          arcIntervals=p2glinks$posCor,arcColors=colorRampPalette(ArchR::ArchRPalettes[[29]]),
          #chipYlim=matrix(c(0,7,0,7,0,7,0,0.8),ncol = 2,byrow = T),
          chipColors=c("#272E6A","#FA7517","#A30005","black"),img_factor=0.9,annIntervals=anno,
          plotOrder=list(scores=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE ,anno=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=FALSE),
          plotRatios=list(unitHeight=100, scores=2, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))




p <- VlnPlot(object, features = c('Cplx2'),cols=cluster_colors[c(1,3,4)],idents = c('AST','iN_1','iN_2'),pt.size = 0,combine = F)
p1 <- p[[1]]  + theme(legend.position='none',axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+ xlab('') + ylab('')+ scale_x_discrete(limits = rev) + coord_flip()+ ggtitle("")
pdf('plots/figures/Figure3I_2.pdf',height=4,width=3)
p1
dev.off()


plotMisha(main_f=main_f,targetGene='Clstn2',outDir='plots/figures/',out_f='Figure3M_PmutNgn2',upstream=4e5,downstream=0.9e5,
          window_scale=2.2,pointCEX=0.5, binSize=5e3,conditions=score_tracks[2:4],
          chipNames=c('Astro','GFP','Ngn2','PmutNgn2','Ngn2 ChIP','PmutNgn2 ChIP'),
          chipTracksToExtract=c("scATAC.repro_Astro","scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2",
                                "chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2"),
          chipYlim=matrix(c(0,3.5,0,3.5,0,3.5,0,3.5,0,0.7,0,0.7),ncol = 2,byrow = T),
          chipColors=c("#18C0C9","#215801","#E6CA17","#D51F26","black","black"),img_factor=1.2,
          plotOrder=list(scores=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE ,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=100, scores=2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))

##Supplementary Figures

res <- FigureS3A(archr_obj=archr_obj,group.by='Cond',which_features = c('PmutNgn2','Ngn2'),
                tss.coords=read.table('data/mm10_TSS.bed'),
                features=c("Rfx4:-43009","Fhl1:-49776","Grin2a:534762","Robo2:-187267","Cdkn1c:59587","Rnd2:17320","Auts2:23090","Gria1:203029","Slit1:9189"),
                sample_colors=cluster_colors[3:4],anno.size=12,out_f='FigureS3A_test',width=6,height=6)



res <- Figure3A(archr_obj=archr_obj,group.by='Cond',which_features = c('PmutNgn2','Ngn2'),
                tss.coords=read.table('data/mm10_TSS.bed'),
                features=c("Cdkn1c:59587","Rnd2:17320","Auts2:23090","Gria1:203029","Slit1:9189"),
                sample_colors=sample_colors,anno.size=12,out_f='Figure3D',width=6,height=6)

Figure3G(object=object,assay = 'chromvar',out_f='FigureS3D_1',motifs=c('TEAD3'),point.size=1,width=5,height=5,plot_filled=T,stroke=0.1,theme = theme_border)
Figure3G(object=object,assay = 'chromvar',out_f='FigureS3D_2',motifs=c('FOS::JUN'),point.size=1,width=5,height=5,plot_filled=T,stroke=0.1,theme = theme_border)
Figure3G(object=object,assay = 'chromvar',out_f='FigureS3D_3',motifs=c('TGIF2'),point.size=1,width=5,height=5,plot_filled=T,stroke=0.1,theme = theme_border)
Figure3G(object=object,assay = 'chromvar',out_f='FigureS3D_4',motifs=c('ZBTB18'),point.size=1,width=5,height=5,plot_filled=T,stroke=0.1,theme = theme_border)
Figure3G(object=object,assay = 'chromvar',out_f='FigureS3D_5',motifs=c('NEUROD2'),point.size=1,width=5,height=5,plot_filled=T,stroke=0.1,theme = theme_border)

motif_footprinting(archr_obj=archr_obj,flank=500,motif='TCF3',regions=NULL,
                   normalize_by='Subtract',out_f='FigureS3E_2',cols=cluster_colors[c(1,3,4)],
                   group.by='Clusters',which_clusters=c('AST','iN_1','iN_2'),plot_bias=F,anno.size=12,key.size=4,width=5,height=4.5)

motif_footprinting(archr_obj=archr_obj,flank=500,motif='TCF3',regions=NULL,
                   normalize_by='Subtract',out_f='FigureS3E',cols=sample_colors,
                   group.by='Cond',which_clusters=c('Astro','GFP','Ngn2','PmutNgn2'),plot_bias=F,anno.size=12,key.size=4,width=5,height=4.5)

FigureS3B(object,features='RFX4',out_f='FigureS3F_1',cols = sample_colors,width=5.5,height=4.5)
FigureS3B(object,features='NEUROG2(var.2)',out_f='FigureS3F_2',cols = sample_colors,width=5.5,height=4.5)
FigureS3B(object,features='TCF3',out_f='FigureS3F_3',cols = sample_colors,width=5.5,height=4.5)
FigureS3B(object,features='TEAD3',out_f='FigureS3F_4',cols = sample_colors,width=4,height=4)
FigureS3B(object,features='TGIF2',out_f='FigureS3F_5',cols = sample_colors,width=4,height=4)
