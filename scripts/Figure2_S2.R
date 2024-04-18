library(Seurat)
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
library(BSgenome.Mmusculus.UCSC.mm10)

## Set main folder as current dir ####

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

source('scripts/aux_functions.R')
source('scripts/figures/config.R')
source('scripts/figures/plot_functions.R')


object <- readRDS('results/SC/seurat_IDs.RDS')
cell_markers <- read.table('data/SC/cell_markers.tsv')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))

Figure2B <- function(object,cols,out_f='Figure2B',point.size,anno.size,key.size,height,width,plot_filled=F,theme=NULL,stroke=0.1,rows=1){
  plot.data <- extractFeatures(object,features='seurat_clusters',cells_toInclude = c('all'),cells_toExclude = 'none')
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=labels)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = labels), pch = I(21),size = point.size,stroke=stroke,alpha=1) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0,(1.03+0.025*rows)),plot.margin = margin(height/10,0.1,0.1,0.1,unit='inches'), legend.box = "horizontal",legend.spacing.x = unit(0.1, 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = rows)) 
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  if(!is.null(theme)){
    print(p + theme)}
  else {
    print(p)
  }
  dev.off()
}

Figure2C <- function(object,cols,out_f,point.size,anno.size,key.size,height,width,plot_filled=F,theme=NULL,stroke=0.1,rows=1){
  plot.data <- extractFeatures(object,features='orig.ident',cells_toInclude = c('all'),cells_toExclude = 'none')
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") + scale_colour_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke,alpha=1) + xlab("UMAP1") + ylab("UMAP2") + scale_fill_manual(name='',values=as.character(cols))  
    p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0,(1.03+0.025*rows)),plot.margin = margin(height/10,0.1,0.1,0.1,unit='inches'), legend.box = "horizontal",legend.spacing.x = unit(0.1, 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = rows)) 
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  if(!is.null(theme)){
    print(p + theme)}
  else {
    print(p)
  }
  dev.off()
}

Figure2D <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=reps,y=1,fill=labels)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='cluster',values = as.character(cols))
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') + theme(legend.position = 'none')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure2E <- function(object,features,cols,out_f,height,width){
  p <- DotPlot(object=object, features = features) + RotatedAxis()
  p <-p + scale_color_gradientn(colours = cols) + xlab('') + ylab('')
  p <- p + scale_y_discrete(limits=rev(levels(object))) + scale_x_discrete(limits=features)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure2F <- function(object,features,cols,point.size,height,width,plot_filled=F,anno.size=10,theme=NULL,direction='vertical',stroke=0.1,out_f){
  plot_list <- list()
  for (feature in features){
    plot.data <- extractFeatures(object,features=feature,cells_toInclude = c('all'),cells_toExclude = 'none')
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(color=cols[1],size = point.size,data = plot.data[plot.data$feature==0,]) + geom_point(aes(color=feature),size = point.size,data = plot.data[plot.data$feature>0,]) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_color_gradientn(colours=cols,name='',breaks=c(0,max(plot.data$feature)),labels=c(0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colours=cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(min(plot.data$feature,na.rm=T),round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1.05, 0.98),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 7.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))
      p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
    }
    if(!is.null(theme)){
      p <- p+theme
      p <- p + annotate("text",x = -Inf, y = Inf, label = feature,size=anno.size,hjust = -0.25, vjust = 1.5)
    } else {
      p <- p + ggtitle(feature)
    }
    plot_list[[feature]] <- p
  }
  p <- Reduce("+", plot_list) 
  if (direction=='vertical'){
    height=length(features)*height
    p_layout <- plot_layout(nrow=length(features),ncol=1)
  } else {
    width=length(features)*width
    p_layout <- plot_layout(ncol=length(features),nrow=1)
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p+p_layout)
  dev.off()
}

Figure2G <- function(object,pgraph=NULL,feature,cols,point.size,line.size,out_f,width,height,plot_filled=F,theme=NULL,stroke=0.1){
  plot.data <- extractFeatures(object,features=feature,cells_toInclude = c('all'),cells_toExclude = 'none')
  if(!plot_filled){
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2,color=feature)) + geom_point(size = point.size) + xlab("UMAP1") + ylab("UMAP2") 
    p <- p + scale_color_gradientn(colours = cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(round(min(plot.data$feature,na.rm=T)),round(max(plot.data$feature,na.rm=T))))
    p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.95, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.4), vjust = 8))
    p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
  } else {
    p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=stroke) + xlab("UMAP1") + ylab("UMAP2") 
    p <- p + scale_fill_gradientn(colours = cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(round(min(plot.data$feature,na.rm=T)),round(max(plot.data$feature,na.rm=T))))
    p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1.02, 0.98),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 7.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))
    p <- p + guides(fill = guide_colourbar(barwidth = 3, barheight = 1))
  }
  if (!is.null(pgraph)){
    edge_df <- readRDS(pgraph)
    p <- p + geom_segment(aes_string(x="source_prin_graph_dim_1",y="source_prin_graph_dim_2",xend="target_prin_graph_dim_1",yend="target_prin_graph_dim_2"),
                          size=line.size,color=I('black'),linetype="solid",na.rm=TRUE,data=edge_df)
  }
  if(!is.null(theme)){
    p <- p+theme
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure2H <- function(mat,object,features,cols,cluster_cols,cond_colors,anno.size,out_f,width,height){
  ra = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%(features)), labels = row.names(mat)[which(row.names(mat)%in%unique(c(features)))],labels_gp = gpar(fontsize = anno.size)))
  col.list <- cluster_cols
  names(col.list) <- levels(object)
  idents_indx <- Idents(object)
  idents_indx <- idents_indx[match(colnames(mat),names(idents_indx))]
  
  col.list2 <- cond_colors
  names(col.list2) <- unique(object$orig.ident)
  idents_indx2 <- object$orig.ident
  idents_indx2 <- idents_indx2[match(colnames(mat),names(idents_indx2))]
  
  ha = HeatmapAnnotation(cluster = idents_indx,cond=idents_indx2 ,col = list(cluster=col.list,cond=col.list2),show_legend = T,show_annotation_name=F,annotation_legend_param=list(cluster=list(direction='vertical'),cond=list(direction='vertical')))

  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  hm <- Heatmap(mat, name = "Pseudotime",cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, right_annotation = ra,col = cols,top_annotation = ha,
                heatmap_legend_param=list(labels = c(0,100),at=c(0,1), direction = "vertical",title='% Max',legend_height = unit(height/6, "inch")))
  draw(hm, merge_legends=T, ht_gap = unit(height/2, "inch")) 
  dev.off()
}

Figure2I <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data$feature[is.na(plot.data$feature)] <- 0
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = as.character(cols)) + geom_boxplot(width=0.1, fill="white",outlier.shape = NA) + theme(legend.position=c('none')) + xlab('') + ylab(features) + ggtitle('')
  return(p)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure2J <- function(object,features,out_f,cols,point.size,alpha,anno.size,key.size,height,width,plot_filled=F,direction='horizontal'){
  plot_list <- list()
  for (feature in features){
    plot.data <- extractFeatures(object,features=c(feature,'pseudotime'),cells_toInclude = c('all'),cells_toExclude = 'none')
    plot.data <- plot.data[!is.na(plot.data$pseudotime),]
    plot.data$labels <- droplevels(plot.data$labels)
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=pseudotime,y=feature,color=reps)) + geom_point(size = point.size,alpha=alpha) + ylab("Expression") + xlab("Pseudotime") + scale_colour_manual(name='',values=as.character(cols)) + geom_smooth(aes(x=pseudotime,y=feature),inherit.aes = F,colour='black')  
      p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) + ggtitle(feature)
    } else {
      p <- ggplot(plot.data, aes(x=pseudotime,y=feature)) + geom_point(aes(fill = labels),colour='black', pch = I(21),size = point.size,alpha=alpha) + ylab("Expression") + xlab("Pseudotime") + scale_fill_manual(name='',values=as.character(cols)) + geom_smooth(aes(x=pseudotime,y=feature),inherit.aes = F,colour='black')
      p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position="top", legend.box = "horizontal",legend.spacing.x = unit(width/(length(levels(object))*20), 'inch')) + guides(fill=guide_legend(override.aes = list(size = key.size),nrow = 1)) + ggtitle(feature)
    }
    plot_list[[feature]] <- p
  }
  if(direction=='vertical'){
    p <- Reduce("+", plot_list)
    p <- p + guide_area() + plot_layout(guides="collect",heights = c(rep(9,length(features)),1),ncol=1,nrow=length(features)+1)
    height=length(features)*height
  } else {
    p <- Reduce("|", plot_list) 
    p <- guide_area()/p + plot_layout(guides="collect",heights = c(1, 9))
    width=length(features)*width
  }
  
  
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + theme_cowplot())
  dev.off()
}

Figure2K <- function(object,clusters,out_f,height,width){
  plot_list <- list()
  for (cluster in clusters){
    marker_genes <- FindMarkers(object,ident.1 = cluster,ident.2 = NULL,only.pos = T,group.by = 'cluster',min.pct = 0.25,method='MAST')
    res <- enrichGO_wrapper(row.names(marker_genes[marker_genes$p_val_adj<=0.05,]))
    plot_list[[cluster]] <- res
  }
 # marker_genes <- FindMarkers(object,ident.1 = 'PmutNgn2',ident.2 = 'Ngn2',only.pos = T,group.by = 'orig.ident',min.pct = 0.25,method='MAST')
  #res <- enrichGO_wrapper(row.names(marker_genes[marker_genes$p_val_adj<=0.05,]))
  #plot_list[['PmutNgn2vsNgn2']] <- res
  res <- merge_result(plot_list)
  p <- dotplot(res, showCategory=11,font.size=12)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure2L <- function(res,features,point.size=1,anno.size=4,sample_colors,out_f,height,width,ylim=NULL,xlim=NULL){
  res <- res[!is.na(res$gene_symbol),]
  p <- ggplot(res,aes(x=log2FoldChange,y=-log10(padj))) + geom_point(color='grey',size = point.size) + geom_point(fill=sample_colors[3],alpha=0.8, pch = I(21),size = point.size,data = res[res$padj<0.05&res$log2FoldChange<0,]) + geom_point(fill=sample_colors[4],alpha=0.8, pch = I(21),size = point.size,data = res[res$padj<0.05&res$log2FoldChange>0,])  
  p <- p + xlab(expression(Log[2]~Fold~Change)) + ylab(expression(-Log[10]~(P)))
  p <- p + ggrepel::geom_text_repel(
    data = res, size = anno.size,seed = 42,
    box.padding =0.8, min.segment.length = 0,max.iter = 10000,max.overlaps = Inf,
    aes(x=log2FoldChange,y=-log10(padj),color=NULL,label=ifelse(gene_symbol%in%features, as.character(gene_symbol), "")),force=10)
  if(!is.null(ylim)){
    p <- p+ ylim(ylim)
  }
  if(!is.null(xlim)){
    p <- p+ xlim(xlim)
  }
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}  
  

Figure2M <- function(res,fc,direction,orgDB=org.Mm.eg.db,plot_what='GO',qvalue.cutoff=0.05,sample_colors,out_f,height,width){
  require(clusterProfiler)
  require(enrichplot)
  df <- res[res$padj<=0.05,]
  if(direction=='up'){
    df <- df[df$log2FoldChange>=fc,]
  } else if(direction=='down'){
    df <- df[df$log2FoldChange<=fc,]
  } else {
    df <- df[abs(df$log2FoldChange)>=fc,]
  }
  geneList <- as.numeric(df$log2FoldChange)
  names(geneList) <- as.character(df$gene_symbol)
  geneList <- geneList[order(geneList,decreasing = T)]
  names(geneList) <- bitr(names(geneList), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = orgDB)$ENTREZID
  if(plot_what=='GO'){
    ego <- enrichGO(gene         = unique(names(geneList)),
                    OrgDb         = orgDB,
                    keyType       = 'ENTREZID',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = qvalue.cutoff)
    #p <- simplify(ego)
  } else if(plot_what=='GSEA'){
    ego <- gseGO(geneList     = geneList,
                  OrgDb        = orgDB,
                  ont          = "BP",
                  keyType       = 'ENTREZID',
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH")
    #p <- simplify(ego)
  } else if(plot_what=='KEGG'){
   ego <- enrichKEGG(gene        = names(geneList),
                   organism     = 'mmu',
                   pvalueCutoff = 0.05)
  }
  return(list(ego=ego,geneList=geneList))
}

#### Supplementary Figure 2

FigureS2A <- function(object_f,point.size,sample_colors,min_nCount_ATAC,max_nCount_ATAC,out_f,min_nCount_RNA,max_nCount_RNA,height,width){
  object <- readRDS(object_f)
  mat <- data.frame(UniqueFragments=object$nCount_ATAC,RNA_UMIs=object$nCount_RNA,TSS_enrichment=object$TSS.enrichment,NuSignal=object$nucleosome_signal,ID=object$orig.ident)
  p <-  ggplot(mat, aes( x = UniqueFragments, y = RNA_UMIs,color=ID )) + geom_point(size = point.size) + scale_color_manual(name='',values = sample_colors)             # scale_color_gradientn(colours=tim.colors(12))
  p<- p + xlab('ATAC Fragments') + ylab('RNA UMIs') + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) 
  p <- p + geom_vline(xintercept=c(min_nCount_ATAC,max_nCount_ATAC),linetype = 2) + geom_hline(yintercept=c(min_nCount_RNA,max_nCount_RNA),linetype = 2) 
  p <-  p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.75, 1))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS2B <- function(object_f,point.size,sample_colors,min_TSS,max_TSS,out_f,min_NuSignal,max_NuSignal,height,width){
  object <- readRDS(object_f)
  mat <- data.frame(UniqueFragments=object$nCount_ATAC,RNA_UMIs=object$nCount_RNA,TSS_enrichment=object$TSS.enrichment,NuSignal=object$nucleosome_signal,ID=object$orig.ident)
  p <-  ggplot(mat, aes( x = TSS_enrichment, y = NuSignal,color=ID )) + geom_point(size = point.size) + scale_color_manual(name='',values = sample_colors)             # scale_color_gradientn(colours=tim.colors(12))
  p<- p + xlab('TSS Enrichment') + ylab('Nucleosome Signal') 
  p <- p + geom_vline(xintercept=c(min_TSS,max_TSS),linetype = 2) + geom_hline(yintercept=c(min_NuSignal,max_NuSignal),linetype = 2) 
  p <-  p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.75, 1))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS2C <- function(archr_obj_path,logFile,cols,out_f,height,width){
  require(ArchR)
  require(doParallel)
  addArchRThreads(threads = 1)
  archr_obj <- loadArchRProject(path=archr_obj_path,force=T)
  df <- plotFragmentSizes(ArchRProj = archr_obj,maxSize = 1000,returnDF = TRUE,groupBy = 'Cond',logFile = logFile )
  plotDF <- data.frame(df)
  p <- ggplot(plotDF, aes(fragmentSize, fragmentPercent,color=group)) + 
    geom_line(size = 0.75) +
    xlab("Fragment Size (bp)") +
    ylab("Percentage of Fragments") + theme(legend.position=c(0.8, 0.95)) + 
    scale_color_manual(values = cols,name='') +
    scale_y_continuous(limits = c(0, max(plotDF$fragmentPercent)*1.05), expand = c(0,0)) +
    scale_x_continuous(limits = c(min(plotDF$fragmentSize), 500),expand = expansion(mult = c(0,0.05)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS2D <- function(object,cols,out_f,height,width){
  p <- TSSPlot(object, assay = 'ATAC', group.by = 'orig.ident', idents = NULL)
  plotDF <- p$data 
  p <- ggplot(plotDF, aes(position, norm.value,color=group)) + 
    geom_line(size = 0.75) +
    xlab("Distance From TSS (bp)") +
    ylab("Mean TSS enrichment") + theme(legend.position=c(0.8, 0.95)) + 
    scale_color_manual(values = cols,name='') +
    scale_y_continuous(limits = c(0, max(plotDF$norm.value)*1.05), expand = c(0,0)) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS2G <- function(object,cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),sample_colors,fpkm_f,out_f,height,width){
  require(edgeR)
  require(LSD)
  object_cpm <- as.data.frame(matrix(NA,nrow=nrow(object),ncol=length(unique(object$orig.ident))))
  colnames(object_cpm) <- unique(object$orig.ident)
  for (s in unique(object$orig.ident)){
    object_sub <- object[,object$orig.ident==s]
    object_cpm[,s] <- cpm(rowSums(object_sub@assays$RNA@counts),log=T,prior.count = 1)
  }
  colnames(object_cpm) <- paste0('pseudo_',colnames(object_cpm))
  row.names(object_cpm) <- row.names(object)
  fpkm_df <- read.table(fpkm_f,header=T,sep='\t')
  fpkm_df[,1:9] <- 2^fpkm_df[,1:9]
  for (s in 1:9){
    fpkm_df[,s] <- cpm(fpkm_df[,s],log=T,prior.count = 1)
  }
  df <- merge(object_cpm,fpkm_df,by.x='row.names',by.y='geneID',sort=F)
  if(!is.null(filter_genes)){
    df <- df[df$Row.names%in%filter_genes,]
  }
  res <- cor(df[,2:14],method = 'pearson')
  distance <- dist(1-res)
  res <- round(res, 2)
  annotation_df <- data.frame(method=factor(c(rep('single_cell',4),rep('bulk',9)),levels=c('single_cell','bulk')),
                              condition=factor(c('Astro','GFP','Ngn2','PmutNgn2',rep('GFP',3),rep('Ngn2',3),rep('PmutNgn2',3)),levels=c('Astro','GFP','Ngn2','PmutNgn2')))
  row.names(annotation_df) <- colnames(res)
  names(sample_colors) <- unique(object$orig.ident)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  hm <- ComplexHeatmap::pheatmap(res,color=cols,breaks=seq(0.5,1,length=101),clustering_method='ward.D2',show_colnames = F,
           annotation_colors=list(method=c(single_cell='darkblue',bulk='darkred'),condition=sample_colors),
           display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.50,0.75,1))
  draw(hm,merge_legend = TRUE)
  dev.off()
}

FigureS2G_alt <- function(df1,df2,cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),sample_colors,fpkm_f,out_f,height,width){
  df <- merge(df1,df2,by='row.names',sort=F)
  if(!is.null(filter_genes)){
    df <- df[df$Row.names%in%filter_genes,]
  }
  df <- df[,c(2,4,5,11,13,15,16,17)]
  res <- cor(df,method = 'pearson')
  distance <- dist(1-res)
  res <- round(res, 2)
  res <- res[1:3,]
  annotation_df <- data.frame(system=factor(c(rep('iN',3),rep('E14_cortex',5)),levels=c('iN','E14_cortex')))
  row.names(annotation_df) <- colnames(res)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  hm <- ComplexHeatmap::pheatmap(res,color=cols,breaks=seq(0.5,1,length=101),clustering_method='ward.D2',show_colnames = T,cluster_rows=F,cluster_cols=F,
                                 annotation_colors=list(system=c(iN="#333333",E14_cortex="#CCCCCC")),annotation_legend = F,legend = F,angle_col=c("0"),
                                 display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.50,0.75,1))
  draw(hm,merge_legend = TRUE)
  dev.off()
}

FigureS2H <- function(misha_root,cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),tracks,track_names,bins,blacklist='/home/hpc/bonev/annotations/mm10/mm10_blacklist.bed'){
  require(misha)
  gsetroot(misha_root)
  df <- misha_extract(tracks = tracks,track_names = track_names,regions=gintervals.all(),iterator=bins,blacklist = blacklist)
  df <- df[rowSums(df[,4:(ncol(df)-1)],na.rm=T)>=0.5,]
  res <- cor(df[,4:(ncol(df)-1)],method = 'pearson',use='complete.obs')
  distance <- dist(1-res)
  res <- round(res, 2)
  annotation_df <- data.frame(method=factor(c(rep('single_cell',4),rep('bulk',9)),levels=c('single_cell','bulk')),
                              condition=factor(c('Astro','GFP','Ngn2','PmutNgn2',rep('GFP',3),rep('Ngn2',3),rep('PmutNgn2',3)),levels=c('Astro','GFP','Ngn2','PmutNgn2')))
  row.names(annotation_df) <- colnames(res)
  names(sample_colors) <- unique(object$orig.ident)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  hm <- ComplexHeatmap::pheatmap(res,color=cols,breaks=seq(0.7,1,length=101),clustering_method='ward.D2',show_colnames = F,
           annotation_colors=list(method=c(single_cell='darkblue',bulk='darkred'),condition=sample_colors),
           display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE ,legend_breaks=c(0.70,0.85,1))
  draw(hm,merge_legend = TRUE)
  dev.off()
}





###########################
#### Plot Figures #########
###########################

### Figure 2 

Figure2B(object,cols=cluster_colors,out_f='Figure2B',point.size=2,anno.size=12,key.size=3,height=6,width=5.5,plot_filled=F,theme = theme_border,rows=1,stroke=0.2)
Figure2C(object,cols=sample_colors,out_f='Figure2C',point.size=1.5,anno.size=12,key.size=3,height=6,width=5.5,plot_filled=F,theme = theme_border,rows=1,stroke=0.2)
Figure2D(object,features='seurat_clusters',cols=cluster_colors,out_f='Figure2D',height=6,width=4)
DefaultAssay(object) <- 'RNA'
Figure2E(object,features=as.vector(cell_markers[,1]),cols=c('grey','red','black'),out_f='Figure2E',height=5,width=8)
Figure2F(object=object,c('Gfap','Sox11','Dcx'),out_f='Figure2F',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
Figure2G(object,pgraph='results/SC/Monocle3_int/ptime_edge_df.RDS',feature='pseudotime',out_f='Figure2G',cols=brewer.pal(9,'Purples'),point.size=2,width=6,height=6,plot_filled=T,theme=theme_border,stroke=0.2,line.size=1)
Figure2H(mat=as.matrix(readRDS('results/SC/Monocle3_int/pd_genes_GAMmat.RDS')),object=object,
         features=cell_markers[,1],
         anno.size=12,out_f='Figure2H',
         cols=rev(brewer.pal(n = 9, name = "RdYlBu")),
         cluster_cols=cluster_colors,cond_colors=sample_colors,width=8,height=10)
Figure2I(object,features='pseudotime',out_f='Figure2I',cols = sample_colors,width=5,height=4)
Figure2J(object,out_f='Figure2J',features=c('Gfap','Sox11','Dcx'),cols=sample_colors,point.size=1,anno.size=12,key.size=4,alpha=0.5,width=6,height=4,plot_filled=F,direction='vertical')
Figure2L(res = read.table('data/RNA/res.PmutNgn2-Ngn2.txt'),ylim=c(0,40),xlim=c(-6,6),
         features <- c('Map2','Brsk2','Reln','Gfap','Yy1','Wnt7b','Clstn2','Bhlhe22','Sox9','Aldoc','Slc1a3'),
         point.size=2,sample_colors = sample_colors,out_f = 'Figure2L',height=5,width=5)
ego <- Figure2M(res = read.table('data/RNA/res.PmutNgn2-Ngn2.txt'),fc = 0.5,plot_what='GO',
              direction = 'no',orgDB = org.Mm.eg.db,qvalue.cutoff = 0.05,sample_colors = sample_colors)
p <- dotplot(simplify(ego$ego), showCategory=10,font.size=12)
pdf('plots/figures/Figure2M.pdf',height=6,width=6)
print(p+ theme(legend.position = c(0.8,0.35)) + ggtitle('PmutNgn2 vs Ngn2'))
dev.off()

#Figure2L
  p1 <- Figure2I(object,features='Cdkn1c',out_f='Figure2L_1',cols = sample_colors,width=4,height=4)
  p2 <- Figure2I(object,features='Reln',out_f='Figure2L_2',cols = sample_colors,width=4,height=4)
  p3 <- Figure2I(object,features='Clstn2',out_f='Figure2L_3',cols = sample_colors,width=4,height=4)
  p <- p1+p2+p3
  p <- Reduce("+", plot_list)
  p <- p + guide_area() + plot_layout(guides="collect",heights = c(rep(9,3),1),ncol=1,nrow=3)
  pdf('plots/figures/Figure2L.pdf',height=12,width=4)
  p
  dev.off()

#Supplementary Figure 2 ##

FigureS2A(object_f='results/SC/preFilter_Seurat.RDS',point.size=0.8,sample_colors=sample_colors,
          min_nCount_ATAC=8000,max_nCount_ATAC=125000,min_nCount_RNA=1000,max_nCount_RNA=30000,
          out_f='FigureS2A',height=6,width=6)
FigureS2B(object_f='results/SC/preFilter_Seurat.RDS',point.size=0.8,sample_colors=sample_colors,
          min_TSS = 1,max_TSS=20,min_NuSignal = 0.2,max_NuSignal = 2,
          out_f='FigureS2B',height=6,width=6)  
FigureS2C(archr_obj_path='data/SC/ArchR/seurat_atac/',logFile='data/SC/ArchR/ArchRLogs/FragmentSize.log',cols=sample_colors,out_f='FigureS2C',height=6,width=7) 
FigureS2D(object,cols=sample_colors,out_f='FigureS2D',height=5,width=6) 
DefaultAssay(object) <- 'GeneActivity'
Figure2E(object,features=as.vector(cell_markers[,1]),cols=c('grey','red','black'),out_f='FigureS2E',height=5,width=8)
DefaultAssay(object) <- 'RNA'
Figure2F(object=object,c('Sox9','Rnd2'),out_f='FigureS2F_1',cols=gene_colors(50),point.size=2,width=6,height=6,plot_filled=T,direction='vertical',anno.size=8,theme=theme_border,stroke=0.2)
DefaultAssay(object) <- 'GeneActivity'
Figure2F(object=object,c('Sox9','Rnd2'),out_f='FigureS2F_2',cols=viridis(50),point.size=2,width=6,height=6,plot_filled=T,direction='vertical',anno.size=8,theme=theme_border,stroke=0.2)
FigureS2G(object,cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),sample_colors=sample_colors,fpkm_f=fpkm_f='data/RNA/log2_norm_counts_annot.txt',out_f='FigureS2G',height=11,width=9)
FigureS2H(misha_root,cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),
          tracks <- c("scATAC.repro_Astro","scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2",
                      "atac.repro_GFP_rep1","atac.repro_GFP_rep2","atac.repro_GFP_rep3",
                      "atac.repro_Ngn2_rep1","atac.repro_Ngn2_rep2","atac.repro_Ngn2_rep3",
                      "atac.repro_PmutNgn2_rep1","atac.repro_PmutNgn2_rep2","atac.repro_PmutNgn2_rep3"),
          track_names <- c('sc_Astro','sc_GFP','sc_Ngn2','sc_PmutNgn2',
                           'GFP_rep1','GFP_rep2','GFP_rep3',
                           'Ngn2_rep1','Ngn2_rep2','Ngn2_rep3',
                           'PmutNgn2_rep1','PmutNgn2_rep2','PmutNgn2_rep3'),
          bins=1e4,blacklist='/home/hpc/bonev/annotations/mm10/mm10_blacklist.bed')
  
source('/home/hpc/bonev/projects/hic/dfg/config.R')

plotMisha(main_f=main_f,targetGene='Rnd2',outDir='plots/figures/',out_f='FigureS2I',upstream=2e4,downstream=7e4,
          window_scale=2.2,pointCEX=0.5, binSize=5e3,conditions=score_tracks[2:4],
          chipNames=c('Astro','GFP','Ngn2','PmutNgn2','GFP','Ngn2','PmutNgn2','GFP','Ngn2','PmutNgn2'),
          chipTracksToExtract=c("scATAC.repro_Astro","scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2",
                                "atac.repro_GFP","atac.repro_Ngn2","atac.repro_PmutNgn2",
                                "rnaseq_RPM.repro_GFP","rnaseq_RPM.repro_Ngn2","rnaseq_RPM.repro_PmutNgn2"),
          chipYlim=matrix(c(0,10,0,10,0,10,0,10,0,4.5,0,4.5,0,4.5,0,30,0,30,0,30),ncol = 2,byrow = T),
          chipColors=c("#18C0C9","#215801","#E6CA17","#D51F26","#215801","#E6CA17","#D51F26","#215801","#E6CA17","#D51F26"),img_factor=1.5,
          plotOrder=list(scores=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE ,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=100, scores=2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))