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
library(seqplots)

## Set main folder as current dir ####

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

source('scripts/figures/config.R')
source('scripts/aux_functions.R')
source('scripts/figures/plot_functions.R')

source('/home/hpc/bonev/projects/hic/dfg/config.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/main_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/aux_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/plot_functions.R')


object <- readRDS('results/SC/Yy1_SCT_Cell_type.RDS')
cell_markers <- read.table('data/SC/cell_markers.tsv')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))

meth_tracks <- c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x")
chip_tracks <- c("chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2")
atac_tracks <- c("scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2")


Figure7F_new <- function(score_f1,score_f2,min_dist,max_dist,labels,cols,p_stats,add_pvalue=T,type_labels=c('NotBound','Bound'),tip.lengths = c(0.03,0.01)){
  df1 <- get(load(score_f1))
  df1 <- as.data.frame(cbind(df1[[1]]$intra[,1:6],df1[[1]]$intra$v_score,df1[[2]]$intra$v_score))
  colnames(df1)[7:ncol(df1)] <- labels
  df1$type <- type_labels[1]
  df2 <- get(load(score_f2))
  df2 <- as.data.frame(cbind(df2[[1]]$intra[,1:6],df2[[1]]$intra$v_score,df2[[2]]$intra$v_score))
  colnames(df2)[7:ncol(df2)] <- labels
  df2$type <- type_labels[2]
  df <- rbind(df1,df2)
  df$dist <- df$start2-df$start1
  df <- df[(df$dist>=min_dist)&(df$dist<=max_dist),-ncol(df)]
 # df_e <- df_e[complete.cases(df_e[,1:2]),]
  df_o <- melt(df[,-c(1:6)],id.vars='type')
  df_o$variable <- factor(df_o$variable,levels=labels)
  df_o$type <- factor(df_o$type,levels=type_labels)
  p <- ggplot(df_o,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
  p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('Hi-C score') + theme(legend.position = "none")  + ggtitle('')
  if(add_pvalue){
    p <- p + stat_compare_means(comparisons = p_stats,label = "p.format",method = 'wilcox',paired = T,tip.length = tip.lengths) 
  }
  return(p+facet_wrap(~type))
}


Figure7E <- function(object,cols,out_f,point.size,anno.size,key.size,height,width,plot_filled=F,theme=NULL,stroke=0.1,rows=1){
  plot.data <- extractFeatures(object,features='orig.ident',cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data <- plot.data[sample(row.names(plot.data),nrow(plot.data)),]
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

Figure7F <- function(object,cols,out_f,point.size,anno.size,key.size,height,width,plot_filled=F,theme=NULL,stroke=0.1,rows=1){
  plot.data <- extractFeatures(object,features='seurat_clusters',cells_toInclude = c('all'),cells_toExclude = 'none')
  plot.data <- plot.data[sample(row.names(plot.data),nrow(plot.data)),]
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



Figure7G <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=reps,y=1,fill=labels)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='cluster',values = as.character(cols))
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') + theme(legend.position = 'none')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure7H <- function(object,features,cols,out_f,height,width){
  p <- DotPlot(object=object, features = features) + RotatedAxis()
  p <-p + scale_color_gradientn(colours = cols) + xlab('') + ylab('')
  p <- p + scale_y_discrete(limits=rev(levels(object))) + scale_x_discrete(limits=features)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

Figure7I <- function(object,features,cols,point.size,height,width,plot_filled=F,anno.size=10,theme=NULL,direction='vertical',stroke=0.1,out_f){
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

Figure7J <- function(object,clusters,group.by = 'cluster',showGO=11,out_f,height,width){
  plot_list <- list()
  for (cluster in clusters){
    marker_genes <- FindMarkers(object,ident.1 = cluster,ident.2 = NULL,only.pos = T,group.by = group.by,min.pct = 0.25,method='MAST')
    res <- enrichGO_wrapper(row.names(marker_genes[marker_genes$p_val_adj<=0.05,]))
    plot_list[[cluster]] <- res
  }
  # marker_genes <- FindMarkers(object,ident.1 = 'PmutNgn2',ident.2 = 'Ngn2',only.pos = T,group.by = 'orig.ident',min.pct = 0.25,method='MAST')
  #res <- enrichGO_wrapper(row.names(marker_genes[marker_genes$p_val_adj<=0.05,]))
  #plot_list[['PmutNgn2vsNgn2']] <- res
  res <- merge_result(plot_list)
  p <- dotplot(res, showCategory=showGO,font.size=12)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + coord_flip())
  dev.off()
  return(p)
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

###Supp Figure 7
FigureS7A <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=reps,fill=labels)) + geom_bar(stat="count",fill=cols,color='black') + scale_fill_manual(values = as.character(cols))
  p <- p + ylab('Number of Cells Passing Filter') + xlab('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p + ggtitle(paste0('Total = ',nrow(plot.data))))
  dev.off()
}

FigureS7B <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = cols) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of UMIs / cell') + ggtitle(paste0('Median = ',median(plot.data$feature)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS7C <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=labels,y=feature,fill=labels)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = as.character(cols)) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of UMIs / cell') + ggtitle('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS7D <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=reps,y=feature,fill=reps)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = cols) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of Genes / cell') + ggtitle(paste0('Median = ',median(plot.data$feature)))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS7E <- function(object,features,cols,out_f,height,width){
  plot.data <- extractFeatures(object,features=features,cells_toInclude = c('all'),cells_toExclude = 'none')
  p <- ggplot(plot.data, aes(x=labels,y=feature,fill=labels)) + geom_violin(trim=T,scale = "width") + scale_fill_manual(values = as.character(cols)) + geom_boxplot(width=0.1, fill="white") + theme(legend.position=c('none')) + xlab('') + ylab('Number of Genes / cell') + ggtitle('')
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}


FigureS7F <- function(object,group.by,which_features,features,fdr.cutoff=0.05,avg_log2FC.cutoff=0.25,point.size=1,cols,anno.size,out_f,width,height,xlim=NULL,ylim=NULL){
  mat <- FindMarkers(object,ident.1 = which_features[1],ident.2 = which_features[2],group.by = group.by,test.use = 'MAST',logfc.threshold = 0)
  p <- ggplot(mat,aes(x=avg_log2FC,y=-log10(p_val_adj))) + geom_point(color='grey',size = point.size) + geom_point(fill=cols[1],alpha=0.8, pch = I(21),size = point.size,data = mat[mat$p_val_adj<fdr.cutoff&mat$avg_log2FC<(avg_log2FC.cutoff*(-1)),]) + geom_point(fill=cols[2],alpha=0.8, pch = I(21),size = point.size,data = mat[mat$p_val_adj<fdr.cutoff&mat$avg_log2FC>avg_log2FC.cutoff,])  
  p <- p + xlab(expression(Log[2]~Fold~Change)) + ylab(expression(-Log[10]~(P)))
  p <- p + ggrepel::geom_text_repel(
    data = mat, size = anno.size,seed = 42,
    box.padding =0.8, min.segment.length = 0,max.iter = 10000,max.overlaps = Inf,
    aes(x=avg_log2FC,y=-log10(p_val_adj),color=NULL,label=ifelse(labels%in%features, as.character(labels), "")),force=10)
  p <- p + ggtitle(paste0(which_features[1],' vs ',which_features[2]),subtitle=paste0(' Up: ',nrow(mat[mat$p_val_adj<fdr.cutoff&mat$avg_log2FC>avg_log2FC.cutoff,]),' Down:',nrow(mat[mat$p_val_adj<fdr.cutoff&mat$avg_log2FC<(avg_log2FC.cutoff*(-1)),])))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS7G <- function(object,which_features,features,point.size=1,scale=T,cols=yy1_colors,cols2=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))),anno.size=12,out_f,width,height){
  df <- FindMarkers(object,ident.1 = which_features,ident.2 = NULL,only.pos = T)
  
  object_cpm <- as.data.frame(matrix(NA,nrow=nrow(object@assays$RNA@counts),ncol=length(unique(object$orig.ident))))
  colnames(object_cpm) <- unique(object$orig.ident)
  for (s in unique(object$orig.ident)){
    object_sub <- object[,object$orig.ident==s]
    object_cpm[,s] <- cpm(rowSums(object_sub@assays$RNA@counts),log=T,prior.count = 1)
  }
  row.names(object_cpm) <- row.names(object@assays$RNA@counts)
  mat <- as.matrix(object_cpm[row.names(object_cpm)%in%row.names(df),])
  mat_m <- melt(mat,measure.vars = levels(object$orig.ident))
  
  if (scale) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
  }
  col.list1 <- cols
  names(col.list1) <- colnames(mat) 
  la2 <- columnAnnotation(foo = anno_text(colnames(mat), location = 0.5,rot=0, just = "center",gp = gpar(fill = col.list1, col = "black", border = "black", fontsize = 12),height = unit(height/50,'inch')))
  if(!is.null(features)){
    ra1 = rowAnnotation(foo = anno_mark(at = which(row.names(mat)%in%features),side='right', labels = row.names(mat)[which(row.names(mat)%in%features)],labels_gp = gpar(fontsize = anno.size)))
  } else {
    ra1 = NULL
  }
  mat_k <- tglkmeans::TGL_kmeans(mat,k = 5,id_column = F,reorder_func = 'median',hclust_intra_clusters = F,seed = 42)
  hm <- Heatmap(mat,name='Normalized Expression',column_title = 'Expression',row_split = mat_k$cluster,show_row_dend = F,cluster_columns = F,show_column_names = F,show_row_names = F,use_raster = TRUE, raster_quality = 10,col = cols2,top_annotation = la2,right_annotation = ra1,
                 heatmap_legend_param=list(direction = "horizontal",title='Expression Z-score'))
  pdf(paste0('plots/figures/',out_f,'_1.pdf'),height=height,width=width)
  draw(hm,heatmap_legend_side='bottom',annotation_legend_side = "bottom")
  dev.off()
  p <- ggplot(mat_m,aes(x=Var2,y=value,fill=Var2)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8) 
  p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('Normalized Expression') + theme(legend.position = "none")
  pdf(paste0('plots/figures/',out_f,'_2.pdf'),height=8,width=6)
  print(p)
  dev.off()
}


###########################
#### Plot Figures #########
###########################

### Figure 7 

tracks=c('CutRun/full/wt_GFP_YY1/rep1/wt_GFP_YY1_rep1/peakcalling/macs2.narrow/wt_GFP_YY1_rep1.cpm.norm.bw',
         'CutRun/full/wt_ngn2_YY1/rep1/peakcalling/macs2.broad/wt_ngn2_YY1_rep1.cpm.norm.bw',
         'CutRun/full/wt_pmut_ngn2_YY1/rep1/peakcalling/macs2.narrow/wt_pmut_ngn2_YY1_rep1.cpm.norm.bw')
peaks='results/ChIP/Ngn2vsGFP_Yy1_LFC2.bed'
peaks_anno <- read.table('CutRun/full/Ngn2vsGFP_Yy1_LFC2.tsv',header=T)

hm_l <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=500,bin=10,clusters=read.table(peaks)$V5,cols=rev(colorpalette('reds')),zlims = matrix(c(0,1,0,1,0,0.6),byrow = T,ncol=2),raster_quality = 5,show_leg = F)
pdf('plots/figures/Figure7B_new.pdf',width=5,height=8)
draw(hm_l$hm_l[[1]] + hm_l$hm_l[[2]] + hm_l$hm_l[[3]],merge_legends=F)
dev.off()

peaks1 <- read.table('results/ChIP/Ngn2vsGFP_Yy1_LFC2.bed',header=F)[,c(1:3,5)]
peaks2 <- read.table('data/ChIP/conservative_peaks.narrowPeak',header=T)[,c(1:3)]
colnames(peaks1) <- c('chrom','start','end','Yy1_direction')
colnames(peaks2) <- c('chrom','start','end')
peaks2 <- peaks2[!duplicated(peaks2),]
df <- gintervals.neighbors(peaks1,peaks2)
df$Ngn2_direction[df$dist<=500] <- 'Yes'
df$Ngn2_direction[df$dist>500] <- 'No'
df$Ngn2_direction <- factor(df$Ngn2_direction,levels=c('No','Yes'))
df$Yy1_direction <- factor(df$Yy1_direction)
levels(df$Yy1_direction)=c('GFP only','shared','Ngn2 only')
df$Yy1_direction <- factor(df$Yy1_direction,levels=rev(levels(df$Yy1_direction)))
df <- df[,c('Yy1_direction','Ngn2_direction')]
p <- ggplot(df, aes(x=Yy1_direction,y=1,fill=Ngn2_direction)) +  geom_bar(position="fill", stat="identity") + scale_fill_grey(name='',start = 0.8,end=0.2) 
p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')
p <- p + theme(legend.position=c(0.65,(1.05)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 0.5),nrow = 1))
pdf('plots/figures/Figure7C_new.pdf',width=5,height=4,useDingbats = F)
print(p+coord_flip())
dev.off()

res <- read.table('CutRun/full/Ngn2vsGFP_Yy1_LFC2.tsv',header=T)
res_hm <- Figure4C(peaks=res[,1:3],logFC=res$logFC,bins=factor(res$direction,levels=c('GFP','shared','Ngn2')),central_bin=2,breaks=NULL,peak_size=200,
                   pwms=readRDS('data/combined_pwm_FPKM1.RDS'),genome_bg=T,genome=BSgenome.Mmusculus.UCSC.mm10,
                   features=c('NEUROG2(var.2)','BHLHE22(var.2)','TGIF2','E2F4','NFIX(var.2)','FOS','FOS::JUN','RFX4','TEAD2','TEAD3'),
                   mcparams=BiocParallel::MulticoreParam(10L),FDR.cutoff=4,enr.cutoff=1,cols=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))))

pdf('plots/figures/Figure7D_new.pdf',width=4,height=8)
draw(res_hm$hm,merge_legends=F,heatmap_legend_side = "bottom")
dev.off()

library(ChIPpeakAnno)

peaks1 <- read.table('results/ChIP/Ngn2vsGFP_Yy1_LFC2.bed')
peaks1 <- peaks1[,c(1:3,5)]
colnames(peaks1) <- c('chrom','start','end','type')
df <- gintervals.neighbors(peaks1,gintervals.load(tss_f))
df$g_type <- 'Distal'
df$g_type[abs(df$dist)<=5000] <- 'Promoter'
df$type <- as.factor(df$type)
levels(df$type) <- c('GFP only','shared','Ngn2 only')
df$type <- factor(df$type,levels=rev(levels(df$type)))

df <- df[,c('type','g_type')]

p <- ggplot(df,aes(x=type,y=1,fill=g_type)) + geom_bar(position="fill", stat="identity") + scale_fill_grey(end = 0.2,start=0.9)
p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') 
p <- p + theme(legend.text = element_text(size=12)) + theme(legend.position=c(0.25,(1.05)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
pdf('plots/figures/Figure7E_new.pdf',width=5,height=4,useDingbats = F)
print(p + coord_flip())
dev.off()

peaks1 <- read.table('CutRun/full/Ngn2vsGFP_Yy1_LFC2.tsv',header=T)[,c(1:3,11)]
peaks2 <- read.table('data/ChIP/Ngn2_peaks5000_distal.bed',header=T)[,c(1:3)]
colnames(peaks1) <- c('chrom','start','end','Yy1_direction')
colnames(peaks2) <- c('chrom','start','end')
peaks2 <- peaks2[!duplicated(peaks2),]
df <- gintervals.neighbors(peaks2,peaks1)
df$Ngn2_direction[df$dist<=100] <- 'Bound'
df$Ngn2_direction[df$dist>100] <- 'NotBound'
df$Ngn2_direction <- factor(df$Ngn2_direction,levels=c('NotBound','Bound'))
df$Yy1_direction <- factor(df$Yy1_direction,levels=c('GFP','shared','Ngn2','PmutNgn2','NotBound'))
df <- df[,-c(4:6,8)]
mat <- misha_extract(tracks=c(atac_tracks),track_names = c('GFP','Ngn2','PmutNgn2'),regions=df,window = 0,iterator=df,mode='avg')
mat$intervalID <- paste0(mat$chrom,mat$start)
df$intervalID <- paste0(df$chrom,df$start)
df <- merge(df,mat[,-c(1:3)],by='intervalID',sort=F)
df <- df[,-c(2:4)]
df$direction <- df$Yy1_direction
df$direction[df$Ngn2_direction=='NotBound'] <- 'NotBound'
df_o <- melt(df[,-c(1:3)],id.vars=c('direction'))
df_o$variable <- factor(df_o$variable,levels=c('GFP','Ngn2','PmutNgn2'))
df_o$direction <- factor(df_o$direction,levels=c('GFP','shared','Ngn2','NotBound'))
p <- ggplot(df_o,aes(x=direction,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = T,width=0.8)
p <- p + scale_fill_manual(name='',values=sample_colors[2:4]) + xlab('') + ylab('Normalized ATAC signal') + theme(legend.position = "none")  + ggtitle('')
p <- p + theme(legend.position=c(0.2,(1.02)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 0.5),nrow = 1))

pdf('plots/figures/Figure7F_new.pdf',width=5,height=5,useDingbats = F)
print(p+coord_cartesian(ylim=c(0,6.3)))
dev.off()


peaks1 <- read.table('CutRun/full/Ngn2vsGFP_Yy1_LFC2.tsv',header=T)[,c(1:3,11)]
peaks2 <- read.table('data/ChIP/Ngn2_peaks5000_distal_Ngn2Motif.bed',header=T)[,c(1:3)]
colnames(peaks1) <- c('chrom','start','end','Yy1_direction')
colnames(peaks2) <- c('chrom','start','end')
peaks2 <- peaks2[!duplicated(peaks2),]
df <- gintervals.neighbors(peaks2,peaks1)
df$Ngn2_direction[df$dist<=100] <- 'Bound'
df$Ngn2_direction[df$dist>100] <- 'NotBound'
df$Ngn2_direction <- factor(df$Ngn2_direction,levels=c('NotBound','Bound'))
df$Yy1_direction <- factor(df$Yy1_direction,levels=c('GFP','shared','Ngn2','PmutNgn2','NotBound'))
df <- df[,-c(4:6,8)]
mat <- misha_extract(tracks=c(meth_tracks),track_names = c('GFP','Ngn2','PmutNgn2'),regions=intervals.centers(df),window = 100,iterator=intervals.centers(df),mode='avg')
mat$intervalID <- paste0(mat$chrom,mat$start)
df$intervalID <- paste0(intervals.centers(df)$chrom,intervals.centers(df)$start)
df <- merge(df,mat[,-c(1:3)],by='intervalID',sort=F)
df <- df[,-c(2:4)]
df$direction <- df$Yy1_direction
df$direction[df$Ngn2_direction=='NotBound'] <- 'NotBound'
df_o <- melt(df[complete.cases(df),-c(1:3)],id.vars=c('direction'))
df_o$variable <- factor(df_o$variable,levels=c('GFP','Ngn2','PmutNgn2'))
df_o$direction <- factor(df_o$direction,levels=c('GFP','shared','Ngn2','NotBound'))
p <- ggplot(df_o,aes(x=direction,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
p <- p + scale_fill_manual(values=sample_colors[2:4]) + xlab('') + ylab('CpG Methylation (%)') + theme(legend.position = "none")  + ggtitle('')
pdf('plots/figures/Figure7G_new.pdf',width=5,height=5,useDingbats = F)
print(p)
dev.off()

pdf('plots/figures/Figure7H_new_YY1_5000_8kb_1.pdf',width=5,height=3.5)
layout(matrix(c(1:5,5),nrow=2,ncol=3,byrow=F),widths = c(4,4,2),heights=c(4,4),respect = T)
params <- plot_aggregateHiC(cells=cells[2],pool=T,intervals1='Ngn2_peaks5000_Yy1Notbound.bed',intervals2='Ngn2_peaks5000_Yy1Notbound.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=8000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(1),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T,plot_mean = F)
params <- plot_aggregateHiC(cells=cells[2],pool=T,intervals1='data/ChIP/Ngn2_peaks5000_Yy1bound.bed',intervals2='data/ChIP/Ngn2_peaks5000_Yy1bound.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=8000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(1),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T,plot_mean = F)
params <- plot_aggregateHiC(cells=cells[3],pool=T,intervals1='data/ChIP/Ngn2_peaks5000_Yy1Notbound.bed',intervals2='data/ChIP/Ngn2_peaks5000_Yy1Notbound.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=8000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(1),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T,plot_mean = F)
params <- plot_aggregateHiC(cells=cells[3],pool=T,intervals1='data/ChIP/Ngn2_peaks5000_Yy1bound.bed',intervals2='data/ChIP/Ngn2_peaks5000_Yy1bound.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=8000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(1),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T,plot_mean = F)
par(mar=c(2.5,1.5,2.5,4))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

p <- Figure7F_new(score_f1='/home/hpc/bonev/projects/hic/dfg/data/cis_decay/Ngn2_peaks5000_Yy1Notbound.bed_Ngn2_peaks5000_Yy1Notbound.bed.1D.5000_score',
                   score_f2='/home/hpc/bonev/projects/hic/dfg/data/cis_decay/Ngn2_peaks5000_Yy1bound.bed_Ngn2_peaks5000_Yy1bound.bed.1D.5000_score',
                   min_dist = 5e4,max_dist = 2e6,labels=c('GFP','Ngn2'),cols=sample_colors[2:3],p_stats=list(c('GFP','Ngn2')))
pdf('plots/figures/Figure7I_new.pdf',height=5,width=4)
print(p)
dev.off()

tracks=c('data/ChIP/Ngn2_chip_D2.bw',
         'CutRun/full/wt_ngn2_YY1/rep1/peakcalling/macs2.broad/wt_ngn2_YY1_rep1.cpm.norm.bw',
         'CutRun/full/wt_ngn2_Cohesin/merge/wt_Ngn2_Rad21.bw',
         'results/SC/macs2/bdg/Ngn2.bw',
         'CutRun/full/wt_ngn2_h3k27ac/rep1/aligned/dedup/wt_Ngn2_H3K27ac.bw')
peaks='results/ChIP/Ngn2_top5000_peaks_Yy1.bed'
peaks_anno <- read.table('results/ChIP/Ngn2_top5000_peaks_Yy1.bed',header=F)

hm_l <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=1000,bin=20,clusters=peaks_anno$V5,cols=rev(colorpalette('reds')),zlims = matrix(c(0,1,0,2,0,0.6,0,4,0,0.6),byrow = T,ncol=2),raster_quality = 5,show_leg = F)
pdf('plots/figures/Figure7J_new.pdf',width=11,height=6)
draw(hm_l$hm_l[[1]] + hm_l$hm_l[[2]] + hm_l$hm_l[[3]]+ hm_l$hm_l[[4]]+ hm_l$hm_l[[5]],merge_legends=F)
dev.off()



Figure7E(object,cols=yy1_colors,out_f='Figure7E',point.size=1.5,anno.size=12,key.size=3,height=6,width=5.5,plot_filled=T,theme = theme_border,rows=1,stroke=0.2)
Figure7F(object,cols=yy1_clusters,out_f='Figure7F',point.size=1.5,anno.size=12,key.size=3,height=6,width=5.5,plot_filled=T,theme = theme_border,rows=2,stroke=0.2)
Figure7G(object,features='seurat_clusters',cols=yy1_clusters,out_f='Figure7G',height=7,width=6)
DefaultAssay(object) <- 'RNA'
Figure7H(object,features=c(as.vector(cell_markers[,1]),'Pdgfra','Foxj1','Sox10'),cols=c('grey','red','black'),out_f='Figure7H',height=5,width=9)
Figure7I(object=object,c('Gfap','Sox11','Dcx'),out_f='Figure7I',cols=gene_colors(50),point.size=1.5,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)
p <- Figure7J(object,clusters=rev(c('Yy1_WT','Yy1_KO','Yy1_WT/Ngn2+','Yy1_KO/Ngn2+')),group.by = 'orig.ident',showGO=6,out_f='Figure7J',height=8,width=16)
  
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

#FigureS7

FigureS7A(object,features='seurat_clusters',cols=yy1_colors,out_f='FigureS7A',height=7,width=6)
FigureS7B(object,features='nCount_RNA',cols=yy1_colors,out_f='FigureS7B',height=7,width=6)
FigureS7C(object,features='nCount_RNA',cols=yy1_clusters,out_f='FigureS7C',height=6,width=8)
FigureS7D(object,features='nFeature_RNA',cols=yy1_colors,out_f='FigureS7D',height=7,width=6)
FigureS7E(object,features='nFeature_RNA',cols=yy1_clusters,out_f='FigureS7E',height=6,width=8)
features <- c('Yy1','Dcx','Sox11','Tubb3','Igfbp1l','Map2','Elavl2','Rnd2','Mapt','Sox4','Sox9','Gfap','Aldoc','Hes6','Id3','Slc1a3','Cdkn1a','Hspe1','Mdm2','Bax')

FigureS7G(object,which_features='iN_3',features=cell_markers[,1],point.size=1,scale=T,cols=yy1_colors,cols2=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))),anno.size=12,out_f='FigureS7G_iN3',width=6,height=8)

Figure7I(object=object,c('Atf5','Hspe1','Hspd1'),out_f='FigureS7G',cols=gene_colors(50),point.size=1.5,width=6,height=6,plot_filled=T,direction='horizontal',anno.size=8,theme=theme_border,stroke=0.2)

p <- VlnPlot(object, features = c('Hspe1'),cols=yy1_colors,pt.size = 0,combine = F,group.by='orig.ident')
p1 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')
pdf('plots/figures/FigureS4G.pdf',height=6,width=4)
p1
dev.off()


library(ChIPpeakAnno)

peaks1 <- import('ChIP_Nov22/Data/Peaks/Ngn2_conservative_peaks.narrowPeak')
peaks2 <- import('ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak')
peaks3 <- import('CutRun/full/wt_ngn2_YY1/rep1/peakcalling/macs2.narrow/wt_ngn2_YY1_rep1_peaks.narrowPeak')
peaks4 <- import('CutRun/full/wt_pmut_ngn2_YY1/rep1/peakcalling/macs2.narrow/wt_pmut_ngn2_YY1_rep1_peaks.narrowPeak')

mcols(peaks1) <- NULL
mcols(peaks2) <- NULL
mcols(peaks3) <- NULL
mcols(peaks4) <- NULL
ol <- findOverlapsOfPeaks(unique(peaks2), unique(peaks4),maxgap = 500)
makeVennDiagram(ol,
                col=c("black", "black"))

