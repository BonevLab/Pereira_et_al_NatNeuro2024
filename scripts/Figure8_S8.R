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
library(tglkmeans)

library(monaLisa)
library(Biostrings)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Hmisc)
library(plyr)
library(dplyr)
library(seqplots)

## Set main folder as current dir ####

library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3)

source('scripts/figures/config.R')
source('scripts/figures/plot_functions.R')
source('scripts/SC/aux_functions.R')
source('scripts/aux_functions.R')

source('/home/hpc/bonev/projects/hic/dfg/config.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/main_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/aux_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/plot_functions.R')


object <- readRDS('results/SC/Yy1_SCT_Cell_type.RDS')
cell_markers <- read.table('data/SC/cell_markers.tsv')

misha_tracks <- c("atac.reproYy1_Yy1plus","atac.reproYy1_Yy1minus","atac.reproYy1_Yy1plus_Ngn2plus","atac.reproYy1_Yy1minus_Ngn2plus")
track_names <- c('YY1_WT','YY1_KO','YY1_WT/Ngn2+','YY1_KO/Ngn2+')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))


Figure8D <- function(poscor_mat_f,peaks_f){
  poscor_mat <- read.table(poscor_mat_f,header=T)
  poscor_mat$cluster <- c('AST','iN_1','iN_2')[max.col(poscor_mat[,c('AST','iN_1','iN_2')])]
  peaks <- read.table(peaks_f,header=T)
  peaks <- peaks[,c(1:3,ncol(peaks))]
  colnames(peaks) <- c('chrom','start','end','peak_type')  
  distal_a <- as.data.frame(stringr::str_split(paste0(poscor_mat$peakName), pattern ='_' , n = 3, simplify = TRUE))
  colnames(distal_a) <- c('chrom','start','end')
  distal_a$gene <- poscor_mat$gene_name
  distal_a$clust <- poscor_mat$cluster
  distal_a$region_type <- 'distal'
  distal_a <- distal_a[!duplicated(distal_a),]
  prom_a <- data.frame(chrom=poscor_mat$gene_chr,start=poscor_mat$gene_start,end=poscor_mat$gene_start+1,gene=poscor_mat$gene_name,clust=poscor_mat$cluster,region_type='prom')
  prom_a <- prom_a[!duplicated(prom_a),]
  mat <- rbind(distal_a,prom_a)
  mat$start <- as.integer(mat$start)
  mat$end <- as.integer(mat$end)
  df <- gintervals.neighbors(mat,peaks)
  df <- df[,-c(1:3,7:9)]
  df$peak_type[df$peak_type=='shared'] <- 'Const'
  df$peak_type[df$peak_type!='Const'] <- 'Diff'
  df$peak_type[abs(df$dist)>=500] <- 'NB'
  df$peak_type <- factor(df$peak_type,levels=c('Diff','Const','NB'))

  p <- ggplot(df, aes(x=clust,y=1,fill=peak_type)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='',values=c("#8C510A","#01665E","grey80"))
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') 
  p <- p + facet_wrap(~region_type) + theme(legend.text = element_text(size=12)) + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
  pdf('plots/figures/Figure8D_new.pdf',width=5,height=5,useDingbats = F)
  print(p)
  dev.off()
  return(df)
}

Figure8E <- function(df,type='distal',logFC_thresh=0.25){
  df$direction[(df$FDR<=0.05)&(df$logFC<=(-logFC_thresh))] <- 'down'
  df$direction[(df$FDR<=0.05)&(df$logFC>logFC_thresh)] <- 'up'
  df <- df[,c("gene","clust",'region_type','peak_type','direction')]
  df <- df[!duplicated(df),]
  df <- df[df$region_type%in%type,]
  df$direction <- factor(df$direction,levels=c('up','n.s.','down'))
  #  mat <- melt(df,id.vars=c('direction','rank','nearestGene'))
  #  mat$direction <- factor(mat$direction,levels=c('Ngn2','shared','PmutNgn2'))
  #  mat$variable <- factor(gsub('vsGFP','',mat$variable),levels=c('Ngn2','PmutNgn2'))
  # p1 <- ggplot(mat,aes(x=variable,y=value,fill=direction)) + geom_boxplot(outlier.size=1,outlier.shape=NA,show.legend = T,width=0.8) 
  # p1 <- p1 + scale_fill_manual(values=cols,name='') + xlab('') + ylab('log2FC') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.35, 0.9))
  #  p1 <- p1 + stat_compare_means(label = "p.format",method='wilcox',paired = F)
  p <- ggplot(df, aes(x=clust,y=1,fill=direction)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='',values = c('red','grey','blue')) 
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')  
  return(p)
}


Figure8G <- function(peaks1_f,peaks2,peaks1_levels=c('WT','shared','KO'),peaks2_levels=c('Yy1WT','shared','Yy1KO')){
  peaks1 <- read.table(peaks1_f,header=T)
  peaks2 <- read.table(peaks2_f,header=F)
  peaks1 <- peaks1[,c(1:3,ncol(peaks1))]
  peaks2 <- peaks2[,c(1:3)]
  colnames(peaks1)[c(1,4)] <- c('chrom','peaks1_direction')
  peaks1$peaks1_direction <- factor(peaks1$peaks1_direction,levels=peaks1_levels)
  colnames(peaks2)[1:3] <- c('chrom','start','end')
  df <- gintervals.neighbors(peaks1,peaks2)
  df$peaks2_direction <- 'NotBound'
  df$peaks2_direction[df$dist<=500] <- 'Bound'
  df <- df[,c('peaks1_direction','peaks2_direction')]
  df <- df[order(df$peaks1_direction,df$peaks2_direction),]
  p <- ggplot(df, aes(x=peaks1_direction,y=1,fill=peaks2_direction)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='',values=c("#737373","#E6E6E6")) 
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')  
  return(p)
}

Figure8G <- function(object_cpm,genes,track_names,cols,comparisons,label.y,out_f,height,width,ylim){
  object_cpm <- readRDS(object_cpm)
  df <- object_cpm[row.names(object_cpm)%in%genes,]
  df <- as.matrix(df[complete.cases(df),])
  colnames(df) <- gsub('pseudo_','',colnames(df))
  colnames(df) <- gsub('Yy1_','YY1_',colnames(df))
  mat <- melt(df)
  mat$Var2 <- factor(mat$Var2,levels = track_names)
  p <- ggplot(mat,aes(x=Var2,y=value,fill=Var2)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8) 
  p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('Normalized Expression') + theme(legend.position = "none")
  p <- p + coord_cartesian(ylim=ylim) + stat_compare_means(comparisons = list(comparisons),label = "p.format",method='wilcox',paired = T,label.y = label.y,tip.length = c(0.02,0.01))
  p <- p + theme(axis.text.x=element_blank())
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}


### Plot Figures
tracks <- c("data/ATAC//WT_RFP_merge.bw","data/ATAC//Homo_RFP_merge.bw","data/ATAC//Homo_GFP_merge.bw","data/ATAC/Homo_Double_merge.bw")
track_names <- c('YY1_WT','YY1_KO','YY1_WT/Ngn2+','YY1_KO/Ngn2+')
peaks='results/ATAC/Ngn2plusYy1KOvsNgn2plusYy1WT_at_Ngn2peaks.bed_2.bed'
peaks_anno <- read.table(peaks,header=F)
hm_l <- seqplots_heatmap(tracks=tracks[3:4],peaks=peaks,window=1000,bin=20,clusters=peaks_anno[,5],cols=rev(colorpalette('reds')),zlims = c(0,4),raster_quality = 5,show_leg = F)
hm_l$hm_l[[1]]@left_annotation <- NULL
pdf('plots/figures/Figure8B_new2.pdf',width=5,height=4)
draw(hm_l$hm_l[[1]] + hm_l$hm_l[[2]] ,merge_legends=F)
dev.off()




peaks1 <- read.table('results/ATAC/Yy1KOvsYY1WT_at_GFP_Yy1peaks.bed')
peaks1 <- peaks1[,c(1:3,5)]
colnames(peaks1) <- c('chrom','start','end','type')
df <- gintervals.neighbors(peaks1,gintervals.load(tss_f))
df$g_type <- 'Distal'
df$g_type[abs(df$dist)<=5000] <- 'Promoter'
df$type <- as.factor(df$type)
levels(df$type) <- c('Yy1 WT only','constant','Yy1 KO only')
df$type <- factor(df$type,levels=rev(levels(df$type)))
df <- df[df$type!='Yy1 KO only',c('type','g_type')]
p <- ggplot(df,aes(x=type,y=1,fill=g_type)) + geom_bar(position="fill", stat="identity") + scale_fill_grey(end = 0.2,start=0.9)
p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') 
p <- p + theme(legend.text = element_text(size=12),axis.text.x = element_blank()) + theme(legend.position=c(0.25,1.05), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
pdf('plots/figures/Figure8C_new2.pdf',width=4.5,height=2.5,useDingbats = F)
print(p + coord_flip() + scale_x_discrete(labels=label_wrap(8))+ylab(''))
dev.off()


tracks <- c("data/ATAC//WT_RFP_merge.bw","data/ATAC//Homo_RFP_merge.bw","data/ATAC//Homo_GFP_merge.bw","data/ATAC/Homo_Double_merge.bw")
track_names <- c('YY1_WT','YY1_KO','YY1_WT/Ngn2+','YY1_KO/Ngn2+')
peaks='results/ATAC/Ngn2plusYy1KOvsNgn2plusYy1WT_at_Ngn2peaks.bed_2.bed'
peaks_anno <- read.table(peaks,header=F)
hm_l <- seqplots_heatmap(tracks=tracks[3:4],peaks=peaks,window=1000,bin=20,clusters=peaks_anno[,5],cols=rev(colorpalette('reds')),zlims = c(0,4),raster_quality = 5,show_leg = F)
hm_l$hm_l[[1]]@left_annotation <- NULL
pdf('plots/figures/Figure8B_new2.pdf',width=5,height=4)
draw(hm_l$hm_l[[1]] + hm_l$hm_l[[2]] ,merge_legends=F)
dev.off()

#Figure8B
p <- Figure8E(df=res,type=c('distal','prom'),logFC_thresh=0.25)
pdf('plots/figures/Figure8E_new.pdf',width=4,height=5)
print(p+ theme(legend.text = element_text(size=12)) + theme(legend.position=c(0.25,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1)))
dev.off()

tracks <- c("CutRun/full/wt_ngn2_cre_Flag/merge/wt_ngn2_cre_Flag.bw","CutRun/full/cKO_ngn2_cre_Flag/merge/cKO_ngn2_cre_Flag.bw")
peaks='CutRun/full/Yy1KOvsWT_Ngn2FlagLFC2_2.bed'
peaks_anno <- read.table(peaks,header=F)
hm_l <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=1000,bin=20,clusters=peaks_anno$V5,cols=rev(colorpalette('reds')),zlims = c(0,1),raster_quality = 5,show_leg = F)
hm_l$hm_l[[1]]@left_annotation <- NULL
pdf('plots/figures/Figure8F_new_2.pdf',width=4.5,height=4)
draw(hm_l$hm_l[[1]]+ hm_l$hm_l[[2]] ,merge_legends=F)
dev.off()

library(seqplots)
res <- getPlotSetArray(tracks=tracks,features='results/beds/YY',refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10)
pdf('plots/figures/FigureS8A.pdf',height=6,width=6)
plotAverage(plotset=res,labels=track_names, xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = yy1_colors[3:4], pointsize = 12)
axis(side = 1,labels=c('-1000','-500','TSS','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=0.065,tick = T)
dev.off()



source('/home/hpc/bonev/projects/hic/dfg/config.R')
p2glinks <- readRDS('results/P2G-Links.RDS')
plotMisha(main_f=main_f,targetGene='Igfbpl1',outDir='plots/figures/',out_f='Figure8I_new',upstream=1.5e5,downstream=2.4e4,
          window_scale=2.2,pointCEX=0.5,conditions=score_tracks[2:4],
          chipNames=c(rep('',8)),
          chipTracksToExtract=c(misha_tracks,"chipseq_RPM.iN_Yy1wtNgn2_Ngn2Flag","chipseq_RPM.iN_Yy1koNgn2_Ngn2Flag","chipseq_RPM.iN_GFP_Yy1","chipseq_RPM.iN_Ngn2_Yy1"),
          arcIntervals=p2glinks$posCor,arcColors=colorRampPalette(ArchR::ArchRPalettes[[29]]),
          chipYlim=matrix(c(0,4.7,0,4.7,0,4.7,0,4.7,0,1.6,0,1.6,0,3,0,3),ncol = 2,byrow = T),
          chipColors=c("#89C75F","#8A9FD1","#7E1416","#D8A767","black","black","black","black"),img_factor=1,
          plotOrder=list(scores=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE ,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=FALSE),
          plotRatios=list(unitHeight=100, scores=2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=1, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))

p <- VlnPlot(object, features = c('Igfbpl1'),cols=yy1_colors,pt.size = 0,combine = F,group.by = 'orig.ident')
p1 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')
pdf('plots/figures/Figure8J.pdf',height=5,width=3.5)
p1 + theme(axis.text.x = element_blank())
dev.off()

###Figure S8
library(seqplots)
res <- getPlotSetArray(tracks=tracks,features="data/mm10_TSS.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10)
pdf('plots/figures/FigureS8A.pdf',height=6,width=6)
plotAverage(plotset=res,labels=track_names, xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = yy1_colors, pointsize = 12)
axis(side = 1,labels=c('-1000','-500','TSS','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=0.065,tick = T)
dev.off()

res <- getPlotSetArray(tracks=tracks,features="results/ATAC/YY1_ATAC_clusters_new_avg_k5.bed",refgenome='mm10',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 20,)
pdf('plots/figures/FigureS8B.pdf',width=8,height=12)
plotHeatmap(res,labels = track_names,sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod="bed_scores",raster = T,indi=F,ggplot=T)
dev.off()



res <- getPlotSetArray(tracks=tracks,features="results/beds/allATAC_ZBTB18Motifs.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10)
pdf('plots/figures/FigureS8B.pdf',height=6,width=6)
plotAverage(plotset=res,labels=track_names, xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = 'Average Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = yy1_colors, pointsize = 12)
axis(side = 1,labels=c('-1000','-500','CTCF ->','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=0.06,tick = T)
dev.off()


df_all=read.table('results/P2G_ClusterBinaryMat_posCor.tsv',header=T)
df_all$cluster <- c('AST','iN_1','iN_2')[max.col(df_all[,c('AST','iN_1','iN_2')])]
df <- df_all[df_all$cluster=='iN_1',]
features <- as.vector(cell_markers[13:23,1])
features2=c("Rbfox3:247336","Sox9:-280433","Sox9:234828","Gfap:-154240","Sox11:277561","Sox11:-233215","Sox4:-435166","Sox4:235153","Hes6:348009","Elavl2:132949","Sema5a:-344033","Sema5a:101208","Dll1:127766",
           "Slc1a3:-105947","Slc1a3:38106","Cdkn1c:59587","Rnd2:17320","Auts2:23090","Auts2:409859","Auts2:-358036","Gria1:203029","Gria1:-455123","Slit1:131953","Slit1:46416")
mat <- as.data.frame(stringr::str_split(paste0(df$peakName), pattern ='_' , n = 3, simplify = TRUE))
colnames(mat) <- c('chrom','start','end')
mat <- unique(mat)
mat <- gintervals(mat[,1],as.numeric(mat[,2]),as.numeric(mat[,3]))
mat <- misha_extract(tracks=misha_tracks,track_names = track_names,regions=mat,window = 0,iterator=mat,mode='avg')
mat$ID <- paste0(mat$chrom,'_',mat$start,'_',mat$end)
res <-  merge(df,mat[,-c(1:3)],by.x='peakName',by.y='ID',sort=F)
mat <- as.matrix(res[,grep('YY1',colnames(res))])
mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
mat_k <- tglkmeans::TGL_kmeans(mat,k = 4,id_column = F,reorder_func = 'hclust',hclust_intra_clusters = F,seed = 42)
res$YY1_cluster <- mat_k$cluster
mat_factor <- factor(mat_k$cluster,levels=c(1,3,4,2))
levels(mat_factor) <- c(1,3,4,2)

ha <- columnAnnotation(foo = anno_text(colnames(mat), location = 0.5,rot=0, just = "center",gp = gpar(fill = yy1_colors, col = "black", border = "black", fontsize = 12),height = unit(0.3,'inch')))
la = rowAnnotation(foo = anno_mark(at = which(res$gene_name%in%features),side='left', labels = res$labels[which(res$gene_name%in%features)],labels_gp = gpar(fontsize = 12)))
la = rowAnnotation(foo = anno_mark(at = which(res$labels%in%features2),side='left', labels = res$labels[which(res$labels%in%features2)],labels_gp = gpar(fontsize = 12)))

hm <- Heatmap(as.matrix(mat),row_title = NULL,row_split = mat_factor,show_row_names = F,cluster_row_slices = F,cluster_rows = F,left_annotation = la,show_row_dend = F,cluster_columns = F,show_column_names = F,use_raster = TRUE, raster_quality = 10, 
              col = colorRamp2(seq(-1.2,1.2,length.out=length(heatmap_colors)), heatmap_colors),top_annotation = ha,
                border=T,heatmap_legend_param=list(at=c(-1.2,0,1.2),direction = "vertical",title='Accessibility Z-scores'))
pdf('plots/figures/FigureS8F_iN1.pdf',height=12,width=8)
draw(hm)
dev.off()

p <- VlnPlot(object, features = c('Sema7a'),assay = 'SCT',cols=yy1_colors,group.by = 'orig.ident',pt.size = 0,combine = F)
p1 <- p[[1]] + theme(legend.position='none',axis.text.x = element_text(angle = 0, hjust = 0.5))+ xlab('')
p1


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

df$cluster <- factor(paste0(df$Yy1_direction,'_',df$Ngn2_direction),
                     levels=c('GFP only_No','GFP only_Yes','shared_No','shared_Yes','Ngn2 only_No','Ngn2 only_Yes'))
tracks <- c("data/ATAC//WT_RFP_merge.bw","data/ATAC//Homo_RFP_merge.bw","data/ATAC//Homo_GFP_merge.bw","data/ATAC/Homo_Double_merge.bw")
track_names <- c('YY1_WT','YY1_KO','YY1_WT/Ngn2+','YY1_KO/Ngn2+')
peaks='results/ChIP/Ngn2vsGFP_Yy1_LFC2.bed'
hm_l <- seqplots_heatmap(tracks=tracks[3:4],peaks=peaks,window=1000,bin=20,clusters=as.numeric(df$cluster),cols=rev(colorpalette('reds')),zlims = c(0,4),raster_quality = 5,show_leg = F)
hm_l$hm_l[[1]]@left_annotation <- NULL
pdf('plots/figures/FigureS8B_Magd.pdf',width=5,height=8)
draw(hm_l$hm_l[[1]] + hm_l$hm_l[[2]] ,merge_legends=F)
dev.off()



pdf('plots/figures/Figure7C_new.pdf',width=5,height=4,useDingbats = F)
print(p+coord_flip())
dev.off()



