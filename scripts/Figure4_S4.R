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

Figure4C <- function(peaks,logFC,bins,central_bin,breaks,pwms,genome_bg=F,genome=BSgenome.Mmusculus.UCSC.mm10,peak_size=200,features=NULL,mcparams=BiocParallel::MulticoreParam(10L),FDR.cutoff=4,enr.cutoff=1,cols,anno.size=8){
  colnames(peaks) <- c('chrom','start','end')
  attributes(bins)$bin0 <- central_bin
  attributes(bins)$breaks <- breaks
  peaks <- makeGRangesFromDataFrame(peaks)
  peakseqs <- getSeq(genome, resize(peaks,width = peak_size,fix = 'center'))
  mcparams <- BiocParallel::MulticoreParam(10L)
  se <- calcBinnedMotifEnrR(seqs = peakseqs, bins = bins,background = ifelse(genome_bg,"genome","zeroBin"), pwmL =toPWM(pwms),BPPARAM = mcparams,genome=genome)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  seSel <- se[sel1&sel2, ]
  df <- assay(seSel,"log2enr")
  row.names(df) <- rowData(seSel)$motif.name
  df <- df[ order(max.col(df, "first")), ]
  col.list=rev(colorpalette('brbg',max(as.numeric(bins))))
  names(col.list) <- 1:max(as.numeric(bins))
  if(!is.null(features)){
    la = rowAnnotation(foo = anno_mark(at = which(row.names(df)%in%features), labels = row.names(df)[which(row.names(df)%in%features)],labels_gp = gpar(fontsize = anno.size)))
  } else {
    la = rowAnnotation(foo = anno_mark(at = 1:nrow(df), labels = row.names(df),labels_gp = gpar(fontsize = anno.size)))
  }
  #SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm,BPPARAM = mcparams)
  #hcl <- hclust(as.dist(1 - SimMatSel), method = "ward.D2")
  ta <- HeatmapAnnotation(cluster = names(col.list) ,col = list(cluster=col.list),show_legend = F,show_annotation_name=F,border=T)
  hm <- Heatmap(df,name='log2enr',cluster_rows = T,show_row_names = F,right_annotation = la,cluster_columns = F,show_column_names = F,show_row_dend = T, col = cols,top_annotation = ta,border=T,
                heatmap_legend_param=list(direction = "horizontal",title='log2enr',legend_width = unit(2, "inch"),xjust=0.5,yjust=0.5))
  sekm <- calcBinnedKmerEnr(seqs = peakseqs, bins = bins, kmerLen = 6,background = ifelse(genome_bg,"genome","zeroBin"),BPPARAM = mcparams,genome=genome,
                            includeRevComp = TRUE)
  sel1 <- apply(assay(sekm, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(sekm, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  sekmSel <- sekm[sel1&sel2, ]
  df2 <- assay(sekmSel,"log2enr")
  row.names(df2) <- rowData(sekmSel)$motif.name
  if(!is.null(features)){
    la = rowAnnotation(foo = anno_mark(at = which(row.names(df)%in%features), labels = row.names(df)[which(row.names(df)%in%features)],labels_gp = gpar(fontsize = anno.size)))
  } else {
    la = rowAnnotation(foo = anno_mark(at = 1:nrow(df), labels = row.names(df),labels_gp = gpar(fontsize = anno.size)))
  }
  sims <- motifKmerSimilarity(x = rowData(sekmSel)$motif.pfm,kmers = rownames(sekmSel),includeRevComp = TRUE)
  hcl <- hclust(as.dist(1 - sims), method = "ward.D2")
  hm2 <- Heatmap(df2,name='log2enr',cluster_rows = hcl,show_row_names = F,cluster_columns = F,show_column_names = F,show_row_dend = T, col = cols,top_annotation = ta,border=T,width=2,
                heatmap_legend_param=list(direction = "horizontal",title='log2enr'))
  sims <- motifKmerSimilarity(x = rowData(se)$motif.pfm,kmers = rownames(sekmSel),includeRevComp = TRUE)
  sims_idx <- apply(sims,MARGIN = 2,function(x){
    return(which(x==max(x)))
  })
  sims <- sims[sims_idx,]
  pre_logo <- rowData(se)$motif.pfm
  pre_logo <- pre_logo[match(row.names(sims),rowData(se)$motif.name)]
  maxwidth <- max(sapply(TFBSTools::Matrix(pre_logo), ncol))
  seqlogoGrobs <- lapply(pre_logo, seqLogoGrob, xmax = maxwidth)
  hmSeqlogo <- columnAnnotation(logo = annoSeqlogo(seqlogoGrobs, which = "column"),
                             annotation_width = unit(1.5, "inch"), 
                             show_annotation_name = FALSE)
  sims <- t(sims)
  sims <- sims[match(row.names(df2),row.names(sims)),]
  hm3 <- Heatmap(sims, 
                 show_row_names = TRUE, row_names_gp = gpar(fontsize = 8),
                 show_column_names = TRUE, column_names_gp = gpar(fontsize = 8),
                 name = "Similarity",width=4,
                 col = colorRamp2(c(0, 1), c("white", "red")),heatmap_legend_param=list(direction = "horizontal"))
  return(list(hm=hm,hm_list=hm2+hm3))
}

Figure4E <- function(res,max_dist,split_by_enrichment=2,split_by_type=FALSE,show_which,rank_names=c('weak','strong')){
  df <- res[res$dist<=max_dist,c('logCPM','direction','type','nearestGene',show_which)]
  colnames(df)[ncol(df)] <- 'diff'
  if (!is.null(split_by_enrichment)){
    df$rank <- as.factor(ntile(x = df$logCPM, split_by_enrichment))
    levels(df$rank) <- rank_names
  } else {
    df$rank <- 1
  }
  df <- df[,-1]
  if(!split_by_type){
    df <- df[,-2]
  }
  df <- df[!duplicated(df),]
  df$direction <- factor(df$direction,levels=c('Ngn2','shared','PmutNgn2'))
  df$diff <- factor(df$diff,levels=c('up','n.s.','down'))
#  mat <- melt(df,id.vars=c('direction','rank','nearestGene'))
#  mat$direction <- factor(mat$direction,levels=c('Ngn2','shared','PmutNgn2'))
#  mat$variable <- factor(gsub('vsGFP','',mat$variable),levels=c('Ngn2','PmutNgn2'))
 # p1 <- ggplot(mat,aes(x=variable,y=value,fill=direction)) + geom_boxplot(outlier.size=1,outlier.shape=NA,show.legend = T,width=0.8) 
 # p1 <- p1 + scale_fill_manual(values=cols,name='') + xlab('') + ylab('log2FC') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.35, 0.9))
#  p1 <- p1 + stat_compare_means(label = "p.format",method='wilcox',paired = F)
  p <- ggplot(df, aes(x=direction,y=1,fill=diff)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='',values = c('red','grey','blue')) 
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')  
  return(p)
}



#### Plot Figures #########

  #Fig 4A
tracks=c('data/ChIP/Ngn2_chip_D2.bw','data/ChIP/PmutNgn2_chip_D2.bw')
peaks='results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.bed'
peaks_anno <- read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.tsv',header=T)

hm_l <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=500,bin=10,clusters=read.table(peaks)$V5,cols=rev(colorpalette('reds')),zlims=c(0,0.8),raster_quality = 5,show_leg = F,show_n=F)
pdf('plots/figures/Figure4A.pdf',width=4,height=8)
draw(hm_l$hm_l$Ngn2_chip_D2 +  hm_l$hm_l$PmutNgn2_chip_D2,merge_legends=F)
dev.off()

tracks=c('results/SC/macs2/bdg/GFP.bw','results/SC/macs2/bdg/Ngn2.bw','results/SC/macs2/bdg/PmutNgn2.bw')
hm_l2 <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=1000,bin=10,clusters=read.table(peaks)$V5,cols=rev(colorpalette('reds')),zlims=c(0,4),idx = hm_l$idx,raster_quality = 5,show_leg = F)
hm_l2$hm_l$GFP@left_annotation <- NULL
pdf('plots/figures/Figure4B.pdf',width=5,height=8)
draw(hm_l2$hm_l$GFP+hm_l2$hm_l$Ngn2+hm_l2$hm_l$PmutNgn2,merge_legends=F)
dev.off()


  #Fig 4C
res <- read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.tsv',header=T)
res_hm <- Figure4C(peaks=res[,1:3],logFC=res$logFC,bins=res$direction,central_bin=2,breaks=NULL,peak_size=200,
         pwms=readRDS('data/combined_pwm_FPKM1.RDS'),genome_bg=F,genome=BSgenome.Mmusculus.UCSC.mm10,
         features=c('NEUROG2(var.2)','BHLHE22(var.2)','REST','TGIF2','NFIB','NFIX(var.2)','NFIC(var.2)','RFX2','HES5','NFKB1','ZFP57'),
         mcparams=BiocParallel::MulticoreParam(10L),FDR.cutoff=4,enr.cutoff=1,cols=colorRamp2(seq(-2,2,length.out=9), rev(colorpalette('rdbu',9))))

pdf('plots/figures/Figure4B.pdf',width=4,height=8)
draw(res_hm$hm,merge_legends=F,heatmap_legend_side = "bottom")
dev.off()

  #Fig 4C
mat <- read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.tsv',header=T)
mat <- mat[,c('direction','type')]
mat$direction <- factor(mat$direction,levels=c('Ngn2','shared','PmutNgn2'))
p <- ggplot(mat, aes(x=direction,y=1,fill=type)) +  geom_bar(position="fill", stat="identity") + scale_fill_grey(start = 0.8,end = 0.2)
p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') 
p <- p + theme(legend.text = element_text(size=12)) + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
pdf('plots/figures/Figure4C.pdf',width=4,height=4,useDingbats = F)
print(p)
dev.off()

  #Fig 4D

peaks_f='results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.bed'

p1 <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Analysis/results_BB/Ngn2_peaks.bed',stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,
                             control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))
p1 <- p1+ ylab('Ngn2 (% of total)') + scale_fill_manual(values=hue_pal(l = 50)(5))+theme(legend.text = element_text(size=7.5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position=c(0.05,1.01),legend.box = "horizontal",plot.margin = unit(c(0,0,0,0), "cm")) + guides(fill=guide_legend(override.aes = list(size = 4),nrow = 1))

p2 <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Analysis/results_BB/shared_peaks.bed',stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,
                              control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))
p2 <- p2+ ylab('shared (% of total)') +scale_fill_manual(values=hue_pal(l = 50)(5))+theme(legend.text = element_text(size=7.5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position=c(0.05,1.01),legend.box = "horizontal",plot.margin = unit(c(0,0,0,0), "cm")) + guides(fill=guide_legend(override.aes = list(size = 4),nrow = 1))

p3 <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Analysis/results_BB/PmutNgn2_peaks.bed',stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,
                              control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))
p3 <- p3 + ylab('PmutNgn2 (% of total)') +scale_fill_manual(values=hue_pal(l = 50)(5))+theme(legend.text = element_text(size=7.5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position=c(0.05,1.01),legend.box = "horizontal",plot.margin = unit(c(0,0,0,0), "cm")) + guides(fill=guide_legend(override.aes = list(size = 4),nrow = 1))
p <- wrap_plots(p1/p2/p3,guides='collect') & theme(legend.position = 'top')

p <- stratifyPeaks_by_Motifs(peaks_f='results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.bed',min.score=9,stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,
                              control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/figures/Figure4D.pdf',width=3.5,height=8.5)
p + facet_wrap(~type,ncol=1) + scale_fill_manual(values=colorpalette('colorblind',5))+theme(legend.text = element_text(size=10),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position=c(0.05,1.07),legend.box = "horizontal",plot.margin = margin(0.5, 0.05, 0.05, 0.05, "inch")) + guides(fill=guide_legend(override.aes = list(size = 4),nrow = 1))
dev.off()

  #Fig 4E

res <- read.table('results/ChIP/PmutNgn2vsNgn2_ChIPvsRNA.tsv',header=T)
p1 <- Figure4E(res,max_dist=5e4,split_by_enrichment=2,split_by_type=FALSE,show_which='Ngn2diff')
p2 <- Figure4E(res,max_dist=5e4,split_by_enrichment=2,split_by_type=FALSE,show_which='PmutNgn2diff')
p1 <- p1 + theme(axis.text.x = element_blank(),legend.text = element_text(size=12)) + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal",plot.margin = margin(0, 0, 0, 0, "cm")) + scale_x_discrete(labels=label_wrap(3))
p2 <- p2 + ylab('') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.text = element_text(size=12),legend.position=c(0.2,(1.04)), legend.box = "horizontal",plot.margin = margin(0, 0, 0, 0, "cm")) + scale_x_discrete(labels=label_wrap(6))
p <- wrap_plots(p1+p2,guides='collect') & theme(legend.position = 'right',legend.justification = 'center')
pdf('plots/figures/Figure4E.pdf',width=8,height=3.5)
print(p)
dev.off()

  #Fig 4F

#tracks=c('CutRun/full/wt_GFP_Cohesin/rep1/wt_GFP_Cohesin_rep1/peakcalling/macs2.narrow/wt_GFP_Cohesin_rep1.cpm.norm.bw',
#         'CutRun/full/wt_ngn2_Cohesin/rep1/peakcalling/macs2.narrow/wt_ngn2_cohesin_rep1.cpm.norm.bw',
#         'CutRun/full/wt_pmut_ngn2_Cohesin/rep1/peakcalling/macs2.narrow/wt_pmut_ngn2_cohesin.cpm.norm.bw')
tracks=c('CutRun/full/wt_GFP_Cohesin/merge/wt_GFP_Rad21.bw',
         'CutRun/full/wt_ngn2_Cohesin/merge/wt_Ngn2_Rad21.bw',
         'CutRun/full/wt_pmut_ngn2_Cohesin/rep1/aligned/dedup/wt_PmutNgn2_Rad21.bw')
hm_l3 <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=1000,bin=10,clusters=read.table(peaks)$V5,cols=rev(colorpalette('reds')),zlims=matrix(c(0,0.4,0,0.4,0,0.2),byrow = T,ncol=2),idx = hm_l$idx,raster_quality = 5,show_leg = F)
pdf('plots/figures/Figure4F.pdf',width=5,height=8)
draw(hm_l3$hm_l[[1]]+hm_l3$hm_l[[2]]+hm_l3$hm_l[[3]],merge_legends=F)
dev.off()
  
  #Fig 4G

tracks=c('CutRun/full/wt_ngn2_h3k27ac/rep1/aligned/dedup/wt_Ngn2_H3K27ac.bw',
         'CutRun/full/wt_pmut_ngn2_h3k27ac/rep1/aligned/dedup/wt_PmutNgn2_H3K27ac.bw')
hm_l4 <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=1000,bin=10,clusters=read.table(peaks)$V5,cols=rev(colorpalette('reds')),zlims=matrix(c(0,0.4,0,0.3),byrow = T,ncol=2),idx = hm_l$idx,raster_quality = 5,show_leg = F)
pdf('plots/figures/Figure4G.pdf',width=3.2,height=8)
draw(hm_l4$hm_l[[1]]+hm_l4$hm_l[[2]],merge_legends=F)
dev.off()

source('/home/hpc/bonev/projects/hic/dfg/config.R')
plotMisha(main_f=main_f,targetGene='Plxna2',outDir='plots/figures/',out_f='Figure4H_alt',upstream=4.5e5,downstream=5e4,
          window_scale=2.2,pointCEX=0.5, binSize=10,conditions=score_tracks[2:4],chipRes=200,
          #chipNames=c('AST','iN_1','iN_2','Ngn2','PmutNgn2','GFP Rad21','Ngn2 Rad21','PmutNgn2 Rad21'),figure_mode=TRUE,
          chipNames=c('','','','','','','','','',''),figure_mode=TRUE,
          chipTracksToExtract=c("scATAC.repro_AST","scATAC.repro_iN_1","scATAC.repro_iN_2","chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2",
                                "chipseq_RPM.iN_GFP_Rad21","chipseq_RPM.iN_Ngn2_Rad21","chipseq_RPM.iN_PmutNgn2_Rad21",
                                "chipseq_RPM.iN_Ngn2_H3K27ac","chipseq_RPM.iN_PmutNgn2_H3K27ac"),
          #arcIntervals=p2glinks$posCor,arcColors=colorRampPalette(ArchR::ArchRPalettes[[29]]),
          chipYlim=matrix(c(0,5,0,5,0,5,0,0.8,0,0.8,0,0.6,0,0.6,0,0.6,0,0.5,0,0.5),ncol = 2,byrow = T),
          chipColors=c("#272E6A","#FA7517","#A30005",rep("black",7)),img_factor=1.2,
          plotOrder=list(scores=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE ,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=FALSE),
          plotRatios=list(unitHeight=100, scores=2, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))




#### Extended Data Figure 4

pdf('plots/figures/FigureS4A.pdf',width=8,height=8)
draw(hm_list,merge_legends=F,heatmap_legend_side = "bottom")
dev.off()

tracks=c(#'data/ChIP/Ngn2_chip_D2.bw',
         #'data/ChIP/PmutNgn2_chip_D2.bw',
         '/home/hpc/bonev/projects/chip/neurog2/Ngn2_12hr_ChIP_Aydin.bw',
         '/home/hpc/bonev/projects/chip/neurog2/Ngn2_48hr_ChIP_Aydin.bw',
         '/home/hpc/bonev/data/ChIPseq/norm_wigs/NPC_Neurog2.bw')
#zlims=round(matrix(c(#0,quantile(res$data$Ngn2_allDatasets_overlap$Ngn2_chip_D2$heatmap,0.95,na.rm=T),
               #0,quantile(res$data$Ngn2_allDatasets_overlap$PmutNgn2_chip_D2$heatmap,0.95,na.rm=T),
#               0,quantile(res$data$Ngn2_allDatasets_overlap$Ngn2_12hr_ChIP_Aydin$heatmap,0.95,na.rm=T),
#               0,quantile(res$data$Ngn2_allDatasets_overlap$Ngn2_48hr_ChIP_Aydin$heatmap,0.95,na.rm=T),
#               0,quantile(res$data$Ngn2_allDatasets_overlap$NPC_Neurog2$heatmap,0.95,na.rm=T)),byrow = T,ncol = 2),2)
hm_l5 <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=500,bin=10,clusters=read.table(peaks)$V5,cols=rev(colorpalette('reds')),zlims=c(0,0.6),idx = hm_l$idx,raster_quality = 5)
pdf('plots/figures/FigureS4C.pdf',width=5,height=8)
draw(hm_l5$hm_l[[1]]+hm_l5$hm_l[[2]]+hm_l5$hm_l[[3]],merge_legends=F)
dev.off()


res <- read.table('results/ChIP/PmutNgn2vsNgn2_ChIPvsRNA.tsv',header=T)
p1 <- Figure4E(res,max_dist=1e5,split_by_enrichment=2,split_by_type=FALSE,show_which='Ngn2diff',rank_names=c('weak','strong'))
p1 <- p1 + facet_wrap(~rank) + theme(legend.position = 'top',legend.justification = 'center')
p2 <- Figure4E(res,max_dist=1e5,split_by_enrichment=2,split_by_type=TRUE,show_which='Ngn2diff',rank_names=c('weak','strong'))
p2 <- p2 + facet_wrap(~type) + theme(legend.position = 'top',legend.justification = 'center') + theme(axis.text.x = element_blank())
p <- wrap_plots(p1/p2,guides='collect') & theme(legend.position = 'top',legend.justification = 'center')
pdf('plots/figures/FigureS4D.pdf',width=4,height=6)
print(p)
dev.off()


peaks=read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.bed')
colnames(peaks)[1:3] <- c('chrom','start','end')
df <- makeGRangesFromDataFrame(peaks[peaks$V5==3,1:3])


p1 <- motif_footprinting(archr_obj=archr_obj,flank=500,motif='NEUROG2(var.2)',regions=makeGRangesFromDataFrame(peaks[peaks$V5==1,1:3]),
                   normalize_by='Subtract',out_f='FigureS4E_1',cols=sample_colors[2:4],
                   group.by='Cond',which_clusters=c("GFP","Ngn2","PmutNgn2"),plot_bias=F,anno.size=12,key.size=4,width=6,height=4.5,returnPlot = T)
p1 <- p1 + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position=c(0.05,1.01),legend.box = "horizontal",plot.margin = unit(c(0,0,0,0), "cm")) + guides(col=guide_legend(override.aes = list(size = 4),nrow = 1)) + ggtitle("") + xlab("") + coord_cartesian(ylim=c(0,6))+ ylab("")
p2 <- motif_footprinting(archr_obj=archr_obj,flank=500,motif='NEUROG2(var.2)',regions=makeGRangesFromDataFrame(peaks[peaks$V5==2,1:3]),
                         normalize_by='Subtract',out_f='FigureS4E_2',cols=sample_colors[2:4],
                         group.by='Cond',which_clusters=c("GFP","Ngn2","PmutNgn2"),plot_bias=F,anno.size=12,key.size=4,width=6,height=4.5,returnPlot = T)
p2 <- p2 + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position=c(0.05,1.01),legend.box = "horizontal",plot.margin = unit(c(0,0,0,0), "cm")) + guides(col=guide_legend(override.aes = list(size = 4),nrow = 1)) + ggtitle("")+ xlab("") + coord_cartesian(ylim=c(0,6))+ ylab("")
p3 <- motif_footprinting(archr_obj=archr_obj,flank=500,motif='NEUROG2(var.2)',regions=makeGRangesFromDataFrame(peaks[peaks$V5==3,1:3]),
                         normalize_by='Subtract',out_f='FigureS4E_3',cols=sample_colors[2:4],
                         group.by='Cond',which_clusters=c("GFP","Ngn2","PmutNgn2"),plot_bias=F,anno.size=12,key.size=4,width=6,height=4.5,returnPlot = T)
p3 <- p3 + theme(legend.position=c(0.05,1.01),legend.box = "horizontal",plot.margin = unit(c(0,0,0,0), "cm")) + guides(col=guide_legend(override.aes = list(size = 4),nrow = 1)) + ggtitle("")+ xlab("") + coord_cartesian(ylim=c(0,6))+ ylab("")

p <- wrap_plots(p1/p2/p3,guides='collect') & theme(legend.position = 'top',legend.justification = 'center',legend.direction = 'horizontal',legend.box = 'horizontal')
pdf('plots/figures/FigureS4E.pdf',width=4,height=8)
print(p)
dev.off()


tracks=c('data/ChIP/Ngn2_chip_D2.bw','data/ChIP/PmutNgn2_chip_D2.bw')
peaks='results/beds/all_posCor_DA.bed'
peaks_anno <- read.table('results/beds/all_posCor_DA.bed',header=F)

hm_l6 <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=1000,bin=10,clusters=peaks_anno$V5,cols=rev(colorpalette('reds')),zlims=c(0,1.2),raster_quality = 10,show_leg = F,type='pf')
pdf('plots/figures/FigureS4F_1.pdf',width=4,height=8)
draw(hm_l6$hm_l[[1]] +  hm_l6$hm_l[[2]],merge_legends=F)
dev.off()

peaks='results/beds/all_posCor_GA.bed'
peaks_anno <- read.table('results/beds/all_posCor_GA.bed',header=F)

hm_l7 <- seqplots_heatmap(tracks=tracks,peaks=peaks,window=2000,bin=20,clusters=peaks_anno$V5,cols=rev(colorpalette('reds')),zlims=c(0,1.2),raster_quality = 5,show_leg = F,type='mf')
pdf('plots/figures/FigureS4F_2.pdf',width=4,height=8)
draw(hm_l7$hm_l[[1]] +  hm_l7$hm_l[[2]],merge_legends=F)
dev.off()



#Fig 4C
mat <- read.table('results/ChIP/Motif_Orientation_res.tsv',header=T)
mat <- mat[,c('direction','type')]
mat$Sample <- factor(gsub('Peaks','',mat$Sample),levels=c('Ngn2','shared','PmutNgn2'))
mat$Motif_Orientation <- factor(mat$Motif_Orientation,levels=c('same','opposite'))
p <- ggplot(mat, aes(x=Sample,y=Percentage,fill=Motif_Orientation)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(values=c('#008b8b','gold2'))
p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') 
p <- p + theme(legend.text = element_text(size=12)) + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
pdf('plots/figures/FigureS4H.pdf',width=4,height=4,useDingbats = F)
print(p)
dev.off()


peaks1 <- read.table('CutRun/full/Ngn2vsGFP_Rad21_LFC2.tsv',header=T)[,c(1:3,11)]
peaks2 <- read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.tsv',header=T)[,c(1:3,11)]
colnames(peaks1) <- c('chrom','start','end','Rad21_direction')
colnames(peaks2) <- c('chrom','start','end','Ngn2_direction')
df <- gintervals.neighbors(peaks1,peaks2)
df$Ngn2_direction[df$dist>500] <- 'NotBound'
df$Ngn2_direction <- factor(df$Ngn2_direction,levels=c('Ngn2','shared','PmutNgn2','NotBound'))
df$Rad21_direction <- factor(df$Rad21_direction,levels=c('Ngn2','shared','GFP'))
df <- df[,c('Rad21_direction','Ngn2_direction')]
p <- ggplot(df, aes(x=Rad21_direction,y=1,fill=Ngn2_direction)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='',values=c("#003C30","#CCCCCC","#543005","#404040"))
p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')
p <- p + theme(legend.position=c(0.02,(1.02)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
pdf('plots/figures/FigureS4M.pdf',width=5,height=6,useDingbats = F)
print(p)
dev.off()




peaks <- read.table('CutRun/full/Ngn2vsGFP_Rad21_LFC2.tsv',header=T)
colnames(peaks)[1] <- 'chrom'
df <- misha_extract(tracks = c("chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2"),regions = peaks[,1:3],window = 0,iterator = peaks[,1:3],track_names = c('Ngn2','PmutNgn2'))
pdf('plots/figures/FigureS4I',width=7,height=8)
draw(hm_l3$hm_l[[1]]+hm_l3$hm_l[[2]]+hm_l3$hm_l[[3]]+hm_l3$hm_l[[4]]+hm_l3$hm_l[[5]],merge_legends=F)
dev.off()

mat <- read.table('data/RNA/Neurog2_UTR_normCounts.tsv',header=T)
df_summ <- data_summary(mat, varname="count",groupnames=c('condition'))
p <- ggplot(df_summ, aes(x=condition, y=count, fill=condition)) + scale_fill_manual(values = sample_colors[2:4]) +
  geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Normalized Counts') +
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.2,position=position_dodge(.9))
p <- p+ geom_point(data = mat,aes(x=condition, y=count, fill=condition),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2) + theme(legend.position = 'none')

pdf('plots/figures/FigureS4J.pdf',width=4,height=4)
print(p)
dev.off()

peaks <- read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.tsv',header=T)
colnames(peaks)[1] <- 'chrom'
df <- misha_extract(tracks = c("chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2"),regions = peaks[,1:3],window = 0,iterator = peaks[,1:3],track_names = c('Ngn2','PmutNgn2'))
df$clust <- peaks$direction[match(paste0(df$chrom,df$start),paste0(peaks$chrom,peaks$start))]
df <- df[,c(4,5,7)]
df <- melt(df,id.vars = 'clust')
df$clust <- factor(df$clust,levels=c('Ngn2','shared','PmutNgn2'))
p <- ggplot(df,aes(x=clust,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = T,width=0.8)
p <- p + scale_fill_manual(values=sample_colors[3:4]) + xlab('') + ylab('ChIP Enrichment') + coord_cartesian(ylim=c(0,0.65))
p <- p + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 1),nrow = 1))

pdf('plots/figures/FigureS4K.pdf',width=4,height=4)
print(p)
dev.off()

peaks <- read.table('results/ChIP/PmutNgn2vsNgn2_Csaw_LFC1.tsv',header=T)
colnames(peaks)[1] <- 'chrom'
df <- misha_extract(tracks = c("scATAC.repro_Astro","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2"),regions = peaks[,1:3],window = 0,iterator = peaks[,1:3],track_names = c('GFP','Ngn2','PmutNgn2'))
df$clust <- peaks$direction[match(paste0(df$chrom,df$start),paste0(peaks$chrom,peaks$start))]
df <- df[,c(4:6,8)]
df <- melt(df,id.vars = 'clust')
df$clust <- factor(df$clust,levels=c('Ngn2','shared','PmutNgn2'))
p <- ggplot(df,aes(x=clust,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = T,width=0.8)
p <- p + scale_fill_manual(values=sample_colors[2:4]) + xlab('') + ylab('Normalized ATAC signal') + coord_cartesian(ylim=c(0,3.7))
p <- p + theme(legend.position=c(0.1,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 1),nrow = 1))

pdf('plots/figures/FigureS4N.pdf',width=4,height=4)
print(p)
dev.off()






p2g_da <- read.table('results/beds/all_posCor_DA.bed')[,c(1:3,5)]
colnames(p2g_da) <- c('chrom','start','end','DA_cluster')
df <- misha_extract(tracks = c("chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2"),regions = p2g_da,window = 500,iterator = p2g_da,track_names = c('Ngn2','PmutNgn2'))
df$clust <- p2g_da$DA_cluster[match(paste0(df$chrom,df$start),paste0(p2g_da$chrom,p2g_da$start))]
df <- df[,c(4,5,7)]
df <- melt(df,id.vars = 'clust')
df$clust <- factor(df$clust)
levels(df$clust) <- c('AST','iN_1','iN_2')
p1 <- ggplot(df,aes(x=clust,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = T,width=0.8)
p1 <- p1 + scale_fill_manual(name='',values=sample_colors[3:4]) + xlab('') + ylab('ChIP Enrichment') + coord_cartesian(ylim=c(0,0.6))
p1 <- p1 + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal",axis.text.x = element_blank(),axis.ticks.x=element_blank()) + guides(fill=guide_legend(override.aes = list(size = 1),nrow = 1))

p2g_ga <- read.table('results/beds/all_posCor_GA.bed')[,c(1:3,5)]
colnames(p2g_ga) <- c('chrom','start','end','GA_cluster')
df <- misha_extract(tracks = c("chipseq_RPM.iN_Ngn2_D2","chipseq_RPM.iN_PmutNgn2_D2"),regions = p2g_ga,window = 500,iterator = p2g_ga,track_names = c('Ngn2','PmutNgn2'))
df$clust <- p2g_ga$GA_cluster[match(paste0(df$chrom,df$start),paste0(p2g_ga$chrom,p2g_ga$start))]
df <- df[,c(4,5,7)]
df <- melt(df,id.vars = 'clust')
df$clust <- factor(df$clust)
levels(df$clust) <- c('AST','iN_1','iN_2')
p2 <- ggplot(df,aes(x=clust,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = T,width=0.8)
p2 <- p2 + scale_fill_manual(name='',values=sample_colors[3:4]) + xlab('') + ylab('ChIP Enrichment') + coord_cartesian(ylim=c(0,0.8))
p2 <- p2 + theme(legend.position=c(0.2,(1.04)), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 1),nrow = 1))
p <- wrap_plots(p1/p2,guides='collect') & ylab(NULL) & ylab(NULL) & theme(legend.position = 'top',legend.justification = 'center',legend.direction = 'horizontal',legend.box = 'horizontal',plot.margin = margin(5.5, 5.5, 5.5, 0))

pdf('plots/figures/FigureS4L.pdf',width=4,height=8)
print(p)
dev.off()
