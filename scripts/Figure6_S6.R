library(ggplot2)
library(readr)
library(reshape2)
library(ggpubr)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
library(Hmisc)
library(Matrix)
library(SummarizedExperiment)
require(edgeR)
require(LSD)
library(ggrepel)
library(ggpointdensity)
library(plyr)
library(dplyr)
library(seqplots)

set.seed(123)
source('scripts/figures/config.R')
source('scripts/figures/plot_functions.R')
source('scripts/SC/aux_functions.R')
source('scripts/aux_functions.R')

source('/home/hpc/bonev/projects/hic/dfg/config.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/main_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/aux_functions.R')
source('/home/hpc/bonev/projects/hic/dfg/scripts/plot_functions.R')

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

p2glinks <- readRDS('results/P2G-Links.RDS')
tracks <- c('data/Methylation/GFP_methylation_10x_CpG.bw','data/Methylation/Ngn2_methylation_10x_CpG.bw','data/Methylation/PmutNgn2_methylation_10x_CpG.bw')
meth_tracks <- c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x")

#########################

Figure6BC <- function(df_all,sig.mat_all=NULL,plot_density=F,point.size=1,clusters=c('NSC','IPC','PN'),domain='intraTAD',features=NULL,anno.size=12,height,width,cols,out_f,boxplot_only=F){
  df_all$cluster <- c('AST','iN_1','iN_2')[max.col(df_all[,c('AST','iN_1','iN_2')])]
  plot_list <- list()
  features1 <- c()
  for (cluster in clusters){ 
    df <- df_all[grep(cluster,df_all$cluster),]
    if(!is.null(domain)){
      df <- df[df$domain==domain,]
    }
    yscore=ifelse(cluster=='AST','GFP_HiC','Ngn2_HiC')
    if(cluster=='iN_2'){
      yscore='PmutNgn2_HiC'
    }
    df_b <- melt(df[,c('GFP_HiC','Ngn2_HiC','PmutNgn2_HiC')])
    df_b$variable <- gsub('_HiC','',df_b$variable)
    df_b$variable <- factor(df_b$variable,levels=c('GFP','Ngn2','PmutNgn2'))
    p1 <- ggplot(df_b,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p1 <- p1 + scale_fill_manual(values=cols) + xlab('') + ylab('Hi-C Score') + theme(legend.position = "none")
    plot_list[[paste0(cluster,'_sign')]] <- p1 + stat_compare_means(comparisons = list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')),label = "p.format",method='wilcox',paired = T)
    for (xcluster in c('GFP_HiC','Ngn2_HiC','PmutNgn2_HiC')[grep(yscore,c('GFP_HiC','Ngn2_HiC','PmutNgn2_HiC'),invert=T)]){
      xscore <- xcluster
      if (plot_density){
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_pointdensity(shape=19,size=point.size,alpha=1) + scale_color_gradientn(colours = c('black','red','orange'),name='')
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + xlim(0,100) + ylim(0,100) + theme(legend.position='none')
      } else {
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='darkblue') 
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='red') + xlim(0,100) + ylim(0,100)
      }
      df <- df[order(df[,yscore]-df[,xscore],df[,yscore],decreasing=T),]
      if (is.null(features)){
        df1 <- as.data.frame(df %>% group_by(gene_name) %>% dplyr::slice(2))
        df1 <- df1[order(df1[,yscore]-df1[,xscore],df1[,yscore],decreasing=T),]
        df1 <- df1[df1[,yscore]>0&df1[,xscore]>0,]
        features1 <- c(head(df1$labels,5),features1)
      } else {
        features1 <- features
      }
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.85, 0.05),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_blank())+ guides(color = guide_colorbar(barwidth = 3, barheight = 1))
      p <- p + geom_text_repel(
        data = df, size = anno.size,box.padding = 0.5,segment.alpha = 0.5, min.segment.length = 0,max.iter = 20000,max.overlaps = 10000,
        aes(color=NULL,label=ifelse(labels%in%features, as.character(labels), "")),force=10)
      if(!boxplot_only){plot_list[[paste0(cluster,'vs',xcluster)]] <- p}
    }
  }
  p <- Reduce("+", plot_list) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height*(length(clusters)),width=width,useDingbats=FALSE)
  if(!boxplot_only){
    print(p + theme_cowplot()+ plot_layout(nrow=length(clusters),ncol=length(clusters),widths = c(1,rep(3,length(clusters)-1))))
  } else {
    print(p + theme_cowplot()+ plot_layout(nrow=length(clusters),ncol=1,heights = rep(3,length(clusters))))
  }
  dev.off()
  return(features1)
}

Figure6FG <- function(df_all,sig.mat_path=NULL,plot_density=F,point.size=1,clusters=c('NSC','IPC','PN'),domain='intraTAD',features=NULL,anno.size=12,height,width,cols,out_f,boxplot_only=F){
  df_all$cluster <- c('AST','iN_1','iN_2')[max.col(df_all[,c('AST','iN_1','iN_2')])]
  plot_list <- list()
  #features1 <- c("Eomes:-150098")
  for (cluster in clusters){ 
    df <- df_all[grep(cluster,df_all$cluster),]
    if(!is.null(domain)){
      df <- df[df$domain==domain,]
    }
    yscore=ifelse(cluster=='AST','distal.GFP_CpG','distal.Ngn2_CpG')
    if(cluster=='iN_2'){
      yscore='distal.PmutNgn2_CpG'
    }
    df_b <- melt(df[!duplicated(df$peakName),c('distal.GFP_CpG','distal.Ngn2_CpG','distal.PmutNgn2_CpG')])
    df_b$variable <- gsub('distal.','',df_b$variable)
    df_b$variable <- gsub('_CpG','',df_b$variable)
    df_b$variable <- factor(df_b$variable,levels=c('GFP','Ngn2','PmutNgn2'))
    p1 <- ggplot(df_b,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p1 <- p1 + scale_fill_manual(values=cols) + xlab('') + ylab('% Methylation') + theme(legend.position = "none")
    plot_list[[paste0(cluster,'_sign')]] <- p1 + stat_compare_means(comparisons = list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')),label = "p.format",method='wilcox',paired=T)+ scale_y_continuous(breaks = c(0,25,50,75,100),labels = c(0,25,50,75,100))
    for (xcluster in c('distal.GFP_CpG','distal.Ngn2_CpG','distal.PmutNgn2_CpG')[grep(yscore,c('distal.GFP_CpG','distal.Ngn2_CpG','distal.PmutNgn2_CpG'),invert=T)]){
      xscore <- xcluster
      if (plot_density){
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_pointdensity(shape=19,size=point.size,alpha=1) + scale_color_gradientn(colours = c('black','red','orange'),name='')
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + xlim(0,100) + ylim(0,100)
      } else {
        p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='darkblue') 
        p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='red') + xlim(0,100) + ylim(0,100)
      }
      df <- df[order(df[,yscore]-df[,xscore],df[,yscore],decreasing=F),]
      if (is.null(features)){
        df1 <- as.data.frame(df %>% group_by(gene_name) %>% dplyr::slice(2))
        df1 <- df1[order(df1[,yscore]-df1[,xscore],decreasing=F),]
        features1 <- head(df1$labels,20)
      } else {
        features1 <- features
      }
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.85, 0.05),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_blank())+ guides(color = guide_colorbar(barwidth = 3, barheight = 1))
      p <- p + geom_text_repel(
        data = df, size = anno.size,box.padding = 0.8,segment.alpha = 0.5, min.segment.length = 0,max.iter = 20000,max.overlaps = 10000,
        aes(color=NULL,label=ifelse(labels%in%features, as.character(labels), "")),force=10)
      if(!boxplot_only){plot_list[[paste0(cluster,'vs',xcluster)]] <- p1 <- p1 + xlab(xscore) + ylab(yscore)}
    }
  }
  p <- Reduce("+", plot_list) 
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height*(length(clusters)),width=width,useDingbats=FALSE)
  if(!boxplot_only){
    print(p + theme_cowplot()+ plot_layout(nrow=length(clusters),ncol=length(clusters),widths = c(1,rep(3,length(clusters)-1))))
  } else {
    print(p + theme_cowplot()+ plot_layout(nrow=length(clusters),ncol=1,heights = rep(3,length(clusters))))
  }
  dev.off()
}


FigureS6B <- function(df_pos,df_neg,df_no,cluster,yscore,domain='intraTAD',cols){
  df_pos$cluster <- c('AST','iN_1','iN_2')[max.col(df_pos[,c('AST','iN_1','iN_2')])]
  df_pos$type <- 'posCor'
  df_neg$cluster <- c('AST','iN_1','iN_2')[max.col(df_neg[,c('AST','iN_1','iN_2')])]
  df_neg$type <- 'negCor'
  df_no$cluster <- c('AST','iN_1','iN_2')[max.col(df_no[,c('AST','iN_1','iN_2')])]
  df_no$type <- 'noCor'
  df_all <- rbind(df_pos,df_neg,df_no)
  df <- df_all[grep(cluster,df_all$cluster),]
  if(!is.null(domain)){
    df <- df[df$domain==domain,]
  }
  df_b <- melt(df[,c('type',yscore)])
  df_b$type <- factor(df_b$type,levels=c('posCor','negCor','noCor'))
  p1 <- ggplot(df_b,aes(x=type,y=value,fill=type)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =F,width=0.8) 
  p1 <- p1 + scale_fill_manual(name='',values=cols) + xlab('') + ylab('Hi-C Score') + theme(legend.position = "none")
  p1 <- p1 + stat_compare_means(comparisons = list(c('posCor','negCor'),c('posCor','noCor')),label = "p.format",method='wilcox')
  return(p1)
}

FigureS8JK <- function(mpra_annot,sig_mpra,res,df_all,clusters=c('NSC','IPC','PN1','PN2')){
  df <- mpra_annot[mpra_annot$type=='WTposCor'&mpra_annot$cluster%in%clusters,]
  df <- df[df$MPRA_name%in%sig_mpra$MPRA_name,]
  df$MPRA_cluster <- sig_mpra$cluster[match(df$MPRA_name,sig_mpra$MPRA_name)]
  df$MPRA_cluster <- factor(df$MPRA_cluster,levels=c(1:7))
  df <- merge(df,res[,c('MPRA_name','mad.score_NSC','mad.score_IPC','mad.score_PN')],by='MPRA_name',sort=F)
  df$domain <- df_all$domain[match(df$labels,df_all$labels)]
  hic_df <- as.data.frame(df[,c('MPRA_name','MPRA_cluster','NSCscore','IPCscore','PNscore')])
  meth_df <- as.data.frame(df[,c('MPRA_name','MPRA_cluster','distal.E14_NSC_10x','distal.E14_IPC_10x','distal.E14_PN_10x')])
  meth_df <- meth_df[!duplicated(meth_df),]
  
  meth_mat <- melt(meth_df,id.vars = 1:2)
  meth_mat$MPRA_cluster <- factor(meth_mat$MPRA_cluster,levels=rev(levels(meth_mat$MPRA_cluster)))
  meth_mat$variable <- factor(meth_mat$variable,levels=rev(levels(meth_mat$variable)))
  p <- ggplot(meth_mat,aes(x=MPRA_cluster,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8) + coord_flip() + scale_fill_manual(name='',values=rev(c("#F51111","#16A810","#4861F0")))
  p_meth <- p + ylab('%CpG Methylation') +theme(legend.position='none') + geom_hline(yintercept = median(meth_mat$value[meth_mat$MPRA_cluster==7],na.rm=T),col='black',lty=2)
  
  hic_mat <- melt(hic_df,id.vars = 1:2)
  hic_mat$MPRA_cluster <- factor(hic_mat$MPRA_cluster,levels=rev(levels(hic_mat$MPRA_cluster)))
  hic_mat$variable <- factor(hic_mat$variable,levels=rev(levels(hic_mat$variable)))
  p <- ggplot(hic_mat,aes(x=MPRA_cluster,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8) + coord_flip() + scale_fill_manual(name='',values=rev(c("#F51111","#16A810","#4861F0")))
  p_hic <- p + ylab('Hi-C score') +theme(legend.position='none') + geom_hline(yintercept = median(hic_mat$value[hic_mat$MPRA_cluster==7],na.rm=T),col='black',lty=2)
  return(list(p_meth=p_meth,p_hic=p_hic))
}


#### Plot Figures ####

#Figure6B-C
Figure6BC(df_all=read.table('results/P2G_ClusterBinaryMat_posCor.tsv',header=T),out_f='Figure6B',
          clusters=c('AST','iN_1','iN_2'),domain='intraTAD',anno.size=5,plot_density=T,boxplot_only = T,
          features=c('Slc1a3:-70573',"Slc1a3:-105947","Aqp4:-17675","Sox9:206868","Sema5a:-344033","Rbfox3:116442"),point.size=1.5,cols=sample_colors[2:4],
          height=5,width=4)

#Figure6D - Using TSS as anchors
res1 <- getPlotSetArray(tracks=tracks,features="results/beds/P2G_AST_distalSummits.bed",refgenome='BSgenome.Mmusculus.UCSC.mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 50)
res2 <- getPlotSetArray(tracks=tracks,features="results/beds/P2G_iN1_distalSummits.bed",refgenome='BSgenome.Mmusculus.UCSC.mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 50)
res3 <- getPlotSetArray(tracks=tracks,features="results/beds/P2G_iN2_distalSummits.bed",refgenome='BSgenome.Mmusculus.UCSC.mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 50)
pdf('plots/figures/Figure6C_1.pdf',height=6,width=6)
plotAverage(plotset=res1, labels = c('GFP','Ngn2','PmutNgn2'), xlim = NULL,
            ylim = c(18,85), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12,
            colvec = sample_colors[2:4], pointsize = 12, ln.v = FALSE, ln.h = NULL)
dev.off()
pdf('plots/figures/Figure6C_2.pdf',height=6,width=6)
plotAverage(plotset=res2, labels = c('GFP','Ngn2','PmutNgn2'), xlim = NULL,
            ylim = c(18,85), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12,
            colvec = sample_colors[2:4], pointsize = 12, ln.v = FALSE, ln.h = NULL)
dev.off()
pdf('plots/figures/Figure6C_3.pdf',height=6,width=6)
plotAverage(plotset=res3, labels = c('GFP','Ngn2','PmutNgn2'), xlim = NULL,
            ylim = c(18,85), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12,
            colvec = sample_colors[2:4], ln.v = FALSE, ln.h = NULL, pointsize = 12)
axis(side = 1,labels=c('-1KB','0','+1KB'),at=c(-1000,0,1000),pos=18,tick = T)
dev.off()

Figure6FG(df_all=read.table('results/P2G_ClusterBinaryMat_posCor.tsv',header=T),out_f='Figure6D',
          clusters=c('AST','iN_1','iN_2'),domain='intraTAD',anno.size=5,boxplot_only = T,
          features=c("Rbfox3:247336","Neurod1:57029","Rbfox3:496513","Elavl4:73759","Gria1:-455123","Gria1:203029","Robo1:139351","Igfbpl1:-289837","Mdga1:-345147","Mdga1:-358976"),point.size=1.5,cols=sample_colors[2:4],
          height=5,width=4)

source('/home/hpc/bonev/projects/hic/dfg/config.R')
plotMisha(main_f=main_f,targetGene='Mdga1',outDir='plots/figures/',out_f='Figure6G',upstream=4e5,downstream=6e5,
          window_scale=2.2, binSize=5e3,conditions=score_tracks[2:4],
          chipNames=c('AST','iN_1','iN_2'),figure_mode=TRUE,
          chipTracksToExtract=c("scATAC.repro_AST","scATAC.repro_iN_1","scATAC.repro_iN_2"),
          arcIntervals=p2glinks$posCor,arcColors=colorRampPalette(ArchR::ArchRPalettes[[29]]),
          #chipYlim=matrix(c(0,7,0,7,0,7,0,0.8),ncol = 2,byrow = T),
          chipColors=c("#272E6A","#FA7517","#A30005"),img_factor=1.4,
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=TRUE, rna=FALSE, chip=TRUE ,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=100, scores=1.9, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=1.2, arcs=0.6, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))



### Supplementary Figure 6

#Figure6A
pdf('plots/figures/FigureS6A.pdf',width=6.5,height=5)
layout(matrix(c(1:9,10,10,10),nrow=3,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4,4,4),respect = T)
for (cell in cells[2:4]){
  for (cluster in c('AST','iN1','iN2')){
    params <- plot_aggregateHiC(cells=cell,pool=T,intervals1=paste0(cluster,"_posCor_GA"),intervals2=paste0(cluster,"_posCor_DA"),range_f=40000,filter_f=0,res_f=1000,plot_res=8000,grid_mode='merged',zlim=c(-0.5,0.5),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T,plot_mean = F)
  }
}
par(mar=c(7,2,7,2.5))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()


p <- FigureS6B(df_pos=read.table('results/P2G_ClusterBinaryMat_posCor.tsv',header=T),
               df_neg=read.table('results/P2G_ClusterBinaryMat_negCor.tsv',header=T),
               df_no=read.table('results/P2G_ClusterBinaryMat_noCor.tsv',header=T),
               cluster='AST',yscore='GFP_HiC',
               cols=colorRampPalette(c('white','gray30'))(3))
pdf('plots/figures/FigureS6B_1.pdf',height=4,width=3)
print(p+ylab('GFP Hi-C score')+scale_y_continuous(breaks=c(-100,-50,0,50,100)))
dev.off()

p <- FigureS6B(df_pos=read.table('results/P2G_ClusterBinaryMat_posCor.tsv',header=T),
               df_neg=read.table('results/P2G_ClusterBinaryMat_negCor.tsv',header=T),
               df_no=read.table('results/P2G_ClusterBinaryMat_noCor.tsv',header=T),
               cluster='iN_1',yscore='Ngn2_HiC',
               cols=colorRampPalette(c('white','gray30'))(3))
pdf('plots/figures/FigureS6B_2.pdf',height=4,width=3)
print(p+ylab('Ngn2 Hi-C score')+scale_y_continuous(breaks=c(-100,-50,0,50,100)))
dev.off()

p <- FigureS6B(df_pos=read.table('results/P2G_ClusterBinaryMat_posCor.tsv',header=T),
               df_neg=read.table('results/P2G_ClusterBinaryMat_negCor.tsv',header=T),
               df_no=read.table('results/P2G_ClusterBinaryMat_noCor.tsv',header=T),
               cluster='iN_2',yscore='PmutNgn2_HiC',
               cols=colorRampPalette(c('white','gray30'))(3))
pdf('plots/figures/FigureS6B_3.pdf',height=4,width=3)
print(p+ylab('PmutNgn2 Hi-C score')+scale_y_continuous(breaks=c(-100,-50,0,50,100)))
dev.off()

pdf('plots/figures/FigureS8B.pdf',width=6.5,height=5)
layout(matrix(c(1:9,10,10,10),nrow=3,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4,4,4),respect = T)
for (cell in cells[1:3]){
  for (cluster in c('NSC','IPC','PN')){
    params <- plot_aggregateHiC(cells=cell,pool=T,intervals1=paste0(cluster,"_noCor_GA"),intervals2=paste0(cluster,"_noCor_DA"),range_f=60000,filter_f=0,res_f=1000,plot_res=6000,grid_mode='merged',zlim=c(-0.6,0.6),which_plot=c(1),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
  }
}
par(mar=c(7,2,7,2.5))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

Figure6BC(df_all=read.table('results/P2G_binaryMat_noCor.tsv',header=T),out_f='FigureS8CD',
          clusters=c('NSC','IPC','PN'),domain='intraTAD',anno.size=3,plot_density=F,
          features=c(,'Eomes:-150098','Gas2:-94763','Sox5:-243255',"Clstn2:-335963","Gas1:77418"),point.size=1.5,cols=c("#F51111","#16A810","#4861F0"),
          height=4,width=5.5)

Figure6FG(df_all=read.table('results/P2G_binaryMat_noCor.tsv',header=T),out_f='FigureS8EF',
          clusters=c('NSC','IPC','PN'),domain='intraTAD',anno.size=5,plot_density=T,
          features=c("Eomes:-150098","Id4:-391898","Fezf2:-103165","Clstn2:-335963","Gas1:77418"),point.size=1.5,cols=c("#F51111","#16A810","#4861F0"),
          height=4,width=5)

#Figure S8G
res1 <- getPlotSetArray(tracks=tracks,features="results/beds/NSC_noCor_DA.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
res2 <- getPlotSetArray(tracks=tracks,features="results/beds/IPC_noCor_DA.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
res3 <- getPlotSetArray(tracks=tracks,features="results/beds/PN_noCor_DA.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
pdf('plots/figures/FigureS8G_1.pdf',height=6,width=6)
plotAverage(plotset=res1, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12,
            colvec = c("#F51111","#16A810","#4861F0"), pointsize = 12)
dev.off()
pdf('plots/figures/FigureS8G_2.pdf',height=6,width=6)
plotAverage(plotset=res2, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12,
            colvec = c("#F51111","#16A810","#4861F0"), pointsize = 12)
dev.off()
pdf('plots/figures/FigureS8G_3.pdf',height=6,width=6)
plotAverage(plotset=res3, labels = c('NSC','IPC','PN'), xlim = NULL,
            ylim = c(0,100), main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12,
            colvec = c("#F51111","#16A810","#4861F0"), pointsize = 12)
axis(side = 1,labels=c('-2KB','0','+2KB'),at=c(-2000,0,2000),pos=0,tick = T)
dev.off()

Figure6BC(df_all=read.table('results/P2G_binaryMat_negCor.tsv',header=T),out_f='FigureS8H',
          clusters=c('NSC','IPC','PN'),domain='intraTAD',anno.size=3,plot_density=F,boxplot_only = T,
          features=c('Nr2e1:-50936',"Gli3:82631",'Eomes:-150098','Gas2:-94763','Sox5:-243255',"Clstn2:-335963","Gas1:77418"),point.size=1.5,cols=c("#F51111","#16A810","#4861F0"),
          height=4,width=0.8)

Figure6FG(df_all=read.table('results/P2G_binaryMat_negCor.tsv',header=T),out_f='FigureS8I',boxplot_only = T,
          clusters=c('NSC','IPC','PN'),domain='intraTAD',anno.size=3,plot_density=F,
          features=c("Eomes:-150098","Id4:-391898","Fezf2:-103165","Clstn2:-335963","Gas1:77418"),point.size=1.5,cols=c("#F51111","#16A810","#4861F0"),
          height=4,width=0.8)