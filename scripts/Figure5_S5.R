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
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(readr)

library(doParallel)
registerDoParallel(cores=8)


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

meth_tracks <- c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x")
object <- readRDS('results/SC/seurat_IDs.RDS')

Figure5F <- function(score_f,obs_f=NULL,min_dist,max_dist,TAD,labels,cols,p_stats,add_pvalue=T,ylim=NULL){
  if(is.null(obs_f)){
    df <- get(load(score_f))
    if(TAD=='intra'){
      df_e <- cbind(df[[2]]$intra$v_score,df[[3]]$intra$v_score,df[[4]]$intra$v_score)
      rownames(df_e) <- as.numeric(df[[2]]$intra$start2-df[[2]]$intra$start1)
    } else if(TAD=='inter'){
      df_e <- cbind(df[[2]]$inter$v_score,df[[3]]$inter$v_score,df[[4]]$inter$v_score)
      rownames(df_e) <- as.numeric(df[[2]]$inter$start2-df[[2]]$inter$start1)
    }
    colnames(df_e) <- labels
    df_e <- df_e[(as.numeric(rownames(df_e))>=min_dist)&(as.numeric(rownames(df_e))<=max_dist),]
    #df_e <- df_e[complete.cases(df_e),]
    #colnames(df_e) <- names(df)
    df_o <- melt(df_e)
    df_o$Var2 <- factor(df_o$Var2,levels=labels)
    p <- ggplot(df_o,aes(x=Var2,y=value,fill=Var2)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
    p <- p + scale_fill_manual(values=cols) + xlab('') + ylab('Hi-C score') + theme(legend.position = "none")  + ggtitle('')
    if(add_pvalue){
      p <- p + stat_compare_means(comparisons = p_stats,label = "p.format",method = 'wilcox',paired = T,tip.length = c(0.03,0.01)) 
    }
    return(p)
  } else{
    df <- get(load(obs_f))
    dist = 10^seq(4,8,by=0.1)
    min_ind <- findInterval(min_dist,dist)+1
    max_ind <- findInterval(max_dist,dist)
    if(TAD=='intra'){
      df_obs <- cbind(df[[1]][[1]]$obs_intra,df[[1]][[2]]$obs_intra,df[[3]][[1]]$obs_intra,df[[3]][[2]]$obs_intra,df[[3]][[3]]$obs_intra,df[[4]][[1]]$obs_intra,df[[4]][[2]]$obs_intra,df[[5]][[1]]$obs_intra,df[[5]][[2]]$obs_intra)
      df_exp <- cbind(df[[1]][[1]]$exp_intra,df[[1]][[2]]$exp_intra,df[[3]][[1]]$exp_intra,df[[3]][[2]]$exp_intra,df[[3]][[3]]$exp_intra,df[[4]][[1]]$exp_intra,df[[4]][[2]]$exp_intra,df[[5]][[1]]$exp_intra,df[[5]][[2]]$exp_intra)
    } else if(TAD=='inter') {
      df_obs <- cbind(df[[5]][[1]]$obs_inter,df[[5]][[2]]$obs_inter,df[[5]][[3]]$obs_inter,df[[6]][[1]]$obs_inter,df[[6]][[2]]$obs_inter,df[[6]][[3]]$obs_inter,df[[1]][[1]]$obs_inter,df[[1]][[2]]$obs_inter,df[[1]][[3]]$obs_inter,df[[3]][[1]]$obs_inter,df[[3]][[2]]$obs_inter,df[[3]][[3]]$obs_inter,df[[4]][[1]]$obs_inter,df[[4]][[2]]$obs_inter,df[[4]][[3]]$obs_inter)
      df_exp <- cbind(df[[5]][[1]]$exp_inter,df[[5]][[2]]$exp_inter,df[[5]][[3]]$exp_inter,df[[6]][[1]]$exp_inter,df[[6]][[2]]$exp_inter,df[[6]][[3]]$exp_inter,df[[1]][[1]]$exp_inter,df[[1]][[2]]$exp_inter,df[[1]][[3]]$exp_inter,df[[3]][[1]]$exp_inter,df[[3]][[2]]$exp_inter,df[[3]][[3]]$exp_inter,df[[4]][[1]]$exp_inter,df[[4]][[2]]$exp_inter,df[[4]][[3]]$exp_inter)
    } else{
      df_obs1 <- cbind(df[[5]][[1]]$obs_intra,df[[5]][[2]]$obs_intra,df[[5]][[3]]$obs_intra,df[[6]][[1]]$obs_intra,df[[6]][[2]]$obs_intra,df[[6]][[3]]$obs_intra,df[[1]][[1]]$obs_intra,df[[1]][[2]]$obs_intra,df[[1]][[3]]$obs_intra,df[[3]][[1]]$obs_intra,df[[3]][[2]]$obs_intra,df[[3]][[3]]$obs_intra,df[[4]][[1]]$obs_intra,df[[4]][[2]]$obs_intra,df[[4]][[3]]$obs_intra)
      df_exp1 <- cbind(df[[5]][[1]]$exp_intra,df[[5]][[2]]$exp_intra,df[[5]][[3]]$exp_intra,df[[6]][[1]]$exp_intra,df[[6]][[2]]$exp_intra,df[[6]][[3]]$exp_intra,df[[1]][[1]]$exp_intra,df[[1]][[2]]$exp_intra,df[[1]][[3]]$exp_intra,df[[3]][[1]]$exp_intra,df[[3]][[2]]$exp_intra,df[[3]][[3]]$exp_intra,df[[4]][[1]]$exp_intra,df[[4]][[2]]$exp_intra,df[[4]][[3]]$exp_intra)
      df_obs2 <- cbind(df[[5]][[1]]$obs_inter,df[[5]][[2]]$obs_inter,df[[5]][[3]]$obs_inter,df[[6]][[1]]$obs_inter,df[[6]][[2]]$obs_inter,df[[6]][[3]]$obs_inter,df[[1]][[1]]$obs_inter,df[[1]][[2]]$obs_inter,df[[1]][[3]]$obs_inter,df[[3]][[1]]$obs_inter,df[[3]][[2]]$obs_inter,df[[3]][[3]]$obs_inter,df[[4]][[1]]$obs_inter,df[[4]][[2]]$obs_inter,df[[4]][[3]]$obs_inter)
      df_exp2 <- cbind(df[[5]][[1]]$exp_inter,df[[5]][[2]]$exp_inter,df[[5]][[3]]$exp_inter,df[[6]][[1]]$exp_inter,df[[6]][[2]]$exp_inter,df[[6]][[3]]$exp_inter,df[[1]][[1]]$exp_inter,df[[1]][[2]]$exp_inter,df[[1]][[3]]$exp_inter,df[[3]][[1]]$exp_inter,df[[3]][[2]]$exp_inter,df[[3]][[3]]$exp_inter,df[[4]][[1]]$exp_inter,df[[4]][[2]]$exp_inter,df[[4]][[3]]$exp_inter)
      df_obs <- df_obs1+df_obs2
      df_exp <- df_exp1+df_exp2
    }
    res <- data.frame(condition=c(rep(labels[1],3),rep(labels[2:4],each=2)),rep=paste0(rep('rep_',9),c(1:3,1:2,1:2)),obs_exp=log2(colSums(df_obs[min_ind:max_ind,],na.rm=T)/colSums(df_exp[min_ind:max_ind,],na.rm=T)))
    res$condition <- factor(res$condition,levels=labels)
    df_summ <- data_summary(res, varname="obs_exp",groupnames=c('condition'))
    p <- ggplot(df_summ, aes(x=condition, y=obs_exp, fill=condition)) + scale_fill_manual(values = cols) +
      geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('log2(obs/exp)') +
      geom_errorbar(aes(ymin=obs_exp-sd, ymax=obs_exp+sd), width=.2,position=position_dodge(.9))
    p <- p + geom_point(data = res,aes(x=condition, y=obs_exp, fill=condition,shape=rep),position = position_jitterdodge(dodge.width =.2,jitter.width = 0.1,seed = 42),size=2) 
    if(add_pvalue){
      p <- p + stat_compare_means(data = res,comparisons = list(p_stats),label = "p.format",method = 't.test',paired = T,label.y = df_summ$obs_exp[2]*1.05,tip.length = c(0.05,0.01)) 
    }
  }
}


Figure5H <- function(object,features,cols,point.size,out_f,height,width,plot_filled=F,anno.size=10,theme=NULL,direction='vertical'){
  plot_list <- list()
  for (feature in features){
    plot.data <- extractFeatures(object,features=feature,cells_toInclude = c('all'),cells_toExclude = 'none')
    if(!plot_filled){
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(color=cols[1],size = point.size,data = plot.data[plot.data$feature==0,]) + geom_point(aes(color=feature),size = point.size,data = plot.data[plot.data$feature>0,]) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_color_gradientn(colours=cols,name='',breaks=c(0,max(plot.data$feature)),labels=c(0,round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1, 1),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 8))
      p <- p + guides(color = guide_colourbar(barwidth = 3, barheight = 1))
    } else {
      p <- ggplot(plot.data, aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(fill = feature), pch = I(21),size = point.size,stroke=0.1) + xlab("UMAP1") + ylab("UMAP2") 
      p <- p + scale_fill_gradientn(colours=cols,name='',breaks=c(min(plot.data$feature,na.rm=T),max(plot.data$feature,na.rm=T)),labels=c(min(plot.data$feature,na.rm=T),round(max(plot.data$feature))))
      p <- p + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(1.05, 0.98),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_text(hjust = c(1.5,-0.5), vjust = 7.5))
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



FigureS7F <- function(cor_file,cells=c('NSC','IPC','PN'),cols=colorRampPalette(rev(colorpalette('ylorrd')))(100),cluster_cols=c("#F51111","#16A810","#4861F0"),out_f,width,height){
  res <- read.table(cor_file)
  row.names(res) <- paste0(rep(cells,each=3),'_',1:3)
  colnames(res) <- row.names(res)
  distance <- dist(1-res)
  res <- round(res, 2)
  annotation_df <- data.frame(cell_Type=factor(rep(cells,each=3),levels=cells))
  colnames(annotation_df) <- c('Condition')
  row.names(annotation_df) <- colnames(res)
  pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
  pheatmap(res,color=cols,breaks=seq(0.85,1,length=101),clustering_distance_rows=distance,clustering_distance_cols=distance,clustering_method='ward.D2',show_colnames = F,
           display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",legend_breaks=c(0.85,0.9,0.95,1))
  dev.off()
}

FigureS5B<-function(meth_tracks,height,width,out_f,cluster_cols,peaks,window_size) {
  res <- extract_lin(regions=peaks,window=window_size,tracks=meth_tracks)
  colnames(res)[4:12] <- c(paste0('GFP rep',1:3),paste0('Ngn2 rep',1:2),paste0('PmutNgn2rep',1:2))
  res<-na.omit(res)
  df<-as.matrix(res[,4:12])
  res<-cor(df,method='pearson')
  distance <- dist(1-res)
  res <- round(res, 2)
  cells<-c('NSC','IPC','PN')
  annotation_df <- data.frame(cell_Type=factor(rep(cells,each=3),levels=cells))
  colnames(annotation_df) <- c('Condition')
  row.names(annotation_df) <- colnames(res)
  cols=colorRampPalette(rev(colorpalette('ylorrd')))(100)
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(pheatmap(res,color=cols,breaks=seq(0.85,1,length=101),clustering_distance_rows=distance,clustering_distance_cols=distance,clustering_method='ward.D',show_colnames = F,annotation_colors=list(Condition=c(NSC=cluster_cols[1],IPC=cluster_cols[2],PN=cluster_cols[3])),
                 angle_col = 0,display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.85,0.9,0.95,1))
  )
  dev.off()
}

FigureS5B <- function(cells,levels_names,out_f,height=4,width=3){
  res_comp <- comp_boundaries(type='eigen',cells=cells, write_borders=F)
  res_comp$condition <- gsub('JD_','',row.names(res_comp))
  res_comp$condition <- gsub('_D2_G0G1','',res_comp$condition)
  df <- melt(res_comp)
  df$condition <- factor(df$condition,levels=levels_names)
  df$variable <- factor(df$variable,levels=c('A-B','B-A'))
  p <- ggplot(df, aes(x=condition,y=value,fill=variable)) +  geom_bar(position="stack", stat="identity",colour="black") + scale_fill_manual(name='',values = c("#F0F0F0","#525252"))
  p <- p + ylab('Compartment transitions') + xlab('') + theme(legend.position=c(0.75,0.95))
  pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}

FigureS7IJ <- function(borders,tss,extra_tracks,ins_tracks,features=NULL,width,height,plot_what='ins',plot_cells=c('NSC','PN'),point.size=1,anno.size=4,out_f,mode='scatter',cols,ylab_n){
  for (track in c(ins_tracks,extra_tracks)){
    gvtrack.create(paste0('v_',track),track,'avg')
  }
  if(!is.null(extra_tracks)){
    df_borders <- gextract(c(paste0('v_',ins_tracks,'*(-1)'),paste0('v_',extra_tracks)),intervals = borders,iterator=borders)
  } else {
    df_borders <- gextract(c(paste0('v_',ins_tracks,'*(-1)')),intervals = borders,iterator=borders)
  }
  
  df_borders$labels <- paste0(df_borders$chrom,':',df_borders$start,'-',df_borders$end)
  df_borders$genes <- gintervals.neighbors(df_borders,tss)$geneName
  if(mode=='scatter'){
    xscore <- colnames(df_borders)[intersect(grep(plot_cells[1],colnames(df_borders)),grep(plot_what,colnames(df_borders)))]
    yscore <- colnames(df_borders)[intersect(grep(plot_cells[2],colnames(df_borders)),grep(plot_what,colnames(df_borders)))]
    p <- ggplot(df_borders, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_point(size=point.size,alpha=1,col='darkblue') 
    p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='red') + xlab(plot_cells[1]) + ylab(plot_cells[2])
    p <- p + geom_text_repel(
      data = df_borders, size = anno.size,box.padding = 0.8, min.segment.length = 0,max.iter = 10000,
      aes(label=ifelse(genes%in%features, as.character(genes), "")),force=10)
    pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
    print(p)
    dev.off()
  } else{
    colnames(df_borders)[4:(3+length(plot_cells))] <- plot_cells
    df_borders <- df_borders[,c('labels',plot_cells)]
    df <- melt(df_borders,id.vars = 'labels')
    p <- ggplot(df,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
    p <- p + scale_fill_manual(values=cols) + xlab('') + ylab(ylab_n) + theme(legend.position = "none")
    return(p)
  }
}


FigureS5D <- function(mat_f,col.names,samples=sample_names[2:4],cols=sample_colors[2:4],ylab='Number of TADs',out_f,height,width){
  mat <- read.table(mat_f,header=T)
  colnames(mat) <- col.names
  mat$condition <- samples
  mat_m <- melt(mat)
  df_summ <- data_summary(mat_m, varname="value",groupnames=c('condition'))
  p <- ggplot(df_summ, aes(x=condition, y=value, fill=condition)) + scale_fill_manual(name='',values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab(ylab) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))
  p <- p+ geom_point(data = mat_m,aes(x=condition, y=value, fill=condition),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2) + theme(legend.position = 'none')
  pdf(paste0('plots/figures/',out_f),height=height,width=width)
  print(p)
  dev.off()
}

FigureS5E <- function(mat_fs,samples=sample_names[2:4],FDR.cutoff=0.1,cols=sample_colors[2:4],ylab='Number of Loops',out_f,height,width,returnData=FALSE){
  loop_list <- list()
  for (i in 1:length(mat_fs)){
    mat <- read.table(mat_fs[i],header=T)
    mat <- mat[!duplicated(mat),]
    mat$SAMPLE <- samples[i]
    loop_list[[i]] <- mat
  }
  df <- do.call("rbind", loop_list)
  df[,1:8] <- sapply(df[,1:8],as.numeric)
  df <- df[df$FDR<=FDR.cutoff,]
  df <- df[!rowSums(is.na(df))>=1,]
  mat <- as.data.frame(table(df$SAMPLE))
  colnames(mat) <- c('condition','loops')
  p <- ggplot(mat, aes(x=condition, y=loops, fill=condition)) + scale_fill_manual(name='',values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab(ylab)
  pdf(paste0('plots/figures/',out_f),height=height,width=width)
  print(p+theme(legend.position = 'none'))
  dev.off()
  if(returnData){
    return(df)
  }
}

FigureS5F <- function(res,peaks_fs,conditions=sample_names[2:4],max_dist=1000,out_f,height,width){
  res[,1] <- paste0('chr',res[,1])
  res[,4] <- paste0('chr',res[,4])
  left_a <- res[,c(1:3,9)]
  colnames(left_a) <- c('chrom','start','end','condition')
  right_a <- res[,c(4:6,9)]
  colnames(right_a) <- c('chrom','start','end','condition')
  peaks1 <- rbind(left_a,right_a)
  peaks1 <- peaks1[!duplicated(peaks1),]
  colnames(peaks1) <- c('chrom','start','end','loops_direction')
  peaks1 <- peaks1[complete.cases(peaks1),]
  df_list <- list()
  for (i in 1:length(conditions)){
    peaks1_s <- peaks1[peaks1$loops_direction==conditions[i],]
    peaks2 <- read.table(peaks_fs[i],header=F)[,c(1:3)]
    colnames(peaks2) <- c('chrom','start','end')
    peaks2 <- peaks2[!duplicated(peaks2),]
    df_s <- gintervals.neighbors(peaks1_s,peaks2)
    df_s$peak_direction <- 'Yes'
    df_s$peak_direction[abs(df_s$dist)>max_dist] <- 'No'
    df_s <- df_s[,c('loops_direction','peak_direction')]
    df_list[[i]] <- df_s
  }
  df <- do.call("rbind", df_list)
  
  df$peak_direction <- factor(df$peak_direction,levels=c('Yes','No'))
  df$loops_direction <- factor(df$loops_direction,levels=c('GFP','Ngn2','PmutNgn2'))
  p <- ggplot(df, aes(x=loops_direction,y=1,fill=peak_direction)) +  geom_bar(position="fill", stat="identity") + scale_fill_grey(name='',start = 0.2,end=0.9) 
  p <- p + ylab('Percentage of Loops') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')
  p <- p + theme(legend.position=c(0.5,1.05), legend.box = "horizontal",plot.margin = margin(0.25,0,0,0,unit="inch")) + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
  pdf(paste0('plots/figures/',out_f),height=height,width=width)
  print(p)
  dev.off()
}

FigureS5F <- function(res,peaks_f,max_dist=1000){
  res[,1] <- paste0('chr',res[,1])
  res[,4] <- paste0('chr',res[,4])
  left_a <- res[,c(1:3,9)]
  colnames(left_a) <- c('chrom','start','end','condition')
  right_a <- res[,c(4:6,9)]
  colnames(right_a) <- c('chrom','start','end','condition')
  peaks1 <- rbind(left_a,right_a)
  peaks1 <- peaks1[!duplicated(peaks1),]
  peaks2 <- read.table(peaks_f,header=T)[,c(1:3,11)]
  colnames(peaks1) <- c('chrom','start','end','loops_direction')
  colnames(peaks2) <- c('chrom','start','end','peak_direction')
  peaks1 <- peaks1[complete.cases(peaks1[,2:3]),]
  peaks2 <- peaks2[complete.cases(peaks2),]
  df <- gintervals.neighbors(peaks1,peaks2)
  df$peak_direction[abs(df$dist)>max_dist] <- 'NotBound'
  df$peak_direction <- factor(df$peak_direction,levels=c('Ngn2','shared','PmutNgn2','NotBound'))
  df$loops_direction <- factor(df$loops_direction,levels=c('GFP','Ngn2','PmutNgn2'))
  df <- df[,c('loops_direction','peak_direction')]
  p <- ggplot(df, aes(x=loops_direction,y=1,fill=peak_direction)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='',values=c("#003C30","#CCCCCC","#543005","#404040"))
  p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')
  p <- p + theme(legend.position=c(0.02,1.05), legend.box = "horizontal") + guides(fill=guide_legend(override.aes = list(size = 5),nrow = 1))
  pdf(paste0('plots/figures/',out_f),height=height,width=width)
  print(p)
  dev.off()
}



FigureS5K <- function(input,cols=ATAC_cluster_colors[1:3],width,height){
  input$labels <- paste0(input$chrom,':',input$start,'-',input$end)
  df_ins <- melt(input[,c(grep('labels',colnames(input)),grep('ins',colnames(input)),grep('cluster',colnames(input)))],id.vars = c('labels','cluster'))
  df_ins$variable <- gsub('v_hic.E14_','',df_ins$variable)
  df_ins$variable <- gsub('.ins_250','',df_ins$variable)
  df_ins$variable <- factor(df_ins$variable,levels=c('NSC','IPC','PN'))
  df_meth <- melt(input[,c(grep('labels',colnames(input)),grep('meth',colnames(input)),grep('cluster',colnames(input)))],id.vars = c('labels','cluster'))
  df_meth$variable <- gsub('v_methylation.E14_','',df_meth$variable)
  df_meth$variable <- gsub('_10x','',df_meth$variable)
  df_meth$variable <- factor(df_meth$variable,levels=c('NSC','IPC','PN'))
  df_expr <- melt(input[abs(input$dist)<=5e4,c(grep('labels',colnames(input)),which(colnames(input)%in%c('NSC','IPC','PN')),grep('cluster',colnames(input)))],id.vars = c('labels','cluster'))
  df_expr$variable <- factor(df_expr$variable,levels=c('NSC','IPC','PN'))
  plot_list <- list()
  for (i in 1:max(input$cluster)){
    p1 <- ggplot(subset(df_ins,cluster==i),aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p1 <- p1 + scale_fill_manual(values=cols) + xlab('')  + scale_y_continuous(n.breaks =4,limits = c(1.75,3.5)) + ylab(ifelse(i==1,'Insulation Score','')) + theme(legend.position = "none")
    plot_list[[paste0('ins_',i)]] <- p1 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.signif",method='wilcox',method.args = list(paired=TRUE)) + ggtitle(paste0('Cluster',i))
    p2 <- ggplot(subset(df_meth,cluster==i),aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    p2 <- p2 + scale_fill_manual(values=cols) + xlab('') + scale_y_continuous(breaks =c(0,50,100),limits = c(0,120)) + ylab(ifelse(i==1,'% CpG Methylation','')) + theme(legend.position = "none")
    plot_list[[paste0('meth_',i)]] <- p2 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.signif",method='wilcox',method.args = list(paired=TRUE)) + ggtitle(' ')
    #  p3 <- ggplot(subset(df_expr,cluster==i),aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8) 
    #   p3 <- p3 + scale_fill_manual(values=cols) + xlab('') + ylab('log2(FPKM+1)') + theme(legend.position = "none")
    #  plot_list[[paste0('expr_',i)]] <- p3 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','PN'),c('IPC','PN')),label = "p.signif",method='wilcox',method.args = list(paired=TRUE)) + ggtitle(' ')
  }
  p <- Reduce("+", plot_list) 
  pdf('plots/figures/figureS5K.pdf',width=width*(length(unique(input$cluster))),height=height*2)
  print(p + theme_cowplot()+ plot_layout(ncol=length(unique(input$cluster)),nrow=2,byrow = F))
  dev.off()
  
}


#### Plot Figures ######
#Figure 5A
run_cisDecay(tracks=all_tracks,out_f='plots/figures/Figure5E.pdf',cells=cells[2:5],path=paste0(main_f,'analysis/cis_decay/'),log_base=2,file_f='cisDecay_080122',minDist = 3e3,maxDist=1e8,colors=c(sample_colors[2:4],'black'),alpha=0.3)

#Figure 5B
plotBinned(extra_tracks=meth_tracks[1],fig_name='plots/figures/Figure5B_1.pdf',cells="JD_GFP_D2_G0G1",init_cell="JD_GFP_D2_G0G1",chr='chr3',binSize=2e5,path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),plot_what='obs',balance=T,width=4,height=4.5)
plotBinned(extra_tracks=meth_tracks[2],fig_name='plots/figures/Figure5B_2.pdf',cells="JD_Ngn2_D2_G0G1",init_cell="JD_GFP_D2_G0G1",chr='chr3',binSize=2e5,path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),plot_what='obs',balance=T,width=4,height=4.5)
plotBinned(extra_tracks=meth_tracks[3],fig_name='plots/figures/Figure5B_3.pdf',cells="JD_pmutNgn2_D2_G0G1",init_cell="JD_GFP_D2_G0G1",chr='chr3',binSize=2e5,path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),plot_what='obs',balance=T,width=4,height=4.5)
pdf('plots/figures/Figure5B_scaleBar.pdf',width=0.75,height=4.5)
par(mar=c(0.5,0.5,0.5,1.5))
image.scale(as.matrix(1),zlim=c(0,0.001568), col=colorRampPalette(c("white","orange","red","darkRed"))(1000),axis.pos=4,adj=1,cex.axis=1,axis_lab=c('','',''))
dev.off()

#Figure 5C
zlim=c(-1.5,1.2)
blue_white_pal = colorRampPalette(c("purple", "navy", "blue", "#87FFFF", "white"))
white_red_pal = colorRampPalette(c("white","#FF413D", "black", "orange", "yellow"))
averageTAD_colors <- c(blue_white_pal(length(seq(zlim[1],0.01,by=0.01))),'white',white_red_pal(length(seq(0.01,zlim[2],by=0.01))))
res1 <- averageTrack(tracks="methylation.JD_GFP_D2_G0G1_CpG_10x",regions="hic.JD_GFP_D2_G0G1.ins_250_domains_expanded",bins=100,anchor_middle=F)
res2 <- averageTrack(tracks="methylation.JD_Ngn2_D2_G0G1_CpG_10x",regions="hic.JD_Ngn2_D2_G0G1.ins_250_domains_expanded",bins=100,anchor_middle=F)
res3 <- averageTrack(tracks="methylation.JD_pmutNgn2_D2_G0G1_CpG_10x",regions="hic.JD_pmutNgn2_D2_G0G1.ins_250_domains_expanded",bins=100,anchor_middle=F)
plot_averageTAD(tracks=all_tracks,fig_name='plots/figures/Figure5C.pdf',add_matrix=list(JD_GFP_D2_G0G1=res1,JD_Ngn2_D2_G0G1=res2,JD_pmutNgn2_D2_G0G1=res3),mat_lim=c(60,90),cells=cells[2:4],path='/home/hpc/bonev/projects/hic/dfg/analysis/averageTAD/',file_f='averageTAD_D2only',stats_f=F,plot_what='',zlim=zlim,flip=T,z_colors=averageTAD_colors,height=2.5,width=4.5)

#Figure 5D
plot_rankMatrix(file_f=paste0(main_f,'analysis/compartments/','rankMatrix_250kb_ranks',100),out_f='plots/figures/Figure5D_1.pdf',zlim=c(-1,1),cells=cells[2],col=wide_red_blue_pal(1000),plot_chip=FALSE,plot_scale=F)
plot_rankMatrix(file_f=paste0(main_f,'analysis/compartments/','rankMatrix_250kb_ranks',100),out_f='plots/figures/Figure5D_2.pdf',zlim=c(-1,1),cells=cells[3],col=wide_red_blue_pal(1000),plot_chip=FALSE,plot_scale=F)
plot_rankMatrix(file_f=paste0(main_f,'analysis/compartments/','rankMatrix_250kb_ranks',100),out_f='plots/figures/Figure5D_3.pdf',zlim=c(-1,1),cells=cells[4],col=wide_red_blue_pal(1000),plot_chip=FALSE,plot_scale=T)
#Figure 5E


pdf('plots/figures/Figure5E_new.pdf',width=7,height=4)
layout(matrix(c(1:4),nrow=1,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4),respect = T)
for (cell in cells[2:4]){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='data/ChIP/Ngn2_peaks5000.bed',intervals2='data/ChIP/Ngn2_peaks5000.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(2.5,1.5,2.5,4))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#Figure 5F
p <- Figure5F(score_f='/home/hpc/bonev/projects/hic/dfg/data/cis_decay/Ngn2_peaks5000.bed_Ngn2_peaks5000.bed.1D.5000_score',min_dist = 5e4,max_dist = 2e6,TAD='intra',labels=c('GFP','Ngn2','PmutNgn2'),cols=sample_colors[2:4],p_stats=list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')))
pdf('plots/figures/Figure5F.pdf',height=5,width=4)
print(p)
dev.off()
#Figure 5G
require(seqplots)
meth_bw <- c('data/Methylation/GFP_methylation_10x_CpG.bw','data/Methylation/Ngn2_methylation_10x_CpG.bw','data/Methylation/PmutNgn2_methylation_10x_CpG.bw')
res <- getPlotSetArray(tracks=meth_bw,features="data/ChIP/Ngn2_peaks5000_distal.bed",refgenome='mm10',type = 'mf',add_heatmap=T,xmin=1000,xmax=1000,bin = 50)
pdf('plots/figures/Figure5G.pdf',height=5,width=5)
plotAverage(plotset=res, labels = c('GFP','Ngn2','PmutNgn2'), xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = '%CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = sample_colors[2:4], pointsize = 12)
axis(side = 1,labels=c('-1000','-500','0','+500','+1000'),at=c(-1000,-500,0,500,1000),pos=21,tick = T)
dev.off()



df <- read.table("data/ChIP/Ngn2_peaks5000_distal.bed",header=F,sep='\t')
df <- gintervals(df$V1,df$V2,df$V3)

mat <- misha_extract(tracks=meth_tracks,track_names = c('GFP','Ngn2','PmutNgn2'),regions=df,window = 100,iterator=df,mode='avg')
df_o <- melt(mat[,4:6])
df_o$variable <- factor(df_o$variable,levels=c('GFP','Ngn2','PmutNgn2'))
p <- ggplot(df_o,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
p <- p + scale_fill_manual(values=sample_colors[2:4]) + xlab('') + ylab('% CpG Methylation') + theme(legend.position = "none")  + ggtitle('')
if(add_pvalue){
  p <- p + stat_compare_means(comparisons = list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')),label = "p.format",method = 'wilcox',paired = T,tip.length = c(0.03,0.01)) 
}

pdf('plots/figures/Figure5H.pdf',height=5,width=4)
print(p + scale_y_continuous(breaks = c(0,25,50,75,100),labels = c(0,25,50,75,100)))
dev.off()

source('/home/hpc/bonev/projects/SC/Noack_NatNeuro2021/scripts/figures/plot_functions.R')
plotMisha(main_f=main_f,targetGene='chr5,125128565,125129645',outDir='plots/figures/',out_f='Figure5I',upstream=200,downstream=200,
          chipNames=c('','','',''),plot.dgram=F,score_tracks=score_tracks[2:4],
          chipTracksToExtract=c("scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2","chipseq_RPM.iN_Ngn2_D2"),
          chipColors=c(sample_colors,'black'),imgPlotScale=3,
          methTracksToExtract=c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x"),methColors=c("#F51111","#16A810","#4861F0"),methNames=c('','',''),
          plotOrder=list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=TRUE,MPRA=FALSE,axis=FALSE,scATAC=FALSE, genes=FALSE,ideogram=FALSE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.6,MPRA=0.2,meth=0.6, domains=0.15, genes=0.7, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

#Figure 5G
plotMisha(object=atac.object,targetGene='Flrt2',out_f='Figure5G',upstream=4e6,downstream=2.5e6,
          chipNames=c('NSC','IPC','PN'),window_scale=2,pointCEX=0.5,chipRes =200,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC','scATAC.E14_PN2'),
          chipColors=ATAC_cluster_colors[c(1,2,4)],scoreTrackToExtract =score_tracks[1:3],
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))
pdf('plots/figures/Figure5G_scaleBar.pdf',height=6,width=1.2)
par(mar=c(0.5,0.5,0.5,4))
image.scale(as.matrix(1),zlim=c(-100,100), col=col.scores[1:201],axis.pos=4,adj=1,cex.axis = 1.2)
dev.off()

#Figure 5H
Figure5H(object=rna.object,c('Gas1','Flrt2'),cols=gene_colors(50),point.size=2,out_f='Figure5H',width=5,height=5,plot_filled=T,direction='vertical',anno.size=8,theme=theme_border)

### Supplementary Figure 5

FigureS5D(mat_f="results/HiC/Domain_number.txt",col.names=paste0('rep',1:3),samples=sample_names[2:4],cols=sample_colors[2:4],ylab='Number of TADs',out_f='FigureS5D_new.pdf',height=4,width=3.5)
res <- FigureS5E(mat_fs=c("results/HiC/D2_GFP_loops.tsv","results/HiC/D2_Ngn2_loops.tsv","results/HiC/D2_PmutNgn2_loops.tsv"),
                 samples=sample_names[2:4],FDR.cutoff=0.1,cols=sample_colors[2:4],ylab='Number of Loops',out_f='FigureS5E_newBla.pdf',height=4,width=3.5,returnData=T)
FigureS5F(res=res,peaks_fs=c('CutRun/full/wt_GFP_Cohesin/merge/wt_GFP_Rad21_peaks.narrowPeak','CutRun/full/wt_ngn2_Cohesin//merge/wt_Ngn2_Rad21_peaks.narrowPeak','CutRun/full/wt_pmut_ngn2_Cohesin/rep1/peakcalling/macs2.narrow/wt_pmut_ngn2_cohesin_peaks.narrowPeak'),
          conditions=sample_names[2:4],max_dist=10000,out_f='FigureS5F_new.pdf',height=4,width=3.5)
FigureS5F(res=res,peaks_fs=c('ChIP_Nov22/Data/Peaks/Ngn2_conservative_peaks.narrowPeak','ChIP_Nov22/Data/Peaks/Ngn2_conservative_peaks.narrowPeak','ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak'),
          conditions=sample_names[2:4],max_dist=10000,out_f='FigureS5F_new.pdf',height=4,width=3.5)

  
FigureS5B(cells=cells[2:4],levels_names=c('GFP','Ngn2','pmutNgn2'),out_f='FigureS5B',height=4,width=3)

library(seqplots)
tracks <- c('data/Methylation/GFP_methylation_10x_CpG.bw','data/Methylation/Ngn2_methylation_10x_CpG.bw','data/Methylation/PmutNgn2_methylation_10x_CpG.bw')
res_meth <- getPlotSetArray(tracks=tracks,features="results/beds/SC_allATAC_CTCF.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10)
pdf('plots/figures/FigureS5C.pdf',height=5,width=5)
plotAverage(plotset=res_meth, labels = c('GFP','Ngn2','PmutNgn2'), xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = sample_colors[2:4], pointsize = 12)
axis(side = 1,labels=c('-1kb','CTCF->','+1kb'),at=c(-1000,0,1000),pos=10,tick = T)
dev.off()

res_meth <- getPlotSetArray(tracks=tracks,features="results/",refgenome='mm10',type = 'af',add_heatmap=F,xmin=2500,xmax=2500,xanchored = 5000,bin = 10)
pdf('plots/figures/FigureS5D.pdf',height=5,width=6)
plotAverage(plotset=res_meth, labels = c('GFP','Ngn2','PmutNgn2'), xlim = NULL,
            ylim = NULL, main = NULL, xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',
            cex.main = 20, cex.legend = 10, ln.v = FALSE, ln.h = NULL,
            colvec = sample_colors[2:4], pointsize = 12)
axis(side = 1,labels=c('-2.5kb','TSS','TTS','+2.5kb'),at=c(-2000,0,5000,7000),pos=9.8,tick = T)
dev.off()

#Figure 5E
pdf('plots/figures/FigureS5E.pdf',width=7,height=4)
layout(matrix(c(1:4),nrow=1,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4),respect = T)
for (cell in cells[2:4]){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='data/ChIP/conservative_peaks.narrowPeak',intervals2='data/ChIP/conservative_peaks.narrowPeak',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(2.5,1.5,2.5,4))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()
#Figure S5F
p <- Figure5F(score_f='/home/hpc/bonev/projects/hic/dfg/data/cis_decay/conservative_peaks.narrowPeak_conservative_peaks.narrowPeak.1D.5000_score',min_dist = 5e4,max_dist = 2e6,TAD='intra',labels=c('GFP','Ngn2','PmutNgn2'),cols=sample_colors[2:4],p_stats=list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')))
pdf('plots/figures/FigureS5F.pdf',height=5,width=4)
print(p)
dev.off()
#FigureS5G
pdf('plots/figures/FigureS5G_new.pdf',width=7,height=4)
layout(matrix(c(1:4),nrow=1,ncol=4,byrow=F),widths = c(4,4,4,2),heights=c(4),respect = T)
for (cell in cells[2:4]){
  params <- plot_aggregateHiC(cells=cell,pool=T,intervals1='data/ChIP/PmutNgn2/PmutNgn2_top5000.bed',intervals2='data/ChIP/PmutNgn2/PmutNgn2_top5000.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),plot_scale = F,interval1_name = '',interval2_name = '',add_plot=T)
}
par(mar=c(2.5,1.5,2.5,4))
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

#Figure S5H
p <- Figure8J(score_f='/home/hpc/bonev/projects/hic/dfg/data/cis_decay/PmutNgn2_top5000.bed_PmutNgn2_top5000.bed.1D.5000_score',min_dist = 5e4,max_dist = 2e6,TAD='intra',labels=c('GFP','Ngn2','PmutNgn2'),cols=sample_colors[2:4],p_stats=list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')))
pdf('plots/figures/FigureS5H_new.pdf',height=5,width=4)
print(p)
dev.off()

df <- read.table("data/ChIP/PmutNgn2/PmutNgn2_top5000_Distal_matched.bed",header=F,sep='\t')
df <- gintervals(df$V1,df$V2,df$V3)
mat <- misha_extract(tracks=meth_tracks,track_names = c('GFP','Ngn2','PmutNgn2'),regions=df,window = 0,iterator=df,mode='avg')
df_o <- melt(mat[,4:6])
df_o$variable <- factor(df_o$variable,levels=c('GFP','Ngn2','PmutNgn2'))
p <- ggplot(df_o,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8)
p <- p + scale_fill_manual(values=sample_colors[2:4]) + xlab('') + ylab('% CpG Methylation') + theme(legend.position = "none")  + ggtitle('')
p <- p + stat_compare_means(comparisons = list(c('GFP','Ngn2'),c('Ngn2','PmutNgn2')),label = "p.format",method = 'wilcox',paired = T,tip.length = c(0.03,0.01)) 
pdf('plots/figures/FigureS5I_new.pdf',height=5,width=4)
print(p + scale_y_continuous(breaks = c(0,25,50,75,100),labels = c(0,25,50,75,100)))
dev.off()


#Figure S5L
df <- read.table("data/ChIP/Ngn2_peaks5000_distal_Ngn2Motif.bed",header=F,sep='\t')
df <- gintervals(df$V1,df$V2,df$V3)

mat <- misha_extract(tracks=meth_tracks,track_names = c('GFP','Ngn2','PmutNgn2'),regions=df,window = 50,iterator=df,mode='avg')


p <- ggplot(data = mat, aes(GFP)) + 
  geom_histogram(bins = 10,fill=sample_colors[2],color='black',position = position_dodge(10),
                 aes(y = stat(width*density*100))) + ylab("Percentage of all Ngn2 peaks") + xlab("% Methylation in GFP")
pdf('plots/figures/FigureS5L.pdf',height=5,width=4)
p
dev.off()
#Figure S5M
p <- ggplot(mat, aes( x = GFP, y = Ngn2 )) + geom_pointdensity(shape=19,size=1.5,alpha=1) + scale_color_gradientn(colours = c('darkblue','red','orange'),name='')
p <- p + geom_abline(slope = 1,size=2,intercept = 0,linetype = 2,col='black',) + xlim(0,100) + ylim(0,100) + theme(legend.position = c(0.15,0.85))
pdf('plots/figures/FigureS5M.pdf',height=5,width=5)
p
dev.off()


hic_files <- c('/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/GFP_D2_G0G1/rep1/aligned/inter_30.hic',
               '/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/GFP_D2_G0G1/rep2/aligned/inter_30.hic',
               '/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/GFP_D2_G0G1/rep3/aligned/inter_30.hic',
               '/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/Ngn2_D2_G0G1/rep1/aligned/inter_30.hic',
               '/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/Ngn2_D2_G0G1/rep2/aligned/inter_30.hic',
               '/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/PmutNgn2_D2_G0G1/rep1/aligned/inter_30.hic',
               '/home/jdiwakar/MethylHiC/Day2_Ngn2_RP/juicer/PmutNgn2_D2_G0G1/rep2/aligned/inter_30.hic')
res_df <- correlate_hic(hic_files,5e4)            #in aux_functions.R
dimnames(res_df) <- list(c('GFP_rep1','GFP_rep2','GFP_rep3','Ngn2_rep1','Ngn2_rep2','PmutNgn2_rep1','PmutNgn2_rep2'))
cols=colorRampPalette(rev(colorpalette('ylorrd')))(100)
res_df <- read.table('results/hicRep_50kb.tsv',header=T)
colnames(res_df) <- row.names(res_df)
distance <- dist(1-res)
res <- round(res, 2)
annotation_df <- data.frame(condition=factor(c(rep('GFP',3),rep('Ngn2',2),rep('PmutNgn2',2)),levels=c('GFP','Ngn2','PmutNgn2')))
row.names(annotation_df) <- colnames(res)
s_cols <- sample_colors[2:4]
names(s_cols) <- levels(annotation_df$condition)
pdf(paste0('plots/figures/',out_f,'.pdf'),width=width,height=height)
hm <- ComplexHeatmap::pheatmap(res,name = 'SCC',color=cols,breaks=seq(0.8,1,length=101),clustering_method = 'ward.D2',show_colnames = F,
                               annotation_colors=list(condition=s_cols),
                               display_numbers = TRUE, number_color = "black",fontsize_number = 10,border_color = "black",annotation_col=annotation_df, annotation_names_col = FALSE,legend_breaks=c(0.8,0.9,1))
draw(hm,merge_legend = TRUE)
dev.off()
dev.off()
