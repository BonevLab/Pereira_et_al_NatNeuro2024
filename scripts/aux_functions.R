intervals.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  return(inv)
}
intervals.2d.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  inv[,5:6] = floor((inv[,5]+inv[,6])/2)
  inv[,6] = inv[,5] + 1
  return(inv)
}

intervals.expand <- function(inv,expansion=100){
  inv[,2]<-inv[,2]-expansion
  inv[,3]<-inv[,3]+expansion
  return(gintervals.force_range(inv))
}

intervals.2d.expand <- function(inv,expansion1, expansion2){
  inv[,2]<-inv[,2]-expansion1
  inv[,3]<-inv[,3]+expansion1
  inv[,5]<-inv[,5]-expansion2
  inv[,6]<-inv[,6]+expansion2
  return(gintervals.force_range(inv))
}


intervals.size <- function(inv){
  return(sum(inv[,3]-inv[,2]))
}

intervals.normalize <- function(inv, size) {
  centers <- intervals.centers(inv)
  centers$end <- centers$end-1;
  return(intervals.expand(centers, floor(size/2)))
}

misha_extract <- function(tracks,regions=gintervals.all(),track_names=NULL,window=NULL,iterator,mode='avg',blacklist=NULL){
  for (track in tracks){
    #gvtrack.rm(paste0('v_',track))
    gvtrack.create(paste0('v_',track),track,mode)
    if(!is.null(window)){
      gvtrack.iterator(paste0('v_',track),sshift = -window/2,eshift = window/2)
    }
  }
  if(!grepl('global.percentile',mode)){
    df <- gextract(paste0('v_',tracks),regions,colnames = track_names,iterator= iterator)
  } else{
    df <- gextract(paste0('-log2(1-v_',tracks,')'),regions,colnames = track_names,iterator= iterator)
  }
  # df <- df[complete.cases(df),]
  if(!is.null(blacklist)){
    bl <- read.table(blacklist)
    bl <- gintervals(bl[,1],bl[,2],bl[,3])
    df_overlap <- gintervals.neighbors(df,bl,na.if.notfound = T)
    df <- df[!(df_overlap$dist<1e3),]
  }
  return(df)
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

centeredTF_plot<- function(pwm,chip_f,motif.name,genome,cutoff=5e-5,tracks,min, max, bin,lab,con,ylim,colours,save_plot,bed_dir,out_file,label,ylab,height=4,width=4){
  if(!file.exists(paste0(bed_dir,motif.name,'.bed'))){
    motif.names <- capitalize(tolower(as.vector(name(pwm))))
    if(!is.null(chip_f)){
      chip <- read.table(chip_f)[,1:3]
      colnames(chip) <- c('chrom','start','end')
      df <- makeGRangesFromDataFrame(chip)
    }
    motif_pos <- unlist(matchMotifs(pwm[[which(motif.names==motif.name)[1]]], df, genome = genome, out = "positions",p.cutoff = cutoff )) #added the which(motif.names==motif.name)[1] sinc ePax6 is double in the rds
    motif_6 <- data.frame(seqnames=seqnames(motif_pos),starts=start(motif_pos),ends=end(motif_pos),names=c(rep(".", length(motif_pos))),scores=elementMetadata(motif_pos[,1]),strands=strand(motif_pos))
    
    write.table(motif_6,paste0(bed_dir,motif.name,'.bed'),quote=F,sep='\t',col.names=F,row.names=F)
  }
  
  res_CG_CTCF <- getPlotSetArray(tracks=tracks,features=paste0(bed_dir,motif.name,'.bed'),
                                 refgenome=genome,type = 'mf',add_heatmap=F,xmin=min,xmax=max,bin = bin, ignore_strand = F)
  
  
  if(isTRUE(save_plot)){pdf(file=out_file,height=height,width=width)} 
  plotAverage(plotset=res_CG_CTCF, labels = lab, xlim = NULL,
              ylim =ylim, main = NULL, xlab = "", ylab = ylab,
              plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
              legend = T, legend_ext = F, legend_pos = 'topright',
              legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
              cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = colours)+
    axis(side = 1,labels=c(paste0('-',min/1000,'kb'),paste0(label),paste0('+',min/1000,'kb')),at=c(-min,0,max),pos = ylim[1],tick = T)
  if(isTRUE(save_plot)){dev.off()}
}

correlate_hic <- function(hic_files,bin,hic_names,h=10,dBPMax=2e6,temp_dir){
  res_df <- matrix(NA,nrow = length(hic_files),ncol=length(hic_files),dimnames = list(hic_names,hic_names))
  for (i in 1:length(hic_files)){
    for (j in 1:length(hic_files)){
      if (i==j){
        res_df[i,j] <- 1
      } else if (is.na(res_df[i,j])&is.na(res_df[j,i])){
        out_f <- tempfile(tmpdir = temp_dir)
        command_f <- paste0("/work/shared/mamba/miniforge3/envs/cooltools/bin/hicrep --binSize ",as.integer(bin)," --h ",h," --dBPMax ",as.integer(dBPMax)," ",hic_files[i]," ",hic_files[j]," ",out_f)
        system(command_f,ignore.stderr = T)
        df <- read.table(out_f,sep='\t',comment.char = "#")
        res_df[i,j] <- mean(df[,1],na.rm=T)
        res_df[j,i] <- mean(df[,1],na.rm=T)
        system(paste0("rm ",out_f))
      }
    }
  }
  return(res_df)
}

generate_rna_quantile <- function(tss_f,fpkm_f,quantiles_thresh=c(0,0.25,0.50,0.75,1),out_name='ES_TSS',out_f=paste0(main_dir,'data/mm10/beds/')){
  tss <-read.table(tss_f)
  fpkm <-read.table(fpkm_f,header=T)
  colnames(fpkm)[1:3] <- paste0('rep',1:3)
  inactive_tss <- tss[tss$V4%in%fpkm$gene_name[fpkm$means<1],]
  write.table(inactive_tss,paste0(out_f,out_name,'_inactive','.bed'),quote=F,row.names=F,col.names=F,sep='\t')
  fpkm <- fpkm[fpkm$means>=1,]
  quantiles <- quantile(fpkm$means,quantiles_thresh)
  for (i in 1:(length(quantiles)-1)){
    fpkm_quants <- fpkm[fpkm$means>quantiles[i]&fpkm$means<=quantiles[i+1],]
    temp_tss <- tss[tss$V4%in%fpkm_quants$gene_name,]
    write.table(temp_tss,paste0(out_f,out_name,'_Q',i,'.bed'),quote=F,row.names=F,col.names=F,sep='\t')
  }
}


filter_smf <- function(df,binSize,filter_dist,ncalls=1,filter_interv_2d=NULL, shuffleSet = "none", dist=0,filter_interv_l=NULL,ldist=0,filter_interv_r=NULL,rdist=0,out_res_coord=FALSE,shufflebS_lim=NULL){
  lcoord_mat <- df[,grepl('call_coord',colnames(df))&grepl('\\.x',colnames(df))]
  rcoord_mat <- df[,grepl('call_coord',colnames(df))&grepl('\\.y',colnames(df))]
  linterv_mat <- df[,grepl('interv_coord',colnames(df))&grepl('\\.x',colnames(df))]
  rinterv_mat <- df[,grepl('interv_coord',colnames(df))&grepl('\\.y',colnames(df))]
  lcall_mat <- df[,grepl('^call',colnames(df))&!grepl('coord',colnames(df))&grepl('\\.x',colnames(df))]
  rcall_mat <- df[,grepl('^call',colnames(df))&!grepl('coord',colnames(df))&grepl('\\.y',colnames(df))]
  if((ncol(lcoord_mat)!=ncol(linterv_mat))|(ncol(rcoord_mat)!=ncol(rinterv_mat))){
    lncol <- min(ncol(lcoord_mat),ncol(linterv_mat),ncol(lcall_mat))
    rncol <- min(ncol(rcoord_mat),ncol(rinterv_mat),ncol(rcall_mat))
    lcoord_mat <- lcoord_mat[,1:lncol]
    linterv_mat <- linterv_mat[,1:lncol]
    lcall_mat <- lcall_mat[,1:lncol]
    rcoord_mat <- rcoord_mat[,1:rncol]
    rinterv_mat <- rinterv_mat[,1:rncol]
    rcall_mat <- rcall_mat[,1:rncol]
  }
  lcoord_mat <- (lcoord_mat-linterv_mat)[,1:ncol(lcall_mat)]
  rcoord_mat <- (rcoord_mat-rinterv_mat)[,1:ncol(rcall_mat)]
  linterv_mat <- linterv_mat[,1:ncol(lcall_mat)]
  rinterv_mat <- rinterv_mat[,1:ncol(rcall_mat)]
  if(shuffleSet=='R1') {
      ldist_idx <- as.matrix(abs(lcoord_mat)>500 & abs(lcoord_mat)<shufflebS_lim)
      ldist_idx[is.na(ldist_idx)] <- FALSE
      #new here
      rdist_idx <- as.matrix(abs(rcoord_mat)<=binSize)
      rdist_idx[is.na(rdist_idx)] <- FALSE
    } else if (shuffleSet == 'R2') {
      ldist_idx <- as.matrix(abs(lcoord_mat)<=binSize)
      ldist_idx[is.na(ldist_idx)] <- FALSE
      #new here
      rdist_idx <- as.matrix(abs(rcoord_mat)>500 & abs(rcoord_mat)<shufflebS_lim)
      rdist_idx[is.na(rdist_idx)] <- FALSE
    } else {
    ldist_idx <- as.matrix(abs(lcoord_mat)<=binSize)
    ldist_idx[is.na(ldist_idx)] <- FALSE
    rdist_idx <- as.matrix(abs(rcoord_mat)<=binSize)
    rdist_idx[is.na(rdist_idx)] <- FALSE
  }
  lcall_mat[!ldist_idx] <- NA
  rcall_mat[!rdist_idx] <- NA 
  lcoord_mat[!ldist_idx] <- NA
  rcoord_mat[!rdist_idx] <- NA
  linterv_mat[!ldist_idx] <- NA
  rinterv_mat[!rdist_idx] <- NA
  
  binSize_idx <- intersect(as.numeric(which(rowSums(ldist_idx)>=1)),as.numeric(which(rowSums(rdist_idx)>=1)))
  
  filter_dist_idx <- which(abs(rowMeans(as.matrix(rinterv_mat),na.rm=T)-rowMeans(as.matrix(linterv_mat),na.rm=T))>filter_dist[1]&abs(rowMeans(as.matrix(rinterv_mat),na.rm=T)-rowMeans(as.matrix(linterv_mat),na.rm=T))<=filter_dist[2])
  #filter_dist_idx <- which(rowMins(as.matrix(abs(rinterv_mat-linterv_mat)),na.rm=T)>filter_dist[1]&rowMaxs(as.matrix(abs(rinterv_mat-linterv_mat)),na.rm=T)<=filter_dist[2])
  #filter_dist_idx <- which((df$interv_distance>=filter_dist[1])&(df$interv_distance<filter_dist[2]))
  ncalls_idx <- which((as.numeric(rowSums(ldist_idx))>=ncalls)&(as.numeric(rowSums(rdist_idx))>=ncalls))
  idx <- intersect(binSize_idx,filter_dist_idx)
  idx <- intersect(idx,ncalls_idx)
#  #NEW#- sampling total intersect as expand.grid function below takes too long to run with full dataset ## REVIEW
#  if(shuffleSet != "None"){
#    sample(idx, 4000)
#  }else{
#    idx}
  res_coords <- data.frame(chrom1=NA,start1=NA,end1=NA,chrom2=NA,start2=NA,end2=NA,index=NA)[-1,]
  if(out_res_coord){
    for (i in idx){
      temp_df <- expand.grid(as.numeric(as.vector(linterv_mat[i,])),as.numeric(as.vector(rinterv_mat[i,])))
      temp_df <- unique(temp_df[complete.cases(temp_df),])
      temp_df <- data.frame(chrom1=df$chr.x[i],start1=temp_df$Var1,end1=temp_df$Var1+1,chrom2=df$chr.x[i],start2=temp_df$Var2,end2=temp_df$Var2+1,index=i)
      res_coords <- rbind(res_coords,temp_df)
    }
    res_coords <- unique(res_coords)
  }
  if (!is.null(filter_interv_2d)){
    colnames(filter_interv_2d) <- c('chrom1','start1','end1','chrom2','start2','end2')
    df_intersect <- gintervals.neighbors(res_coords,filter_interv_2d)
    filter_interv_idx <- which((abs(df_intersect$dist1)<=dist)&(abs(df_intersect$dist2)<=dist))
    filter_interv_idx <- res_coords$index[filter_interv_idx]
    idx <- intersect(idx,filter_interv_idx)
    res_coords <- res_coords[res_coords$index%in%idx,]
  }
  if (!is.null(filter_interv_l)){
    temp <- res_coords[,1:3]
    colnames(temp) <- c('chrom','start','end')
    colnames(filter_interv_l) <- c('chrom','start','end')
    df_intersect <- gintervals.neighbors(temp,filter_interv_l)
    filter_interv_idx <- which(abs(df_intersect$dist)<=ldist)
    filter_interv_idx <- res_coords$index[filter_interv_idx]
    idx <- intersect(idx,filter_interv_idx)
    res_coords <- res_coords[res_coords$index%in%idx,]
  }
  if (!is.null(filter_interv_r)){
    temp <- res_coords[,4:6]
    colnames(temp) <- c('chrom','start','end')
    colnames(filter_interv_r) <- c('chrom','start','end')
    df_intersect <- gintervals.neighbors(temp,filter_interv_r)
    filter_interv_idx <- which(abs(df_intersect$dist)<=rdist)
    filter_interv_idx <- res_coords$index[filter_interv_idx]
    idx <- intersect(idx,filter_interv_idx)
    res_coords <- res_coords[res_coords$index%in%idx,]
  }
  res_ratios <- cbind(rowMeans(lcall_mat[idx,],na.rm=T),rowMeans(rcall_mat[idx,],na.rm=T))
  res <- list(res_ratios=res_ratios,res_coords=res_coords,
              lcoord_mat=lcoord_mat[idx,],rcoord_mat=rcoord_mat[idx,],
              lcall_mat=lcall_mat[idx,],rcall_mat=rcall_mat[idx,],
              linterv_mat=linterv_mat[idx,],rinterv_mat=rinterv_mat[idx,],
              idx=idx)
  return(res)
}

smfCoveragePlot <- function(
  res,
  pair_smf_f,
  binSize,
  mat,
  window = 1,
  cluster_cols=hue_pal()(length(unique(mat$cluster))),
  meth_colors=c('grey','black'),
  max_binSize=500,
  window2=5,
  point.size=4,
  plot_binSize=2*binSize
) {
  lcoord_mat <- pair_smf_f[,grepl('call_coord',colnames(pair_smf_f))&grepl('\\.x',colnames(pair_smf_f))]
  rcoord_mat <- pair_smf_f[,grepl('call_coord',colnames(pair_smf_f))&grepl('\\.y',colnames(pair_smf_f))]
  linterv_mat <- pair_smf_f[,grepl('interv_coord',colnames(pair_smf_f))&grepl('\\.x',colnames(pair_smf_f))]
  rinterv_mat <- pair_smf_f[,grepl('interv_coord',colnames(pair_smf_f))&grepl('\\.y',colnames(pair_smf_f))]
  lcall_mat <- pair_smf_f[,grepl('^call',colnames(pair_smf_f))&!grepl('coord',colnames(pair_smf_f))&grepl('\\.x',colnames(pair_smf_f))]
  rcall_mat <- pair_smf_f[,grepl('^call',colnames(pair_smf_f))&!grepl('coord',colnames(pair_smf_f))&grepl('\\.y',colnames(pair_smf_f))]
  if((ncol(lcoord_mat)!=ncol(linterv_mat))|(ncol(rcoord_mat)!=ncol(rinterv_mat))){
    lncol <- min(ncol(lcoord_mat),ncol(linterv_mat),ncol(lcall_mat))
    rncol <- min(ncol(rcoord_mat),ncol(rinterv_mat),ncol(rcall_mat))
    lcoord_mat <- lcoord_mat[,1:lncol]
    linterv_mat <- linterv_mat[,1:lncol]
    lcall_mat <- lcall_mat[,1:lncol]
    rcoord_mat <- rcoord_mat[,1:rncol]
    rinterv_mat <- rinterv_mat[,1:rncol]
    rcall_mat <- rcall_mat[,1:rncol]
  }
  lcoord_mat <- lcoord_mat-linterv_mat
  rcoord_mat <- rcoord_mat-rinterv_mat
  lcall_mat <- lcall_mat[res$idx,]
  rcall_mat <- rcall_mat[res$idx,]
  lcoord_mat <- lcoord_mat[res$idx,1:ncol(lcall_mat)]
  rcoord_mat <- rcoord_mat[res$idx,1:ncol(rcall_mat)]
  
  reads1 <- data.frame(chrom='chr1',
                       start=as.numeric(c(t(as.matrix(lcoord_mat))))+max_binSize+1,
                       end=as.numeric(c(t(as.matrix(lcoord_mat))))+max_binSize+1,
                       read=as.character(rep(pair_smf_f$read_name[res$idx],each=ncol(lcoord_mat))),
                       count=as.numeric(c(t(as.matrix(lcall_mat)))),
                       group=rep(mat$cluster,each=ncol(lcoord_mat)))
  reads1 <- reads1[complete.cases(reads1),]
  
  reads2 <- data.frame(chrom='chr1',
                       start=as.numeric(c(t(as.matrix(rcoord_mat))))+max_binSize+1,
                       end=as.numeric(c(t(as.matrix(rcoord_mat))))+max_binSize+1,
                       read=as.character(rep(pair_smf_f$read_name[res$idx],each=ncol(rcoord_mat))),
                       count=as.numeric(c(t(as.matrix(rcall_mat)))),
                       group=rep(mat$cluster,each=ncol(rcoord_mat)))
  reads2 <- reads2[complete.cases(reads2),]
  cells.per.group <- table(mat$cluster)
  reads1 <- makeGRangesFromDataFrame(reads1,keep.extra.columns = T)
  reads2 <- makeGRangesFromDataFrame(reads2,keep.extra.columns = T)
  
  region <- GRanges(seqnames='chr1', ranges=IRanges(1,2*max_binSize+1), strand=NULL)
  sliding_region <- slidingWindows(region,window,window)[[1]]
  mat1 <- binReads(sliding_region,reads1,window,max_binSize)
  mat2 <- binReads(sliding_region,reads2,window,max_binSize)
  df.long1 <- melt(as.matrix(mat1))
  df.long1$group <- reads1$group[match(df.long1$Var1,reads1$read)]
  df.long2 <- melt(as.matrix(mat2))
  df.long2$group <- reads2$group[match(df.long2$Var1,reads2$read)]
  cells_group <- ''
  plotLabel <- as.list(cells_group)
  plot_labeller <- function(variable,value){
    return(plotLabel[value])
  }
  df_sub1 <- subset(df.long1, value != 0)
  df_sub1$value[df_sub1$value==(-1)] <- 0
  df_sub1$color <- df_sub1$group
  df_sub1 <- df_sub1[sample(row.names(df_sub1),nrow(df_sub1)),] %>% arrange(group) %>% mutate(Var1=factor(Var1, levels=unique(Var1)))
  row.names(df_sub1) <- 1:nrow(df_sub1)
  
  df_sub2 <- subset(df.long2, value != 0)
  df_sub2$value[df_sub2$value==(-1)] <- 0
  df_sub2$color <- df_sub2$group
  df_sub2$Var1 <- factor(df_sub2$Var1,levels=levels(df_sub1$Var1))
  df_sub2 <- df_sub2[order(as.numeric(df_sub2$Var1)),]
  row.names(df_sub2) <- 1:nrow(df_sub2)
  
  p1 <- ggplot(data = df_sub1, mapping = aes(x = Var2, y = Var1, fill = as.factor(value),size=point.size),order=Var1) +
    geom_tile() + scale_x_continuous(expand = c(0,0),limits=c(-plot_binSize,plot_binSize),breaks = breaks_extended(3)) +
    scale_fill_manual(name='',values=meth_colors,na.value = 'white')+
    facet_grid(rows=vars(group), scales="free_y",space='free_y',labeller =plot_labeller,margins=F, switch="y") +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y=element_blank(),
      axis.text.x=element_text(size=rel(0.5)),
      panel.spacing = unit(0.010, "lines"),
      legend.position = 'none',
      strip.text.y = element_text(hjust=-10),
      strip.background = element_rect(
        color="black", fill=cluster_cols, size=1, linetype="solid"
      ),
      plot.margin = unit(c(0,0,0,0), "cm")
    )
  g <- ggplot_gtable(ggplot_build(p1))
  stripr <- which(grepl('strip-l', g$layout$name))
  fills <- cluster_cols
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  p1 <- g
  plotlist <- list(p1)
  rel.heights=c(1)    
  p_la <- suppressWarnings(cowplot::plot_grid(
    plotlist = plotlist,
    ncol = 1,
    axis = 'btlr',
    rel_heights = rel.heights,
    align = 'v',
    greedy = FALSE
  ))
  
  p1 <- ggplot(data = df_sub2, mapping = aes(x = Var2, y = Var1, fill = as.factor(value),size=point.size),order=Var1) +
    geom_tile() + scale_x_continuous(expand = c(0,0),limits=c(-plot_binSize,plot_binSize), breaks = breaks_extended(3)) +
    scale_fill_manual(name='',values=meth_colors,na.value = 'white')+
    facet_grid(rows=vars(group), scales="free_y",space='free_y',labeller =plot_labeller,margins=F, switch="both") +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y=element_blank(),
      axis.text.x=element_text(size=rel(0.5)),
      panel.spacing = unit(0.010, "lines"),
      legend.position = 'none',
      strip.text.y = element_text(hjust=-10),
      strip.background = element_rect(
        color="white", fill="white", size=1, linetype="solid"
      ),
      plot.margin = unit(c(0,0,0,0), "cm")
    )
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
  plotlist <- list(p1)
  rel.heights=c(1)    
  p_ra <- suppressWarnings(cowplot::plot_grid(
    plotlist = plotlist,
    ncol = 1,
    axis = 'btlr',
    rel_heights = rel.heights,
    align = 'v',
    greedy = FALSE
  ))
  # p <- p + theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  message('Succesfully extracted methylation calls')
  mat1 <- as.matrix(mat1)
  mat1[mat1==0] <- NA
  mat1[mat1==(-1)] <- 0
  mat2 <- as.matrix(mat2)
  mat2[mat2==0] <- NA
  mat2[mat2==(-1)] <- 0
  sliding_region <- slidingWindows(region,window2,window2)[[1]]
  lmat2=as.matrix(binReads(sliding_region,reads1,window2,max_binSize))
  lmat2[lmat2==0] <- NA
  lmat2[lmat2==(-1)] <- 0
  rmat2=as.matrix(binReads(sliding_region,reads2,window2,max_binSize))
  rmat2[rmat2==0] <- NA
  rmat2[rmat2==(-1)] <- 0
  return(list(p_la=p_la,p_ra=p_ra,lmat=mat1,rmat=mat2,lmat2=lmat2,rmat2=rmat2))
}


binReads <- function(query, fragments,window,binSize=max_binSize){
  #Count By Fragments Insertions
  inserts <- GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), end(fragments)), RG = mcols(fragments))
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], "RG.read"]
  overlapDF$count <- mcols(inserts)[overlapDF[, 2], "RG.count"]
  overlapDF$count[overlapDF$count==0] <- -1
  overlapDF$count[overlapDF$count>1] <- 1
  overlapDF$count[overlapDF$count<(-1)] <- -1
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  
  overlapTDF <- as.data.frame(as.data.frame(overlapTDF) %>% group_by(queryHits,name) %>% summarise(count=round(mean(count),0),id=id) %>% arrange(queryHits,id))
  
  #Summarize
  sparseM <- t(Matrix::sparseMatrix(
    i = overlapTDF[, 'queryHits'], 
    j = overlapTDF[, 'id'],
    x = overlapTDF[, 'count'], 
    dims = c(length(query), length(unique(fragments$read)))))
  row.names(sparseM) <- unique(overlapTDF$name)
  colnames(sparseM) <- seq(-binSize,binSize,by = window)
  sparseM@x[sparseM@x>1] <- 1
  sparseM@x[sparseM@x<(-1)] <- -1
  return(sparseM)
}

extract_hic_score <- function(res,plot_res,pair_smf_f,tracks,expand=c(-5e3,5e3),add_shuffle=FALSE,shift_r=0){
  options(gmax.data.size=5e7)
  options(gmultitasking=F)
  for (i in 1:length(tracks)){
   # gvtrack.rm(paste0('v_',tracks[i]))
    gvtrack.create(paste0('v_',tracks[i]),tracks[i],'avg')
    gvtrack.iterator.2d(paste0('v_',tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
  }
  res_coords <- unique(res$res_coords)
  row.names(res_coords) <- 1:nrow(res_coords)
  for (i in 1:nrow(res_coords)){
    if(res_coords[i,2]>res_coords[i,5]){
      i_temp <- res_coords[i,1:3]
      res_coords[i,1:3] <- res_coords[i,4:6]
      res_coords[i,4:6] <- i_temp
    }
  }
  res_coords <- res_coords[!duplicated(res_coords[,-ncol(res_coords)]),]
  scores <- gextract(paste0('v_',tracks),intervals = res_coords,iterator = res_coords,band=-c(5e6,max(expand)))
  df_id <- data.frame(cluster=plot_res$cluster_id,index=row.names(pair_smf_f[res$idx,]))
  scores$index <- res_coords$index[match(scores$intervalID,row.names(res_coords))]
  scores <- merge(scores,df_id,by='index')
  if(add_shuffle){
    c_scores <- shuffle_regions(regions=scores[,-1],range=shift_r)
    c_scores <- gextract(paste0('v_',tracks),intervals = c_scores,iterator = c_scores,band=-c(5e6,max(expand)))
    c_scores$cluster <- 'control'
    c_scores$index <- NA
    scores <- rbind(scores,c_scores)
  }
  return(scores)
}

extract_linear_score <- function(res,plot_res,pair_smf_f,tracks,expand=c(-100,100),add_shuffle=FALSE,shift_r){
  options(gmax.data.size=5e7)
  options(gmultitasking=F)
  for (i in 1:length(tracks)){
    gvtrack.create(paste0('v_',tracks[i]),tracks[i],'avg')
    gvtrack.iterator(paste0('v_',tracks[i]), sshift=min(expand), eshift=max(expand))
  }
  res_coords <- unique(res$res_coords)
  res_coords <- res_coords[!duplicated(res_coords[,-ncol(res_coords)]),]
  row.names(res_coords) <- 1:nrow(res_coords)
  res_coords_l <- res_coords[,1:3]
  colnames(res_coords_l) <- c('chrom','start','end')
  scores_l <- gextract(paste0('v_',tracks),intervals = res_coords_l,iterator = res_coords_l)
  res_coords_r <- res_coords[,4:6]
  colnames(res_coords_r) <- c('chrom','start','end')
  scores_r <- gextract(paste0('v_',tracks),intervals = res_coords_r,iterator = res_coords_r)
  scores <- merge(scores_l,scores_r,by='intervalID')
  
  df_id <- data.frame(cluster=plot_res$cluster_id,index=row.names(pair_smf_f[res$idx,]))
  scores$index <- res_coords$index[match(scores$intervalID,row.names(res_coords))]
  scores <- merge(scores,df_id,by='index')
  if(add_shuffle){
    l_scores <- shuffle_regions(regions=scores[,3:5],range=shift_r)
    l_scores <- gextract(paste0('v_',tracks),intervals = l_scores,iterator = l_scores)
    r_scores <- shuffle_regions(regions=scores[,(ncol(scores_l)+2):(ncol(scores_l)+4)],range=shift_r)
    r_scores <- gextract(paste0('v_',tracks),intervals = r_scores,iterator = r_scores)
    c_scores <- merge(l_scores,r_scores,by='intervalID')
    c_scores$index <- NA
    c_scores$cluster <- 'control'
    scores <- rbind(scores,c_scores)
  }
  return(scores)
}

shuffle_regions <- function(regions,range){
  if(colnames(regions)[1]=='chrom1'){
   grid <- regions[,1:6] 
   r_shift <- sample(range,nrow(grid),replace=T)
   grid[,c(2:3,5:6)] <- grid[,c(2:3,5:6)] + r_shift
   grid <- gintervals.force_range(unique(grid))
   grid <- gintervals.canonic(grid)
   
  } else {
    grid <- regions[,1:3] 
    colnames(grid) <- c('chrom','start','end')
    r_shift <- sample(range,nrow(grid),replace=T)
    grid[,c(2:3)] <- grid[,c(2:3)] + r_shift
    grid <- gintervals.force_range(unique(grid))
    grid <- gintervals.canonic(grid)
  }
  return(grid)
}


getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)
  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
    seqnames=paste0("chr",seqnames(geneGTF1)),
    ranges=ranges(geneGTF1),
    strand=strand(geneGTF1),
    gene_name=geneGTF1$gene_name,
    gene_id=geneGTF1$gene_id
  ) %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}

extract_misha_ep_pairs <- function(enhancer_bed,tss_bed,tss_full,interval_window,min_dist,max_dist,expand,domains_f,score_tracks,hic_names,other_tracks,other_names){
  tss <- read.table(tss_bed)
  tss <- gintervals(tss[,1],as.numeric(tss[,2]),as.numeric(tss[,3]))
  tss <- gintervals.neighbors(tss_full,tss)
  tss <- tss[tss$dist==0,1:6]
#  colnames(tss) <- c('chrom','start','end','geneName','score','strand')
  enh <- read.table(enhancer_bed)
  enh <- gintervals(enh[,1],enh[,2],enh[,3])
  enh <- intervals.centers(enh)
  domains <- gintervals.load(domains_f)
  dom_2d = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
  
  for (i in 1:length(score_tracks)){
    gvtrack.create(paste0('v_',score_tracks[i]),score_tracks[i],'max')
    gvtrack.iterator.2d(paste0('v_',score_tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
  }
  
  for (i in 1:length(other_tracks)){
    gvtrack.create(paste0('v_',other_tracks[i]),other_tracks[i],'avg')
    gvtrack.iterator(paste0('v_',other_tracks[i]), sshift=-interval_window/2, eshift=interval_window/2)
  }
  
  ########################
  
  grid1 <- unique(construct.grid(enh,tss,min_dist,max_dist))
  grid1 <- gintervals.canonic(grid1)
  grid2 <- unique(construct.grid(tss,enh,min_dist,max_dist))
  grid2 <- gintervals.canonic(grid2)
  intra_grid1 <- gintervals.intersect(grid1,dom_2d)
  intra_grid2 <- gintervals.intersect(grid2,dom_2d)
  inter_grid1 <- grid1[!(paste0(grid1[,1],':',grid1[,2],'-',grid1[,4],':',grid1[,5])%in%paste0(intra_grid1[,1],':',intra_grid1[,2],'-',intra_grid1[,4],':',intra_grid1[,5])),]
  inter_grid2 <- grid2[!(paste0(grid2[,1],':',grid2[,2],'-',grid2[,4],':',grid2[,5])%in%paste0(intra_grid2[,1],':',intra_grid2[,2],'-',intra_grid2[,4],':',intra_grid2[,5])),]
  intra_scores1<- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = intra_grid1,iterator = intra_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
  intra_scores2 <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = intra_grid2,iterator = intra_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
  inter_scores1 <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = inter_grid1,iterator = inter_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
  inter_scores2 <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = inter_grid2,iterator = inter_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
  tss_intra_scores1 <- intra_scores1[,4:6]
  colnames(tss_intra_scores1) <- c('chrom','start','end')
  intra_scores1$geneName <- gintervals.neighbors(tss_intra_scores1,tss)$geneName
  tss_intra_scores2 <- intra_scores2[,1:3]
  colnames(tss_intra_scores2) <- c('chrom','start','end')
  intra_scores2$geneName <- gintervals.neighbors(tss_intra_scores2,tss)$geneName
  tss_inter_scores1 <- inter_scores1[,4:6]
  colnames(tss_inter_scores1) <- c('chrom','start','end')
  inter_scores1$geneName <- gintervals.neighbors(tss_inter_scores1,tss)$geneName
  tss_inter_scores2 <- inter_scores2[,1:3]
  colnames(tss_inter_scores2) <- c('chrom','start','end')
  inter_scores2$geneName <- gintervals.neighbors(tss_inter_scores2,tss)$geneName
  intra_scores2 <- intra_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
  colnames(intra_scores2)[1:6] <- colnames(intra_scores1)[1:6]
  inter_scores2 <- inter_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
  colnames(inter_scores2)[1:6] <- colnames(inter_scores1)[1:6]
  intra_scores <- rbind(intra_scores1,intra_scores2)
  inter_scores <- rbind(inter_scores1,inter_scores2)
  intra_scores$domain <- 'intraTAD'
  inter_scores$domain <- 'interTAD'
  scores <- rbind(intra_scores,inter_scores)
  distal_scores <- gextract(paste0('v_',other_tracks),colnames=paste0(other_names,'_distal'),intervals = intervals.centers(enh),iterator = intervals.centers(enh))
  distal_scores$intervalID <- paste0(distal_scores$chrom,':',distal_scores$start)
  tss_scores <- gextract(paste0('v_',other_tracks),colnames=paste0(other_names,'_tss'),intervals = intervals.centers(tss),iterator = intervals.centers(tss))
  tss_scores$intervalID <- paste0(tss_scores$chrom,':',tss_scores$start)
  
  scores$distalIDs <- paste0(scores$chrom1,':',scores$start1)
  scores$tssIDs <- paste0(scores$chrom2,':',scores$start2)
  res <- cbind(scores,distal_scores[match(scores$distalIDs,distal_scores$intervalID),-c(1:3)],tss_scores[match(scores$tssIDs,tss_scores$intervalID),-c(1:3)])
  res <- res[,grep('tssID',colnames(res),invert=T)]
  res <- res[,grep('intervalID',colnames(res),invert=T)]
  row.names(res) <- 1:nrow(res)
  return(res)
}

subset_misha_pairs <- function(res,distal_bed,dist=0){
  res_idx <- res[,1:3]
  colnames(res_idx) <- c('chrom','start','end')
  dist_peaks <- read.table(distal_bed)
  dist_peaks <- gintervals(dist_peaks[,1],dist_peaks[,2],dist_peaks[,3])
  res_idx <- gintervals.neighbors(res_idx,dist_peaks,na.if.notfound = T)
  return(res[res_idx$dist<=dist,])
}

extract_misha_pairs <- function(bed1,bed2,interval_window,min_dist,max_dist,expand,domains_f,score_tracks,hic_names,other_tracks,other_names,intra_only=F){
  interv1 <- read.table(bed1)
  interv1 <- interv1[interv1[,1]%in%gintervals.all()$chrom,]
  interv1 <- gintervals(interv1[,1],interv1[,2],interv1[,3])
  interv1 <- intervals.centers(interv1)
  interv2 <- read.table(bed2)
  interv2 <- interv2[interv2[,1]%in%gintervals.all()$chrom,]
  interv2 <- gintervals(interv2[,1],interv2[,2],interv2[,3])
  interv2 <- intervals.centers(interv2)
  
  domains <- gintervals.load(domains_f)
  dom_2d = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
  
  for (i in 1:length(score_tracks)){
    gvtrack.create(paste0('v_',score_tracks[i]),score_tracks[i],'avg')
    gvtrack.iterator.2d(paste0('v_',score_tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
  }
  
  for (i in 1:length(other_tracks)){
    gvtrack.create(paste0('v_',other_tracks[i]),other_tracks[i],'avg')
    gvtrack.iterator(paste0('v_',other_tracks[i]), sshift=-interval_window/2, eshift=interval_window/2)
  }
  
  ########################
  
  grid1 <- unique(construct.grid(interv1,interv2,min_dist,max_dist))
  grid1 <- gintervals.canonic(grid1)
  ###Test###
 # if(nrow(grid1)>1e4){
 #   grid1 <- grid1[sample(row.names(grid1),1e4,replace = F),]
 # }
  ####
  intra_grid1 <- gintervals.intersect(grid1,dom_2d)
  inter_grid1 <- grid1[!(paste0(grid1[,1],':',grid1[,2],'-',grid1[,4],':',grid1[,5])%in%paste0(intra_grid1[,1],':',intra_grid1[,2],'-',intra_grid1[,4],':',intra_grid1[,5])),]
  intra_scores<- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = intra_grid1,iterator = intra_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
  
  intra_scores$domain <- 'intraTAD'
  
  if(intra_only){
    scores <- intra_scores
  } else {
    inter_scores <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = inter_grid1,iterator = inter_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
    inter_scores$domain <- 'interTAD'
    scores <- rbind(intra_scores,inter_scores)
  }
  
  r1_scores <- gextract(paste0('v_',other_tracks),colnames=paste0(other_names,'_r1'),intervals = interv1,iterator = interv1)
  r1_scores$intervalID <- paste0(r1_scores$chrom,':',r1_scores$start)
  r2_scores <- gextract(paste0('v_',other_tracks),colnames=paste0(other_names,'_r2'),intervals = interv2,iterator = interv2)
  r2_scores$intervalID <- paste0(r2_scores$chrom,':',r2_scores$start)
  
  scores$r1IDs <- paste0(scores$chrom1,':',scores$start1)
  scores$r2IDs <- paste0(scores$chrom2,':',scores$start2)
  res <- cbind(scores,r1_scores[match(scores$r1IDs,r1_scores$intervalID),-c(1:3)],r2_scores[match(scores$r2IDs,r2_scores$intervalID),-c(1:3)])
  res <- res[,grep('intervalID',colnames(res),invert=T)]
  row.names(res) <- 1:nrow(res)
  res <- res[,-c(1:6)]
  res <- res[complete.cases(res[,grep('ID',colnames(res),invert=T)]),]
  if(length(other_names)==2){
    cor_res <- c(cor.test(res[,hic_names],rowMeans(res[,paste0(other_names[1],c('_r1','_r2'))]))$estimate,
                 cor.test(res[,hic_names],rowMeans(res[,paste0(other_names[2],c('_r1','_r2'))]))$estimate,
                 cor.test(res[!duplicated(res$r1IDs),paste0(other_names[1],'_r1')],res[!duplicated(res$r1IDs),paste0(other_names[2],'_r1')])$estimate)
    names(cor_res) <- c(paste0(hic_names,'vs',other_names[1]),paste0(hic_names,'vs',other_names[2]),paste0(other_names[1],'vs',other_names[2]))
    return(cor_res)
  } else {
    return(res)
  }
  
}
# 
# extract_misha_ep_pairs <- function(enhancer_bed,tss_bed,interval_window,min_dist,max_dist,expand,domains_f,score_tracks,hic_names,other_tracks,other_names){
#   tss <- read.table(tss_bed)
#   colnames(tss) <- c('chrom','start','end','geneName','score','strand')
#   enh <- read.table(enhancer_bed)
#   enh <- gintervals(enh[,1],enh[,2],enh[,3])
#   
#   domains <- gintervals.load(domains_f)
#   dom_2d = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
#   
#   for (i in 1:length(score_tracks)){
#     gvtrack.create(paste0('v_',score_tracks[i]),score_tracks[i],'max')
#     gvtrack.iterator.2d(paste0('v_',score_tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
#   }
#   
#   for (i in 1:length(other_tracks)){
#     gvtrack.create(paste0('v_',other_tracks[i]),other_tracks[i],'avg')
#     gvtrack.iterator(paste0('v_',other_tracks[i]), sshift=-interval_window/2, eshift=interval_window/2)
#   }
#   
#   ########################
#   
#   grid1 <- unique(construct.grid(enh,tss,min_dist,max_dist))
#   grid1 <- gintervals.canonic(grid1)
#   grid2 <- unique(construct.grid(tss,enh,min_dist,max_dist))
#   grid2 <- gintervals.canonic(grid2)
#   
#   intra_grid1 <- gintervals.intersect(grid1,dom_2d)
#   intra_grid2 <- gintervals.intersect(grid2,dom_2d)
#   inter_grid1 <- grid1[!(paste0(grid1[,1],':',grid1[,2],'-',grid1[,4],':',grid1[,5])%in%paste0(intra_grid1[,1],':',intra_grid1[,2],'-',intra_grid1[,4],':',intra_grid1[,5])),]
#   inter_grid2 <- grid2[!(paste0(grid2[,1],':',grid2[,2],'-',grid2[,4],':',grid2[,5])%in%paste0(intra_grid2[,1],':',intra_grid2[,2],'-',intra_grid2[,4],':',intra_grid2[,5])),]
#   
#   intra_scores1<- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = intra_grid1,iterator = intra_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
#   intra_scores2 <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = intra_grid2,iterator = intra_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
#   inter_scores1 <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = inter_grid1,iterator = inter_grid1,band = -c(max_dist+max(expand),min_dist-max(expand)))
#   inter_scores2 <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = inter_grid2,iterator = inter_grid2,band = -c(max_dist+max(expand),min_dist-max(expand)))
#   
#   tss_intra_scores1 <- intra_scores1[,4:6]
#   colnames(tss_intra_scores1) <- c('chrom','start','end')
#   intra_scores1$geneName <- gintervals.neighbors(tss_intra_scores1,tss)$geneName
#   tss_intra_scores2 <- intra_scores2[,1:3]
#   colnames(tss_intra_scores2) <- c('chrom','start','end')
#   intra_scores2$geneName <- gintervals.neighbors(tss_intra_scores2,tss)$geneName
#   tss_inter_scores1 <- inter_scores1[,4:6]
#   colnames(tss_inter_scores1) <- c('chrom','start','end')
#   inter_scores1$geneName <- gintervals.neighbors(tss_inter_scores1,tss)$geneName
#   tss_inter_scores2 <- inter_scores2[,1:3]
#   colnames(tss_inter_scores2) <- c('chrom','start','end')
#   inter_scores2$geneName <- gintervals.neighbors(tss_inter_scores2,tss)$geneName
#   
#   
#   intra_scores2 <- intra_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
#   colnames(intra_scores2)[1:6] <- colnames(intra_scores1)[1:6]
#   inter_scores2 <- inter_scores2[,c(4:6,1:3,7:ncol(intra_scores2))]
#   colnames(inter_scores2)[1:6] <- colnames(inter_scores1)[1:6]
#   intra_scores <- rbind(intra_scores1,intra_scores2)
#   inter_scores <- rbind(inter_scores1,inter_scores2)
#   intra_scores$domain <- 'intraTAD'
#   inter_scores$domain <- 'interTAD'
#   scores <- rbind(intra_scores,inter_scores)
#   
#   
#   distal_scores <- gextract(paste0('v_',other_tracks),colnames=paste0(other_names,'_distal'),intervals = intervals.centers(enh),iterator = intervals.centers(enh))
#   distal_scores$intervalID <- paste0(distal_scores$chrom,':',distal_scores$start)
#   tss_scores <- gextract(paste0('v_',other_tracks),colnames=paste0(other_names,'_tss'),intervals = intervals.centers(tss),iterator = intervals.centers(tss))
#   tss_scores$intervalID <- paste0(tss_scores$chrom,':',tss_scores$start)
#   
#   scores$distalIDs <- paste0(scores$chrom1,':',scores$start1)
#   scores$tssIDs <- paste0(scores$chrom2,':',scores$start2)
#   res <- cbind(scores,distal_scores[match(scores$distalIDs,distal_scores$intervalID),-c(1:3)],tss_scores[match(scores$tssIDs,tss_scores$intervalID),-c(1:3)])
#   res <- res[,grep('ID',colnames(res),invert=T)]
#   row.names(res) <- 1:nrow(res)
#   return(res)
# }

RegionStats <- function(object,genome,verbose = TRUE) {
  sequence.length <- width(x = object)
  sequences <- BSgenome::getSeq(x = genome, names = object)
  gc <- Biostrings::letterFrequency(
    x = sequences, letters = 'CG'
  ) / sequence.length * 100
  colnames(gc) <- 'GC.percent'
  dinuc <- Biostrings::dinucleotideFrequency(sequences)
  sequence.stats <- cbind(dinuc, gc, sequence.length)
  return(sequence.stats)
}

MatchRegionStats <- function(meta.feature,query.feature,features.match = c("GC.percent"),n = 10000,verbose = TRUE) {
  if (!inherits(x = meta.feature, what = 'data.frame')) {
    stop("meta.feature should be a data.frame")
  }
  if (!inherits(x = query.feature, what = "data.frame")) {
    stop("query.feature should be a data.frame")
  }
  if (length(x = features.match) == 0) {
    stop("Must supply at least one sequence characteristic to match")
  }
  if (nrow(x = meta.feature) < n) {
    n <- nrow(x = meta.feature)
    warning("Requested more features than present in supplied data.
            Returning ", n, " features")
  }
  # features.choose <- meta.feature[choosefrom, ]
  for (i in seq_along(along.with = features.match)) {
    featmatch <- features.match[[i]]
    if (!(featmatch %in% colnames(x = query.feature))) {
      if (i == "GC.percent") {
        stop("GC.percent not present in meta.features.",
             " Run RegionStats to compute GC.percent for each feature.")
      } else {
        stop(i, " not present in meta.features")
      }
    }
    if (verbose) {
      message("Matching ", featmatch, " distribution")
    }
    density.estimate <- density(
      x = query.feature[[featmatch]], kernel = "gaussian", bw = 1
    )
    weights <- approx(
      x = density.estimate$x,
      y = density.estimate$y,
      xout = meta.feature[[featmatch]],
      yright = 0.0001,
      yleft = 0.0001
    )$y
    if (i > 1) {
      feature.weights <- feature.weights * weights
    } else {
      feature.weights <- weights
    }
  }
  feature.select <- sample.int(
    n = nrow(x = meta.feature),
    size = n,
    prob = feature.weights
  )
  feature.select <- rownames(x = meta.feature)[feature.select]
  return(feature.select)
}

FindMotifs_fisher <- function(features,bg_features,background = 40000,genome_BS=genome,genome='hg38',pwm) {
  names(features) <- paste0(seqnames(features),'_',start(features),'_',end(features))
  names(bg_features) <- paste0(seqnames(bg_features),'_',start(bg_features),'_',end(bg_features))
  features$GC.percent <- as.data.frame(RegionStats(features,genome_BS))$GC.percent
  bg_features$GC.percent <- as.data.frame(RegionStats(bg_features,genome_BS))$GC.percent
  background <- MatchRegionStats(meta.feature = as.data.frame(bg_features),
                                 query.feature  = as.data.frame(features),
                                 n=background)
  motif_matrix <- CreateMotifMatrix(
    features = bg_features,
    pwm = pwm,
    genome = genome,
    sep = c("_", "_")
  )
  motif.names <- capitalize(tolower(as.vector(name(pwm))))
  names(motif.names) <- names(name(pwm))
  query.motifs <- motif_matrix[names(features), ]
  background.motifs <- motif_matrix[background, ]
  query.counts <- colSums(x = query.motifs)
  background.counts <- colSums(x = background.motifs)
  percent.observed <- query.counts / length(x = features) * 100
  percent.background <- background.counts / length(x = background) * 100
  fold.enrichment <- percent.observed / percent.background
  p.list <- c()
  for (i in seq_along(along.with = query.counts)) {
    dat <- data.frame(Background=c(background.counts[[i]],nrow(background.motifs) - background.counts[[i]]),Query=c(query.counts[[i]],length(features)-query.counts[[i]]))
    p.list[[i]] <- fisher.test(dat)$p.value
  }
  results <- data.frame(
    motif = names(x = query.counts),
    observed = query.counts,
    background = background.counts,
    percent.observed = percent.observed,
    percent.background = percent.background,
    fold.enrichment = fold.enrichment,
    pvalue = p.list,
    motif.name = as.vector(x = unlist(x = motif.names[names(x = query.counts)])),
    stringsAsFactors = FALSE
  )
  if (nrow(x = results) == 0) {
    return(results)
  } else {
    return(results[with(data = results, expr = order(pvalue, -fold.enrichment)), ])
  }
}

enrichGO_wrapper <- function(genes,qvalue.cutoff=0.05,orgDB=org.Mm.eg.db,organism='mmu'){
  gene.df <- bitr(genes, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = orgDB)
  ego <- enrichGO(gene         = unique(gene.df$ENTREZID),
                  OrgDb         = orgDB,
                  keyType       = 'ENTREZID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = qvalue.cutoff)
  p <- simplify(ego)
  return(p)
}

enrichedMotifs <- function(peaks,bg_peaks,genome,genome_BS,features=NULL,expr=NULL,cols=c('green','grey80','red'),logFC=log2(1.5),logP=5,point.size=2,anno.size=8,pwm=pwm){
  peaks <- read.table(peaks)[,1:3]
  colnames(peaks) <- c('chrom','start','end')
  bg_peaks <- read.table(bg_peaks)[,1:3]
  colnames(bg_peaks) <- c('chrom','start','end')
  background <- nrow(peaks)
  plot_list <- list()
  motif_list <- list()
  enriched.motifs <- FindMotifs_fisher(
    features = makeGRangesFromDataFrame(peaks),
    bg_features=makeGRangesFromDataFrame(unique(rbind(peaks,bg_peaks))),
    background=background,
    genome_BS=genome_BS,
    genome=genome,
    pwm
  )
  enriched.motifs$fold.enrichment <- log2(enriched.motifs$fold.enrichment)
  enriched.motifs$pvalue <- -log10(enriched.motifs$pvalue)
  enriched.motifs$col <- 'grey80'
  enriched.motifs$col[enriched.motifs$fold.enrichment>logFC&enriched.motifs$pvalue>=logP] <- 'red'
  enriched.motifs$col[enriched.motifs$fold.enrichment<(-logFC)&enriched.motifs$pvalue>=logP] <- 'green'
  enriched.motifs <- enriched.motifs[order(abs(enriched.motifs$fold.enrichment),enriched.motifs$pvalue,decreasing = T),]
  if (!is.null(expr)){
    tf_name <- sapply(strsplit(as.character(enriched.motifs$motif.name), "\\(|\\:|\\."), function(x) x[[1]])
    expr <- expr[match(tf_name,row.names(expr)),]
    enriched.motifs <- enriched.motifs[complete.cases(expr),]
  }
  p <- ggplot(enriched.motifs,aes(x=fold.enrichment,y=pvalue)) + geom_point(aes(fill = col), pch = I(21),size = point.size) 
  p <- p + scale_fill_manual(values = cols,labels=NULL) + xlab(expression(Log[2]~Fold~Change)) + ylab(expression(-Log[10]~(P)))
  if (is.null(features)){
    sel_df <- rbind(head(enriched.motifs[order(enriched.motifs$fold.enrichment,enriched.motifs$pvalue,decreasing=T),],10),head(enriched.motifs[order(enriched.motifs$fold.enrichment*(-1),enriched.motifs$pvalue,decreasing=T),],10))
    features <- sel_df$motif.name
  }
  p <- p + ggrepel::geom_text_repel(
    data = enriched.motifs, size = anno.size,seed = 42,
    box.padding =0.8, min.segment.length = 0,max.iter = 10000,
    aes(x=fold.enrichment,y=pvalue,color=NULL,label=ifelse(motif.name%in%features, as.character(motif.name), "")),force=10)
  enriched.motifs <- enriched.motifs[!is.infinite(enriched.motifs$fold.enrichment),]
  p <- p + theme(legend.position = "none") + xlim(c(-max(abs(enriched.motifs$fold.enrichment),na.rm=T),max(abs(enriched.motifs$fold.enrichment),na.rm=T)))
  return(list(p=p,enriched.motifs=enriched.motifs))
}

GOterm_enrichment<-function(distal_mat,promoter_bed,all_TSS,closest_TSS=FALSE,orgDB=org.Hs.eg.db) {
  promoter_mat <- read.table(promoter_bed)
  colnames(promoter_mat) <- c('chrom','start','end')
  promoter_mat <- gintervals.neighbors(promoter_mat,all_TSS)
  if(!closest_TSS){
    #mat <- mat[mat$domain=='intraTAD',]
    one_gene_per_enh <- ddply(distal_mat,.(distalIDs),function(x){
      score_maxs <- rowMaxs(as.matrix(x[,grep(c('score'),colnames(x))]),na.rm=T)
      score_maxs[is.infinite(score_maxs)] <- NA
      return(x[which(score_maxs==max(score_maxs,na.rm=T)),])
    },.parallel = T)
  } else {
    one_gene_per_enh <- ddply(distal_mat,.(distalIDs),function(x){
      return(x[which(abs(x$start2-x$start1)==min(abs(x$start2-x$start1))),])
    },.parallel = T)
  }
  gene.df <- bitr(c(as.character(one_gene_per_enh$geneName),as.character(promoter_mat$gene_name)), fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = orgDB)
  ego <- enrichGO(gene         = unique(gene.df$ENTREZID),
                  OrgDb         = orgDB,
                  keyType       = 'ENTREZID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  p <- simplify(ego)
  return(list(mat=one_gene_per_enh,p=p))
}

DMR_call<-function(con,files,genome,threshold,difference,min_Cov,qval,high_cov_perc,statistical_method,region,output_path){
  CpG_meth=methRead(location=files,
                    pipeline= "bismarkCoverage",
                    sample.id=list('3DRAM_Pax6_rep1','3DRAM_Pax6_rep2','3DRAM_Tbr2_rep1','3DRAM_Tbr2_rep2'),
                    assembly=genome,
                    treatment=c(0,0,1,1),
                    context=con)
  filtered_CpG_meth=filterByCoverage(CpG_meth,lo.count=min_Cov,hi.perc=high_cov_perc)
  united_CpG_meth=unite(filtered_CpG_meth, destrand=F)
  
  #Calculate counts at gNOME peaks and call differential regions
  tiled_CpG_meth=regionCounts(united_CpG_meth,region)
  #perform statistical test
  if (statistical_method=='default') {
    myDiff=calculateDiffMeth(tiled_CpG_meth,mc.cores=16)
  } else if (statistical_method=='over_disp_ftest') 
  {
    myDiff=calculateDiffMeth(tiled_CpG_meth,mc.cores=16,overdispersion="MN")
  } else if (statistical_method=='over_disp_chi') 
  {
    myDiff=calculateDiffMeth(tiled_CpG_meth,mc.cores=16,overdispersion="MN",test="Chisq")
  }
  #safe calls for later plotting
  write.table(as.data.frame(myDiff),file=paste0(output_path,con,'_all_',statistical_method,'.txt'),quote=F, sep='\t', row.names = F, col.names = T)
  myDiff_hyper=getMethylDiff(myDiff,difference=difference,qvalue=qval,type="hyper")
  export.bed(myDiff_hyper,con =paste0(output_path,con,'_hyper_',qval,'_',difference,'_',statistical_method,'.bed'))
  myDiff_hypo=getMethylDiff(myDiff,difference=difference,qvalue=qval,type="hypo")
  export.bed(myDiff_hypo,con =paste0(output_path,con,'_hypo_',qval,'_',difference,'_',statistical_method,'.bed'))
}

DMR_with_motifs<-function(chip_f,pwm_f,motif.names,motif_name,genome,cutoff,bed_dir) {
  chip <- read.table(chip_f)[,1:3]
  colnames(chip) <- c('chrom','start','end')
  df <- makeGRangesFromDataFrame(chip)
  motif_pos <- unlist(matchMotifs(pwm_f[[which(motif.names==motif_name)[1]]], df, genome = genome, out = "positions",p.cutoff = cutoff )) #added the which(motif.names==motif.name)[1] sinc ePax6 is double in the rds
  motif_6 <- data.frame(seqnames=seqnames(motif_pos),starts=start(motif_pos),ends=end(motif_pos),names=c(rep(".", length(motif_pos))),scores=elementMetadata(motif_pos[,1]),strands=strand(motif_pos))
  motif_name<-gsub('(','',motif_name, fixed = T)
  motif_name<-gsub(')','',motif_name, fixed = T)
  write.table(motif_6,paste0(bed_dir,motif_name,'_',cutoff,'.bed'),quote=F,sep='\t',col.names=F,row.names=F)
  command_f<-paste0('bedtools intersect -a ',chip_f,' -b ',paste0(bed_dir,motif_name,'_',cutoff,'.bed'),' -u > ',paste0(bed_dir,motif_name,'_',cutoff,'_',motif_name,'.bed && bedtools intersect -a ',chip_f,' -b ',paste0(bed_dir,motif_name,'_',cutoff,'.bed'),' -v > ',paste0(bed_dir,motif_name,'_',cutoff,'_no',motif_name,'.bed')))
  system(command_f)
}

repeat_analysis <- function(repeat_bed,interval_window,min_dist,max_dist,expand,domains_f,score_tracks,hic_names,other_tracks,other_names){
  domains <- gintervals.load(domains_f)
  dom_2d = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
  
  for (i in 1:length(score_tracks)){
    gvtrack.create(paste0('v_',score_tracks[i]),score_tracks[i],'max')
    gvtrack.iterator.2d(paste0('v_',score_tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
  }
  
  for (i in 1:length(other_tracks)){
    gvtrack.create(paste0('v_',other_tracks[i]),other_tracks[i],'avg')
    gvtrack.iterator(paste0('v_',other_tracks[i]), sshift=-interval_window/2, eshift=interval_window/2)
  }
  
  ########################
  repeat_bed <- gintervals.canonic(repeat_bed)
  grid <- unique(construct.grid(repeat_bed,repeat_bed,min_dist,max_dist))
  grid <- gintervals.canonic(grid)
  
  intra_grid <- gintervals.intersect(grid,dom_2d)
  inter_grid <- grid[!(paste0(grid[,1],':',grid[,2],'-',grid[,4],':',grid[,5])%in%paste0(intra_grid[,1],':',intra_grid[,2],'-',intra_grid[,4],':',intra_grid[,5])),]

  intra_scores<- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = intra_grid,iterator = intra_grid,band = -c(max_dist+max(expand),min_dist-max(expand)))
  inter_scores <- gextract(paste0('v_',score_tracks),colnames = hic_names,intervals = inter_grid,iterator = inter_grid,band = -c(max_dist+max(expand),min_dist-max(expand)))
  
  intra_scores$domain <- 'intraTAD'
  inter_scores$domain <- 'interTAD'
  scores <- rbind(intra_scores,inter_scores)

  distal_scores <- gextract(paste0('v_',other_tracks),colnames=other_names,intervals = repeat_bed,iterator = repeat_bed)
  distal_scores[,grep('ins',colnames(distal_scores))] <- distal_scores[,grep('ins',colnames(distal_scores))]*(-1)
  
  return(list(hicscores=scores,lin_scores=distal_scores))
}

featureToGR <- function(feature,pattern="_"){
  featureSplit <- stringr::str_split(paste0(feature), pattern =pattern , n = 3, simplify = TRUE)
  gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
  return(gr)
}

motifsPerCluster <- function(df,pwms,cre_size,genome,mcparams,FDR.cutoff,enr.cutoff,out_f,height,width){
  pwms <- toPWM(pwms)
  pwms <- pwms[!duplicated(names(pwms))]
  peaks <- makeGRangesFromDataFrame(df[,2:4])
  names(peaks) <- df$intervalID
  peaks <- resize(peaks,cre_size,fix='center')          #Change coords to actual MPRA coords
  peakseqs <- getSeq(genome, peaks)
  se <- calcBinnedMotifEnr(seqs = peakseqs, bins = factor(df$cluster), motifs = pwms,BPPARAM = mcparams)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  seSel <- se[sel1&sel2, ]
  hcl <- hclust(dist(assay(seSel,"log2enr")), method = "ward.D")
  pdf(out_f,height=height,width=width)
  plotMotifHeatmaps(x = seSel, which.plots = c("log2enr"), width = 1.8,cluster=hcl,
                    show_dendrogram = TRUE, show_seqlogo = TRUE,
                    width.seqlogo = 1.2)
  dev.off()
  return(se)
}

associate_peaks_with_summit <- function(peak_f,summit_f='results/SC/macs2/union_summits.bed',out_f){
  test <- read.table(peak_f)
  test <- test[,1:3]
  colnames(test)[1:3] <- c('chrom','start','end')
  test2 <- read.table(summit_f)
  colnames(test2)[1:3] <- c('chrom','start','end')
  test2 <- test2[test2$chrom%in%gintervals.all()$chrom,]
  df <- gintervals.neighbors(test2,test)
  df$peakID <- paste0(df[,6],df[,7])
  df <- df[df$dist==0,]
  df <- df[order(df$V5,decreasing=T),]
  df <- df[!duplicated(df$peakID),]
  df <- df[order(df[,1],df[,2]),1:3]
  write.table(df,out_f,quote=F,col.names=F,row.names=F,sep='\t')
  
}

extract_MotifPos <- function(pwm,chip_f=NULL,motif.name,genome){
  motif.names <- as.vector(name(pwm))
  if(!is.null(chip_f)){
    chip <- read.table(chip_f)[,1:3]
    colnames(chip) <- c('chrom','start','end')
    df <- makeGRangesFromDataFrame(chip)
  }
  motif_pos <- unlist(matchMotifs(pwm[[which(motif.names==motif.name)]], df, genome = genome, out = "positions",bg='genome') )
  motif_6 <- data.frame(seqnames=seqnames(motif_pos),starts=start(motif_pos),ends=end(motif_pos),names=c(rep(".", length(motif_pos))),scores=elementMetadata(motif_pos[,1]),strands=strand(motif_pos))
  return(motif_6)
}


cor_ChIP_ATAC <- function(peaks_f,tracks,window_f=200,track_names){
  peaks <- read.table(peaks_f)[,1:3]
  colnames(peaks) <- c('chrom','start','end')
  if(!is.null(window_f)){
    peaks <- intervals.normalize(peaks,window_f)
    peaks <- gintervals.canonic(peaks)
  } else {
    peaks <- peaks[!duplicated(peaks),]
  }
  df <- misha_extract(tracks,regions=peaks,track_names=track_names,window = window_f,iterator = peaks)
  df <- df[,-c(1:3)]
  p <- ggplot(df, aes_( x = as.name(track_names[2]), y = as.name(track_names[1])  )) + geom_pointdensity(shape=19,size=1,alpha=1) + scale_color_gradientn(colours = c('darkblue','red','orange'),name='')
  p <- p + coord_cartesian(xlim=c(0,quantile(df[,2],0.999,na.rm=T)),ylim=c(0,quantile(df[,1],0.999,na.rm=T))) + theme(legend.position = "none") + annotate("text",x = 4.2, y = 0.05, label = paste0("r=",round(cor(df[,1],df[,2],use = 'complete.obs',method='pearson'),2)),size=5)
  return(p)
}

stratifyPeaks_by_Motifs <- function(peaks_f,pwm_f='data/combined_pwm_FPKM1.RDS',motif='NEUROG2(var.2)',min.score=10,stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,control_regions,genome,mcparams=BiocParallel::MulticoreParam(10L)){
  pwms <- readRDS(pwm_f)
  pwms <- toPWM(pwms)
  peaks <- import(peaks_f)
  #mcols(peaks) <- NULL
  peaks <- peaks[!duplicated(peaks),]
  peakseqs <- getSeq(genome, peaks)
  names(peakseqs) <- 1:length(peakseqs)
  res <- findMotifHits(query = pwms[which(name(pwms)==motif),],
                       subject = peakseqs,
                       min.score = min.score,
                       method = "matchPWM",
                       BPPARAM = mcparams)
  
  
  ###
  
  res <- table(factor(seqnames(res), levels = names(peakseqs)),
               factor(res$pwmname, levels = name(pwms)))
  if(!is.null(stratify_by)){
    peaks_m <- as.data.frame(peaks)[,-c(4:6)]
    colnames(peaks_m)[1] <- 'chrom'
    df <- misha_extract(tracks=stratify_by,regions=peaks_m,track_names=stratify_by,window = window_f,iterator = peaks_m)
    df$rank <- ntile(df[,4],n_quantile)
    mat <- data.frame(ID=as.character(names(peakseqs)),nMotif=as.numeric(res[,motif]),rank=as.factor(df$rank),type=type_names[peaks_m$score])
    
    c_peaks <- import(control_regions)
    mcols(c_peaks) <- NULL
    c_peaks <- c_peaks[!duplicated(c_peaks)&seqnames(c_peaks)%in%df$chrom,]
    c_peaks_m <- as.data.frame(c_peaks)[,1:3]
    colnames(c_peaks_m)[1] <- 'chrom'
    peaks_overlap <- gintervals.neighbors(c_peaks_m,peaks_m)
    c_peaks <- c_peaks[peaks_overlap$dist>=1000]
    c_peaks <- c_peaks[sample(1:length(c_peaks),nrow(mat[mat$rank==1,]))]
    seqlevels(c_peaks) <- seqlevelsInUse(c_peaks)
    c_peakseqs <- getSeq(genome, c_peaks)
    names(c_peakseqs) <- 1:length(c_peakseqs)
    c_res <- findMotifHits(query = pwms[which(name(pwms)==motif),],
                           subject = c_peakseqs,
                           min.score = min.score,
                           method = "matchPWM",
                           BPPARAM = mcparams)
    c_res <- table(factor(seqnames(c_res), levels = names(c_peakseqs)),
                   factor(c_res$pwmname, levels = name(pwms)))
    c_mat <- data.frame(ID=as.character(names(c_peakseqs)),nMotif=as.numeric(c_res[,motif]),rank=rep("C",nrow(c_res)),type=sample(type_names,size=nrow(c_res),replace = T,prob = table(mat$type)/sum(table(mat$type))))
    mat <- rbind(mat,c_mat)
    mat$ID <- as.character(1:nrow(mat))
    mat$rank <- factor(mat$rank,levels=c("C",1:n_quantile))
    mat$type <- factor(mat$type,levels=type_names)
    mat_m <- melt(mat,id.vars = c('ID','rank','type'))
    mat_m$value[mat_m$value>=4] <- 4
    mat_m$value <- factor(mat_m$value,levels=c(0:4))
    levels(mat_m$value) <- c('0','1','2','3', '4+')
    p <- ggplot(mat_m, aes(x=rank,y=1,fill=value)) +  geom_bar(position="fill", stat="identity")     #+ scale_fill_manual(name='cluster',values = as.character(cols))
    p <- p + ylab('Percentage of Total') + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('')  + theme(legend.title = element_blank())
    return(p)
  } else {
    mcols(peaks) <- as.data.frame.matrix(res)
    return(peaks)
  }
}


stratifyPeaks_by_Motifs_extract <- function(peaks_f,pwm_f='data/combined_pwm_FPKM1.RDS',motif.name='NEUROG2(var.2)',tracks='scATAC.repro_AST',track_names,window_f=200,max_motifs=4,genome,mcparams=BiocParallel::MulticoreParam(10L),cols=sample_colors){
  pwms <- readRDS(pwm_f)
  pwms <- toPWM(pwms)
  peaks <- import(peaks_f)
  mcols(peaks) <- NULL
  peaks <- peaks[!duplicated(peaks),]
  peakseqs <- getSeq(genome, peaks)
  names(peakseqs) <- 1:length(peakseqs)
  motif.names <- as.vector(name(pwms))
  res <- findMotifHits(query = pwms[[which(motif.names==motif.name)]],
                       subject = peakseqs,
                       min.score = 6.0,
                       method = "matchPWM",
                       BPPARAM = mcparams)
  #####
  # df <- as.data.frame(resize(res,width = 1,fix='center'))
  # test <- dlply(df,.(seqnames),function(x){
  #   x_ranges <- x$start
  #   x_ranges <- sort(dist(x_ranges[order(x_ranges)]))
  #   #test <- c(test,x_ranges)
  #   return(x_ranges)
  # })
  # test <- as.vector(unlist(test))
  # 
  ###
  
  res <- table(factor(seqnames(res), levels = names(peakseqs)),
               factor(res$pwmname, levels = name(pwms[[which(motif.names==motif.name)]])))
  if(!is.null(tracks)){
    peaks_m <- as.data.frame(peaks)[,1:3]
    colnames(peaks_m)[1] <- 'chrom'
    df <- misha_extract(tracks=tracks,regions=peaks_m,track_names=track_names,window = window_f,iterator = peaks_m)
    df$n_motif <- res[,1]
    df$n_motif[df$n_motif>max_motifs] <- max_motifs
    mat <- df[,c('intervalID',track_names,'n_motif')]
    mat_m <- melt(mat,id.vars = c('intervalID','n_motif'))
    mat_m$n_motif <- factor(mat_m$n_motif,levels=c(0:max_motifs))
    p <- ggplot(mat_m, aes(x=n_motif,y=value,fill=variable)) + geom_boxplot(outlier.size=1,outlier.shape = NA,width=0.8)     #+ scale_fill_manual(name='cluster',values = as.character(cols))
    p <- p  + scale_fill_manual(name='',values=cols)
    return(p)
  } else {
    mcols(peaks) <- as.data.frame.matrix(res)
    return(peaks)
  }
}

de_novo_Motifs <- function(peaks_f1,peaks_f2,pwm_f='data/combined_pwm_FPKM1.RDS',genome,mcparams=BiocParallel::MulticoreParam(10L),FDR.cutoff=8,enr.cutoff=0.25){
  pwms <- readRDS(pwm_f)
  pwms <- toPWM(pwms)
  peaks1 <- import(peaks_f1)
  mcols(peaks1) <- NULL
  peaks1 <- peaks1[!duplicated(peaks1),]
  if(!is.null(peaks2)){
    peaks2 <- import(peaks_f2)
    mcols(peaks2) <- NULL
    peaks2 <- peaks2[!duplicated(peaks2),]
    peaks <- c(peaks1,peaks2)
    bins <- as.factor(c(rep('Ngn2',length(peaks1)),rep('PmutNgn2',length(peaks2))))
  } else {
    peaks <- peaks1
    bins <- NULL
  }
  peakseqs <- getSeq(genome, peaks)
  names(peakseqs) <- 1:length(peakseqs)
  
  se <- calcBinnedMotifEnrR(seqs = peakseqs,bins = bins,pwmL = pwms,background = 'genome',genome = BSgenome.Mmusculus.UCSC.mm10,BPPARAM = mcparams)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  seSel <- se[sel1&sel2, ]
  
  sekm <- calcBinnedKmerEnr(seqs = peakseqs,bins = bins,kmerLen = 8,background = 'genome',genome = BSgenome.Mmusculus.UCSC.mm10,BPPARAM = mcparams)
  selkm1 <- apply(assay(sekm, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  selkm2 <- apply(assay(sekm, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  sekmSel <- sekm[selkm1&selkm2, ]
  pfmSel <- rowData(seSel)$motif.pfm
  sims <- motifKmerSimilarity(x = pfmSel,
                              kmers = rownames(sekmSel),
                              includeRevComp = TRUE)
  dim(sims)
  
  maxwidth <- max(sapply(TFBSTools::Matrix(pfmSel), ncol))
  seqlogoGrobs <- lapply(pfmSel, seqLogoGrob, xmax = maxwidth)
  hmSeqlogo <- rowAnnotation(logo = annoSeqlogo(seqlogoGrobs, which = "row"),
                             annotation_width = unit(1.5, "inch"), 
                             show_annotation_name = FALSE)
  hm <- Heatmap(sims, 
                show_row_names = TRUE, row_names_gp = gpar(fontsize = 8),
                show_column_names = TRUE, column_names_gp = gpar(fontsize = 8),
                name = "Similarity", column_title = "Selected TFs and enriched k-mers",
                col = colorRamp2(c(0, 1), c("white", "red")), 
                right_annotation = hmSeqlogo)
  
}


motif_footprinting <- function(archr_obj=archr_obj,nTop=NULL,flank=500,motif=NULL,regions=NULL,normalize_by='Subtract',out_f,group.by,which_clusters,cols=cluster_colors(length(levels(atac.object))),plot_bias=T,anno.size=12,key.size=4,width,height,returnPlot=F){
  if (is.character(regions)){
    chip <- read.table(regions)
    colnames(chip)[1:3] <- c('chrom','start','end')
    regions <- makeGRangesFromDataFrame(chip)
  }
  if(!is.null(motif)){
    motifPositions <- getPositions(archr_obj)
    names(motifPositions) <- names(motifPositions)
    name_arch <- names(motifPositions)[grep(gsub('\\(|\\:|\\)','.',motif),names(motifPositions))]
    positions <- motifPositions[name_arch,]
    if(!is.null(regions)){
      subset_positions <- GRangesList(subsetByOverlaps(positions[[1]],regions))
      names(subset_positions) <- names(positions)[1]
      positions <- subset_positions
    }
  } else {
    positions <- GRangesList(regions)
    name_arch <- NULL
  }
  message(length(positions[[1]]))
  seFoot <- suppressMessages(getFootprints(
    ArchRProj = archr_obj, flank = flank,nTop = nTop,
    positions = positions, 
    groupBy = group.by,useGroups = which_clusters,threads = 1,verbose = F,logFile = '/dev/null'
  ))
  p <- ggFootprint(
    seFoot = seFoot,
    name = name_arch[1],
    pal = cols,
    smoothWindow = 50,
    flank = flank,
    flankNorm = 50,
    baseSize = 10,
    normMethod = 'Subtract',
    which_clusters=which_clusters
  )
  
  p1 <- p$ggFoot + theme(legend.text = element_text(size=anno.size),legend.justification=c(0,0), legend.position=c(0.80, 0.75), legend.box = "vertical",legend.spacing.x = unit(0.1, 'inch')) + guides(colour=guide_legend(override.aes = list(size = key.size),ncol = 1)) + ggtitle(ifelse(is.null(motif),'',motif))
  p2 <- p$ggBias + theme(legend.position="none") + ylab('') + scale_y_continuous(n.breaks=3) + xlab('')
  if (plot_bias){
    p <- p1+ p2 + plot_layout(nrow=2,heights=c(4,1))
  } else {
    p <- p1
  }
  if(returnPlot){
    return(p)
  } else {
    pdf(paste0('plots/figures/',out_f,'.pdf'),height=height,width=width)
    print(p)
    dev.off()
  }
}