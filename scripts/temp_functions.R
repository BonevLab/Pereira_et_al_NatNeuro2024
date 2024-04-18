
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

stratifyPeaks_by_Motifs <- function(peaks_f,pwm_f='data/combined_pwm_FPKM1.RDS',motif='NEUROG2(var.2)',stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,control_regions,genome,mcparams=BiocParallel::MulticoreParam(10L)){
  pwms <- readRDS(pwm_f)
  pwms <- toPWM(pwms)
  peaks <- import(peaks_f)
  mcols(peaks) <- NULL
  peaks <- peaks[!duplicated(peaks),]
  peakseqs <- getSeq(genome, peaks)
  names(peakseqs) <- 1:length(peakseqs)
  res <- findMotifHits(query = pwms,
                       subject = peakseqs,
                       min.score = 12.0,
                       method = "matchPWM",
                       BPPARAM = mcparams)
  #####
  df <- as.data.frame(res)
  test <- ddply(df,.(seqnames),function(x){
    
  })
  
  
  ###
  
  res <- table(factor(seqnames(res), levels = names(peakseqs)),
               factor(res$pwmname, levels = name(pwms)))
  if(!is.null(stratify_by)){
    peaks_m <- as.data.frame(peaks)[,1:3]
    colnames(peaks_m)[1] <- 'chrom'
    df <- misha_extract(tracks=stratify_by,regions=peaks_m,track_names=stratify_by,window = window_f,iterator = peaks_m)
    df$rank <- ntile(df[,4],n_quantile)
    mat <- data.frame(ID=as.character(names(peakseqs)),nMotif=as.numeric(res[,motif]),rank=as.factor(df$rank))
    
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
    c_res <- findMotifHits(query = pwms,
                           subject = c_peakseqs,
                           min.score = 6.0,
                           method = "matchPWM",
                           BPPARAM = mcparams)
    c_res <- table(factor(seqnames(c_res), levels = names(c_peakseqs)),
                   factor(c_res$pwmname, levels = name(pwms)))
    c_mat <- data.frame(ID=as.character(names(c_peakseqs)),nMotif=as.numeric(c_res[,motif]),rank=rep("C",nrow(c_res)))
    mat <- rbind(mat,c_mat)
    mat$ID <- as.character(1:nrow(mat))
    mat$rank <- factor(mat$rank,levels=c("C",1:n_quantile))
    mat_m <- melt(mat,id.vars = c('ID','rank'))
    mat_m$value[mat_m$value>=4] <- 4
    mat_m$value <- factor(mat_m$value,levels=c(0:4))
    levels(mat_m$value) <- c('No motifs','1 motif','2 motifs','3 motifs', '4+ motifs')
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
  res <- findMotifHits(query = pwms,
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

p1 <- cor_ChIP_ATAC(peaks_f='data/ChIP/conservative_peaks.narrowPeak',tracks=c('chipseq_RPM.iN_Ngn2_D2','scATAC.repro_AST'),window_f=200,track_names=c('Ngn2_ChIP','AST_ATAC'))
p2 <- cor_ChIP_ATAC(peaks_f='data/ChIP/conservative_peaks.narrowPeak',tracks=c('chipseq_RPM.iN_Ngn2_D2','scATAC.repro_iN_1'),window_f=200,track_names=c('Ngn2_ChIP','iN1_ATAC'))
p3 <- cor_ChIP_ATAC(peaks_f='data/ChIP/conservative_peaks.narrowPeak',tracks=c('chipseq_RPM.iN_Ngn2_D2','scATAC.repro_iN_2'),window_f=200,track_names=c('Ngn2_ChIP','iN2_ATAC'))

p <- p1+p2+p3 + plot_layout(guides = "collect")
pdf('plots/temp/Ngn2ChIPvsATAC_Ngn2Peaks.pdf',height=4,width=10)
p
dev.off()

p1 <- cor_ChIP_ATAC(peaks_f='results/beds/unionATAC_peaks_Ngn2Motifs.bed',tracks=c('chipseq_RPM.iN_Ngn2_D2','scATAC.repro_AST'),window_f=200,track_names=c('Ngn2_ChIP','AST_ATAC'))
p2 <- cor_ChIP_ATAC(peaks_f='results/beds/unionATAC_peaks_Ngn2Motifs.bed',tracks=c('chipseq_RPM.iN_Ngn2_D2','scATAC.repro_iN_1'),window_f=200,track_names=c('Ngn2_ChIP','iN1_ATAC'))
p3 <- cor_ChIP_ATAC(peaks_f='results/beds/unionATAC_peaks_Ngn2Motifs.bed',tracks=c('chipseq_RPM.iN_Ngn2_D2','scATAC.repro_iN_2'),window_f=200,track_names=c('Ngn2_ChIP','iN2_ATAC'))
p <- p1+p2+p3 + plot_layout(guides = "collect")
pdf('plots/temp/Ngn2ChIPvsATAC_Ngn2Motifs.pdf',height=4,width=10)
p
dev.off()


p1 <- cor_ChIP_ATAC(peaks_f='ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak',tracks=c('chipseq_RPM.iN_PmutNgn2_D2','scATAC.repro_AST'),window_f=200,track_names=c('PmutNgn2_ChIP','AST_ATAC'))
p2 <- cor_ChIP_ATAC(peaks_f='ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak',tracks=c('chipseq_RPM.iN_PmutNgn2_D2','scATAC.repro_iN_1'),window_f=200,track_names=c('PmutNgn2_ChIP','iN1_ATAC'))
p3 <- cor_ChIP_ATAC(peaks_f='ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak',tracks=c('chipseq_RPM.iN_PmutNgn2_D2','scATAC.repro_iN_2'),window_f=200,track_names=c('PmutNgn2_ChIP','iN2_ATAC'))
p <- p1+p2+p3 + plot_layout(guides = "collect")
pdf('plots/temp/PmutNgn2ChIPvsATAC_PmutNgn2Peaks.pdf',height=4,width=10)
p
dev.off()

p1 <- cor_ChIP_ATAC(peaks_f='results/beds/unionATAC_peaks_Ngn2Motifs.bed',tracks=c('chipseq_RPM.iN_PmutNgn2_D2','scATAC.repro_AST'),window_f=200,track_names=c('Ngn2_ChIP','AST_ATAC'))
p2 <- cor_ChIP_ATAC(peaks_f='results/beds/unionATAC_peaks_Ngn2Motifs.bed',tracks=c('chipseq_RPM.iN_PmutNgn2_D2','scATAC.repro_iN_1'),window_f=200,track_names=c('Ngn2_ChIP','iN1_ATAC'))
p3 <- cor_ChIP_ATAC(peaks_f='results/beds/unionATAC_peaks_Ngn2Motifs.bed',tracks=c('chipseq_RPM.iN_PmutNgn2_D2','scATAC.repro_iN_2'),window_f=200,track_names=c('Ngn2_ChIP','iN2_ATAC'))
p <- p1+p2+p3 + plot_layout(guides = "collect")
pdf('plots/temp/PmutNgn2ChIPvsATAC_Ngn2Motifs.pdf',height=4,width=10)
p
dev.off()

p <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Data/Peaks/Ngn2_conservative_peaks.narrowPeak',stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,
                                    control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/temp/Ngn2ChIP_Ngn2Motifs_byAST_ATACranks.pdf',height=4,width=5)
print(p+ggtitle("Ngn2 ChIP - AST ATAC"))
dev.off()

p <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Data/Peaks/Ngn2_conservative_peaks.narrowPeak',stratify_by='scATAC.repro_iN_1',window_f=200,n_quantile=4,
                             control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/temp/Ngn2ChIP_Ngn2Motifs_by_iN1_ATACranks.pdf',height=4,width=5)
print(p+ggtitle("Ngn2 ChIP - iN1 ATAC"))
dev.off()

p <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Data/Peaks/Ngn2_conservative_peaks.narrowPeak',stratify_by='scATAC.repro_iN_2',window_f=200,n_quantile=4,
                             control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/temp/Ngn2ChIP_Ngn2Motifs_by_iN2_ATACranks.pdf',height=4,width=5)
print(p+ggtitle("Ngn2 ChIP - iN2 ATAC"))
dev.off()


p <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Analysis/results_BB/shared_peaks.bed',stratify_by='scATAC.repro_AST',window_f=200,n_quantile=4,
                             control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/temp/Ngn2specific_by_AST_ATACranks.pdf',height=4,width=5)
print(p+ggtitle("Ngn2 specific Peaks - AST ATAC"))
dev.off()

p <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak',stratify_by='scATAC.repro_iN_1',window_f=200,n_quantile=4,
                             control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/temp/PmutNgn2ChIP_Ngn2Motifs_by_iN1_ATACranks.pdf',height=4,width=5)
print(p+ggtitle("PmutNgn2 ChIP - iN1 ATAC"))
dev.off()

p <- stratifyPeaks_by_Motifs(peaks_f='ChIP_Nov22/Data/Peaks/PmutNgn2_conservative_peaks.narrowPeak',stratify_by='scATAC.repro_iN_2',window_f=200,n_quantile=4,
                             control_regions='results/SC/macs2/union_peaks.narrowPeak',genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L))

pdf('plots/temp/PmutNgn2ChIP_Ngn2Motifs_by_iN2_ATACranks.pdf',height=4,width=5)
print(p+ggtitle("PmutNgn2 ChIP - iN2 ATAC"))
dev.off()


pdf('plots/temp/PmutNgn2_deNovo_Motifs.pdf',height=6,width=24)
draw(hm)
dev.off()


### 

p <- stratifyPeaks_by_Motifs_extract(peaks_f='ChIP_Nov22/Analysis/results_BB/Ngn2_distalPeaks5000.bed',tracks=c("scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2"),track_names=c('GFP','Ngn2','PmutNgn2'),
                                     window_f=200,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),cols = sample_colors[2:4])
pdf('plots/temp/Ngn2peaks5000d_accessibility_Ngn2motifs.pdf',height=4,width=4.5)
p + theme(legend.position = c(0.7,0.95)) + coord_cartesian(ylim=c(0,3)) + ylab('Accessibility') + xlab('Number of Neurog2 motifs') + ggtitle('Ngn2 peaks')
dev.off()

p <- stratifyPeaks_by_Motifs_extract(peaks_f='ChIP_Nov22/Analysis/results_BB/shared_distalPeaks5000.bed',tracks=c("scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2"),track_names=c('GFP','Ngn2','PmutNgn2'),
                                     window_f=200,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),cols = sample_colors[2:4])
pdf('plots/temp/Sharedpeaks5000d_accessibility_Ngn2motifs.pdf',height=4,width=4.5)
p + theme(legend.position = c(0.7,0.95)) + coord_cartesian(ylim=c(0,3)) + ylab('Accessibility') + xlab('Number of Neurog2 motifs') + ggtitle('Shared peaks')
dev.off()

p <- stratifyPeaks_by_Motifs_extract(peaks_f='ChIP_Nov22/Analysis/results_BB/PmutNgn2_distalPeaks5000.bed',tracks=c("scATAC.repro_GFP","scATAC.repro_Ngn2","scATAC.repro_PmutNgn2"),track_names=c('GFP','Ngn2','PmutNgn2'),max_motifs = 5,
                                     window_f=200,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),cols = sample_colors[2:4])
pdf('plots/temp/PmutNgn2peaks5000d_accessibility_Ngn2motifs.pdf',height=4,width=4.5)
p + theme(legend.position = c(0.7,0.95)) + coord_cartesian(ylim=c(0,3)) + ylab('Accessibility') + xlab('Number of Neurog2 motifs') + ggtitle('PmutNgn2 peaks')
dev.off()

p <- stratifyPeaks_by_Motifs_extract(peaks_f='ChIP_Nov22/Analysis/results_BB/Ngn2_distalPeaks5000.bed',tracks=c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x"),track_names=c('GFP','Ngn2','PmutNgn2'),
                                     window_f=200,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),cols = sample_colors[2:4])
pdf('plots/temp/Ngn2peaks5000d_CpG_Ngn2motifs.pdf',height=4,width=4.5)
p + theme(legend.position = c(0.07,0.85)) + ylab('CpG(%)') + xlab('Number of Neurog2 motifs') + ggtitle('Ngn2 peaks')
dev.off()

p <- stratifyPeaks_by_Motifs_extract(peaks_f='ChIP_Nov22/Analysis/results_BB/shared_distalPeaks5000.bed',tracks=c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x"),track_names=c('GFP','Ngn2','PmutNgn2'),
                                     window_f=200,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),cols = sample_colors[2:4])
pdf('plots/temp/Sharedpeaks5000d_CpG_Ngn2motifs.pdf',height=4,width=4.5)
p + theme(legend.position = c(0.07,0.85)) + ylab('CpG(%)') + xlab('Number of Neurog2 motifs') + ggtitle('Shared peaks')
dev.off()

p <- stratifyPeaks_by_Motifs_extract(peaks_f='ChIP_Nov22/Analysis/results_BB/PmutNgn2_distalPeaks5000.bed',tracks=c("methylation.JD_GFP_D2_G0G1_CpG_10x","methylation.JD_Ngn2_D2_G0G1_CpG_10x","methylation.JD_pmutNgn2_D2_G0G1_CpG_10x"),track_names=c('GFP','Ngn2','PmutNgn2'),
                                     window_f=200,genome=BSgenome.Mmusculus.UCSC.mm10,mcparams=BiocParallel::MulticoreParam(10L),cols = sample_colors[2:4])
pdf('plots/temp/PmutNgn2peaks5000d_CpG_Ngn2motifs.pdf',height=4,width=4.5)
p + theme(legend.position = c(0.07,0.85)) + ylab('CpG(%)') + xlab('Number of Neurog2 motifs') + ggtitle('PmutNgn2 peaks')
dev.off()



