#### If Please cite as Noack et al., Nature Neuroscience 2021 """
#### Created by Boyan Bonev ####


main_f <- '/home/hpc/bonev/projects/repro_dfg/'        #Please set this to the main path where you wish to perform the analysis
library(pals)
### Global Directories ###
sample_dir <- '/home/hpc/bonev/projects/SC/analysis/scATAC/'           # the files from the 10x scATAC run should be located here
sample_names <- c('Astro','GFP','Ngn2','PmutNgn2')

library(Hmisc)
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

### Colors
#archr_colors <- ArchR::ArchRPalettes$stallion
#sample_colors <-  function(x){return(brewer.pal(n = x, name = "Reds"))}
sample_colors <- c("#18C0C9","#215801","#E6CA17","#D51F26")
cluster_colors <- c("#272E6A","#208A42","#FA7517","#A30005","#89288F")
yy1_colors <- c("#89C75F","#8A9FD1","#7E1416","#D8A767")
yy1_clusters <- c("#272E6A","#208A42","#FEE500","#FA7517","#A30005","#89288F","#90D5E4","#0C727C","#3D3D3D")
gene_colors <- function(x){colorRampPalette(rev(colorpalette('reds',10)))(x)}
#cluster_colors <- function(x){cols25(x)}
heatmap_colors <- colorpalette('matlablike',10)
col.pbreaks <<- c(20,35,50,65,75,85,95)        #Original
col.pos <<- palette.breaks(100 , c("lavenderblush2","#f8bfda","lightcoral","red","orange","yellow"), col.pbreaks)
col.nbreaks <<- c(20,35,50,65,75,85,95)
col.neg <<- rev(palette.breaks(100 , c("powderblue", "cornflowerblue", "blue","blueviolet", "#8A2BE2", "#4B0082"), col.nbreaks ))
hic.scores <<- c(col.neg, "lightgrey", col.pos)

### Marker Genes
diff_markers <- c('Hes1','Neurog2','Eomes','Dcx','Mapt','S100b','Tead3','Rbfox3','Tubb3')
layer_markers <- c('Gfap','Apoe','Sox9','Sox11','Vim','Neurod2','Id4')
other_markers <- c('C1qc','Vtn','Reln','Gad2','Gad1','Vip','Sst')
mitotic_markers <- c('Cdk1','Ube2c','Top2a','Hist1h4e','Mki67','Pcna')
tf_markers <- c('Id4','Id1','Tcf4','Hey1','Yy1','Zeb1','Meis2','Smarca2','Tcf12','Ldb1','Lhx2','Tgif2','Lmo3','Lmo4','Isl1','Nr2f1','Smarca5','Chd4')
repro_markers <- c('Dcc','Plxna4','Tubb3','Neurog2','Sema5a')
repro_diff_markers <- c('Epha3','Robo1','Zeb2','Reln','Clstn2')
repro_diff_markers2 <- c('Gfap','Ephx1','Zeb2','Reln','Clstn2')
### Genomic databases ###

#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
