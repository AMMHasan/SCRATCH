# Submit: qsub -N CASCADE -l h_rt=48:00:00 -l mem=40G -l tmpfs=20G -wd /lustre/scratch/scratch/regmpcr/stratosphere/data/LOG/ /lustre/scratch/scratch/regmpcr/bash-qsub/analysis_CASCADE_byPatient.sh
# R --silent --slave --vanilla < /lustre/scratch/scratch/regmpcr/stratosphere/scripts/analysis_CASCADE_byPatient.R -r
# /lustre/scratch/scratch/regmpcr/bash-qsub/analysis_CASCADE_byPatient.sh > cn_CASCADE_byPatient.txt 2>&1 &
###################################################################################################
### Load Libraries
###################################################################################################
require(ACE)
require(arm)
require(arules)
require(arulesViz)
require(basetheme)
require(Biobase)
require(BSgenome.Hsapiens.UCSC.hg19)
require(CGHbase)
require(CGHcall)
require(circlize)
require(CGPfunctions)
require(Clonality)
require(CloneDeMix)
require(cluster)   
require(ComplexHeatmap)
require(cowplot)
require(dendextend)
require(dendsort)
require(DNAcopy)
require(dplyr)
require(FactoMineR) 
require(factoextra) 
require(ggnewscale)
require(ggplot2)
require(ggridges)
require(ggrepel)
require(gplots)
require(ggpubr)
require(ggpmisc)
require(ggsignif)
require(grDevices)
require(gridExtra)
require(Hmisc)
require(imputeTS)
require(magrittr)
require(Mcomp)
require(QDNAseq)
require(matrixStats)
require(NbClust)
require(optparse)
require(PerformanceAnalytics)
require(phangorn)
require(phylogram)
require(plyr)
require(psych)
require(RColorBrewer)
require(reshape2)
require(rlang)
require(smooth)
require(tidyverse) 
require(TTR)
require(vegan)
require(VennDiagram)
require(UpSetR)
options(stringsAsFactors = FALSE)


###################################################################################################
### Create clustering trees
###################################################################################################
fTreeEst <- function(x,data,out='fit') {
  oRet <- NULL
  if (x=='diana') {
    if (out=='fit') oRet <- diana(data)$dc
    if (out=='tree') oRet <- diana(data)
  } else {
    if (out=='fit') oRet <- agnes(data, method = x)$ac
    if (out=='tree') oRet <- agnes(data, method = x)
  }
  return(oRet)
}


###################################################################################################
### Create phylo trees
###################################################################################################
fTreePhylo <- function(data,hamming.dist=FALSE,parsimony=FALSE,add.root=FALSE,return.dist=FALSE) {
  tree <- NULL
  data.dist <- NULL
  if (hamming.dist) {
    data.philo <- data
    data.philo[] <- data[]
    class(data.philo) <- "character"
    data.philo[is.na(data.philo)] <- '?'
    contrast <- matrix(data = c(1,0,0,1,1,1),
                       ncol = 2, byrow = TRUE)
    dimnames(contrast) <- list(c('0','1','?'),c('0','1'))
    if (add.root) {
      data.philo <- cbind(
        data.philo,
        matrix(rep('0',nrow(data.philo)),dimnames=list(rownames(data.philo),c('root')))
      )
    }
    data.philo <- phyDat(
      data = t(data.philo),
      type = 'USER',
      contrast=contrast,
      ambiguity = c("?")
    )
    data.dist <- dist.hamming(data.philo, ratio = TRUE, exclude = "none")
    if (parsimony) {
      tree <- optim.parsimony(NJ(data.dist[]), data.philo)
    } else {
      tree <- NJ(data.dist[])
    }
    if (add.root) {        
      tree <- root(tree,outgroup = "root", resolve.root = TRUE)
    }
  } else {
    try(tree <- NJ(data[]))
  }
  if (return.dist) {
    return(data.dist)
  } else {
    return(tree)
  }
}


###################################################################################################
### Create manual dendrogram
###################################################################################################
fCreateManualTree <- function(l) {
  t <- list()
  t$merge <- matrix(
    c(-1,1:(length(l)-2),
      -1*(2:(length(l)))), 
    nc=2, 
    byrow=FALSE )
  t$height <- rep(0,length(l)-1)
  t$order <- 1:(length(l))
  t$labels <- l
  class(t) <- "hclust" 
  return(t)
}


###################################################################################################
### Created Distance matrix
###################################################################################################
fClustering <- function(a,dist_type='cor.dist') {
  res_dist <- matrix(NA,nrow=ncol(a),ncol=ncol(a))
  rownames(res_dist) <- colnames(a)
  colnames(res_dist) <- colnames(a)
  for (sRow in colnames(a)){
    for (sCol in colnames(a)){
      lRegions <- !is.na(a[,sRow]) & !is.na(a[,sCol])
      if (sum(lRegions) < 3) {
        if (dist_type=='cor.dist') res_dist[sRow,sCol] <- 1.0
        if (dist_type=='cor.p') res_dist[sRow,sCol] <- 1.0
      } else {
        if (sd(a[lRegions,sRow],na.rm=TRUE) == 0 | sd(a[lRegions,sCol],na.rm=TRUE) == 0) {
          if (sd(a[lRegions,sRow],na.rm=TRUE) == sd(a[lRegions,sCol],na.rm=TRUE)) {
            if (dist_type=='cor.dist') res_dist[sRow,sCol] <- 0.0
            if (dist_type=='cor.p') res_dist[sRow,sCol] <- 10^-100
          } else {
            if (dist_type=='cor.dist') res_dist[sRow,sCol] <- 1.0
            if (dist_type=='cor.p') res_dist[sRow,sCol] <- 1.0
          }
        } else {
          if (dist_type=='cor.dist') res_dist[sRow,sCol] <- 1.0 - cor(a[lRegions,sRow],a[lRegions,sCol],use = "complete.obs")
          if (dist_type=='cor.p') res_dist[sRow,sCol] <- max(10^-100,cor.test(a[lRegions,sRow],a[lRegions,sCol],use = "complete.obs")$p.value)
        }
      }
    }
  }
  if (dist_type=='cor.dist') res_dist <- as.dist(res_dist)
  return(res_dist)
}


###################################################################################################
### Add tree distances to HC dist matrix
###################################################################################################
fSetLPdist <- function(data,dist,info,source,color,title) {
  title <- gsub('_fullX_aut','_aut',title,fixed=TRUE)
  title <- gsub('_ARpos_aut','_aut',title,fixed=TRUE)
  title <- gsub('_fullX_chrX','_fullX',title,fixed=TRUE)
  title <- gsub('_ARpos_chrX','_ARpos',title,fixed=TRUE)
  try(if (!title %in% colnames(data)) data[,title] <- NA)
  coladd <- c(C1='gray',C2='gray',C3='gray',C4='gray',C5='gray',C6='gray',C7='gray',C8='gray',C9='gray',C10='gray')
  
  
  ###################################################################################################
  ### Set distance comparison
  ###################################################################################################
  for (sCoreID1 in unique(data[,'Core_ID_1'])) {
    for (sCoreID2 in unique(data[data[,'Core_ID_1']==sCoreID1,'Core_ID_2'])) {
      if (all(c(sCoreID1,sCoreID2) %in% colnames(as.matrix(dist)))) {try({
        data[
          data[,'Core_ID_1']==sCoreID1 &
            data[,'Core_ID_2']==sCoreID2,
          title] <- as.matrix(dist)[sCoreID1,sCoreID2]
      })}
    } 
  }
  
  
  ###################################################################################################
  ### Plot Distance distribution
  ###################################################################################################
  try({
    sTest <- "wilcox.test"
    oDataPlot <- data.frame(
      group = 'ALL',
      type = 'ALL',
      label = 'ALL',
      dist = as.numeric(as.matrix(dist)),
      stringsAsFactors = FALSE
    )
    lComparisons <- list()
    for (sOrgan in names(table(info[colnames(as.matrix(dist)),'Organ']))[
      table(info[colnames(as.matrix(dist)),'Organ'])>2
    ]) {
      oDataPlot <- rbind(oDataPlot,data.frame(
        group = sOrgan,
        type = 'organ',
        label = sOrgan,
        dist = as.numeric(as.matrix(dist)[
          rownames(as.matrix(dist)) %in% info[info[,'Organ']==sOrgan,'Core_ID'],
          colnames(as.matrix(dist)) %in% info[info[,'Organ']==sOrgan,'Core_ID']
        ]),
        stringsAsFactors = FALSE
      ))
      oDataPlot <- rbind(oDataPlot,data.frame(
        group = sOrgan,
        type = 'organ',
        label = paste0(sOrgan,' vs. Other'),
        dist = as.numeric(as.matrix(dist)[
          rownames(as.matrix(dist)) %in% info[info[,'Organ']==sOrgan,'Core_ID'],
          colnames(as.matrix(dist)) %in% info[info[,'Organ']!=sOrgan,'Core_ID']
        ]),
        stringsAsFactors = FALSE
      ))
      lComparisons <- c(lComparisons,list(tmp=c(sOrgan,paste0(sOrgan,' vs. Other'))))
    }
    for (sCluster in names(table(info[colnames(as.matrix(dist)),paste0('clusters_str_',source)]))[
      table(info[colnames(as.matrix(dist)),paste0('clusters_str_',source)])>2
    ]) {
      oDataPlot <- rbind(oDataPlot,data.frame(
        group = paste0('C',sCluster),
        type = 'cluster',
        label = paste0('C',sCluster),
        dist = as.numeric(as.matrix(dist)[
          rownames(as.matrix(dist)) %in% info[info[,paste0('clusters_str_',source)]==sCluster,'Core_ID'],
          colnames(as.matrix(dist)) %in% info[info[,paste0('clusters_str_',source)]==sCluster,'Core_ID']
        ]),
        stringsAsFactors = FALSE
      ))
      oDataPlot <- rbind(oDataPlot,data.frame(
        group = paste0('C',sCluster),
        type = 'cluster',
        label = paste0('C',sCluster,' vs. Other'),
        dist = as.numeric(as.matrix(dist)[
          rownames(as.matrix(dist)) %in% info[info[,paste0('clusters_str_',source)]==sCluster,'Core_ID'],
          colnames(as.matrix(dist)) %in% info[info[,paste0('clusters_str_',source)]!=sCluster,'Core_ID']
        ]),
        stringsAsFactors = FALSE
      ))
      lComparisons <- c(lComparisons,list(tmp=c(paste0('C',sCluster),paste0('C',sCluster,' vs. Other'))))
    }
    oDataPlot <- oDataPlot[oDataPlot[,'dist']>0,]
    oDataTable <- dcast(
      data.frame(table(info[,c('Organ',paste0('clusters_str_',source))])),
      formula = paste0('Organ',intToUtf8(126),'clusters_str_',source) 
    )
    colnames(oDataTable)[colnames(oDataTable)!='Organ'] <- paste0('C',colnames(oDataTable)[colnames(oDataTable)!='Organ'])
    oPlot <- ggplot(oDataPlot,aes(x=label,y=dist,fill=group))
    oPlot <- oPlot + geom_boxplot()
    oPlot <- oPlot + stat_compare_means(
      aes(label = paste0(..method.., 'p =', ..p.format..)),
      label.x = 1,
      method = sTest,
      comparisons = lComparisons)
    oPlot <- oPlot + annotate(
      geom = "table", 
      x = 'contingency', 
      y = 1, 
      label = list(oDataTable), 
      vjust = 1, hjust = 0)
    oPlot <- oPlot + ggtitle(paste(sTest,':',title))
    oPlot <- oPlot + theme_bw()
    oPlot <- oPlot + scale_y_continuous(name=title,limits=c(0,2))
    oPlot <- oPlot + scale_x_discrete(
      name = 'Comparisons across organs / clusters',
      limits=c('ALL',unlist(lComparisons),'contingency','space1','space2'),
      breaks=c('ALL',unlist(lComparisons),'','',''),
      labels=c('ALL',unlist(lComparisons),'','','')
    )
    oPlot <- oPlot + scale_fill_manual(
      name = 'Organ',
      values=c(ALL='black',color,coladd),
      breaks=c('ALL',names(color),names(coladd)),
      labels=c('ALL',names(color),names(coladd))
    )
    plot(oPlot)
  })
  
  
  ###################################################################################################
  ### Return Distance data
  ###################################################################################################
  return(data)
}


###################################################################################################
### Plot tree comparison
###################################################################################################
fCompareTrees <- function(
    tree1,
    tree2,
    treeName1,
    treeName2,
    hamming.data=NULL,
    info = NULL
) {
  marOld <-  par('mar')
  marNew <-  c(8,0,8,0)
  lCTrees <- list()
  lCPlots <- list()
  lCLables <- list()
  nItem <- 0
  for (tree in list(a=tree1,b=tree2)) {
    nItem <- nItem + 1
    if (any(names(tree) == 'order.lab')) {
      lCTrees[[nItem]] <- as.dendrogram(tree)
      lCPlots[[nItem]] <- as.dendrogram(tree)
    }
    if (any(names(tree) == 'tip.label')) {
      lCTrees[[nItem]] <- as.dendrogram(tree)
      lCPlots[[nItem]] <- as.dendrogram(tree)
    }
    lCLables[[nItem]] <- labels(as.dendrogram(lCTrees[[nItem]]))
    lCPlots[[nItem]] <- set(
      lCPlots[[nItem]],
      'labels_colors',
      ifelse(lCLables[[nItem]]=='root','black',
             lTrees[['info']][lCLables[[nItem]],'Organ.Color'])
    ) 
    lCPlots[[nItem]] <- set(
      lCPlots[[nItem]],
      'labels',
      ifelse(lCLables[[nItem]]=='root','root',paste(
        lTrees[['info']][lCLables[[nItem]],'Sample_ID'],
        lTrees[['info']][lCLables[[nItem]],'Organ'])
      )
    )
  }
  par(mfrow=c(2,2),mar=marNew)
  plot(lCPlots[[1]],main=treeName1)
  plot(lCPlots[[2]],main=treeName2)
  plot(lCTrees[[1]],main=treeName1)
  plot(lCTrees[[2]],main=treeName2)
  par(mfrow=c(1,1),mar=marOld)
  if (!is.null(hamming.data) & !is.null(info)) {try({
    draw(Heatmap(
      hamming.data[,lCLables[[1]][lCLables[[1]]!='root']],
      column_title = treeName1,
      name = 'TPs',
      col = c('0'='white','1'='black'),
      row_split = sapply(rownames(hamming.data),function(s) strsplit(s,'_')[[1]][1]),
      cluster_rows = TRUE,
      column_order = lCLables[[1]][lCLables[[1]]!='root'],
      cluster_columns = FALSE,
      show_row_names = FALSE,
      row_title_gp = gpar(fontsize = 8, fontface = "bold"),
      top_annotation = HeatmapAnnotation(
        Organ = lTrees[['info']][lCLables[[1]][lCLables[[1]]!='root'],'Organ'],
        ARcnLAB = oDataStats[lCLables[[1]][lCLables[[1]]!='root'],'ARcnLAB'],
        foo = anno_block(gp = gpar(fill = c('white','gray')),
                         labels_gp = gpar(col = "black", fontsize = 16)),
        col = list(
          Organ = lHMcolOrgan,
          ARcnLAB = lHMcolAR))
    ) + Heatmap(
      hamming.data[,lCLables[[2]][lCLables[[2]]!='root']],
      column_title = treeName2,
      name = 'TPs',
      col = c('0'='white','1'='black'),
      row_split = sapply(rownames(hamming.data),function(s) strsplit(s,'_')[[1]][1]),
      cluster_rows = TRUE,
      column_order = lCLables[[2]][lCLables[[2]]!='root'],
      cluster_columns = FALSE,
      show_row_names = FALSE,
      row_title_gp = gpar(fontsize = 8, fontface = "bold"),
      top_annotation = HeatmapAnnotation(
        Organ = lTrees[['info']][lCLables[[2]][lCLables[[2]]!='root'],'Organ'],
        ARcnLAB = oDataStats[lCLables[[2]][lCLables[[2]]!='root'],'ARcnLAB'],
        foo = anno_block(gp = gpar(fill = c('white','gray')),
                         labels_gp = gpar(col = "black", fontsize = 16)),
        col = list(
          Organ = lHMcolOrgan,
          ARcnLAB = lHMcolAR))
    ))
  })}
}


###################################################################################################
### Set main parameters
###################################################################################################
bCreateSummary <- TRUE
if (file.exists('/Volumes/MahedisData/MH_work_UCL/CASCADE/SCRATCH\'s\ scratch/SCRATCH_diff_binSize/data/TXT')) {
  sPathBED <- '/Volumes/MahedisData/MH_work_UCL/CASCADE/SCRATCH\'s\ scratch/SCRATCH_diff_binSize/data/BED'
  sPathTXT <- '/Volumes/MahedisData/MH_work_UCL/CASCADE/SCRATCH\'s\ scratch/SCRATCH_diff_binSize/data/TXT'
  sPathBAM <- '/Volumes/MahedisData/MH_work_UCL/CASCADE/SCRATCH\'s\ scratch/SCRATCH_diff_binSize/data/BAM'
  sPathRST <- '/Volumes/MahedisData/MH_work_UCL/CASCADE/SCRATCH\'s\ scratch/SCRATCH_diff_binSize/data/BAM'
  sPathOBJ <- 'C:/Data/UCL_stratosphere/data/OBJ/W500000_CA_0.001_all_patients'
  sPathCASCADE <- '/Volumes/MahedisData/MH_work_UCL/CASCADE/SCRATCH\'s\ scratch/SCRATCH_diff_binSize/reports/CASCADE'
} else {
  sPathBED <- '/lustre/scratch/scratch/regmpcr/stratosphere/data/BED'
  sPathTXT <- '/lustre/scratch/scratch/regmpcr/stratosphere/data/TXT'
  sPathBAM <- '/lustre/scratch/scratch/regmpcr/stratosphere/data/BAM'
  sPathRST <- '/lustre/scratch/scratch/regmpcr/stratosphere/data/ROSETREE'
  sPathOBJ <- '/lustre/scratch/scratch/regmpcr/stratosphere/data/OBJ/W500000_CA_0.001_all_patients'
  sPathCASCADE <- '/lustre/scratch/scratch/regmpcr/stratosphere/reports/CASCADE'
}
nPenPloi <- 0.5
nPenCell <- 0.5
nARregionStart <- 66000000
nARregionEnd <- 67500000
lACEchrExlude <- c('Y','MT')
nHMrowgroups <- 5000000
lShifts <- c(0,1,2)
nACEplog10min <- 1
nACEcellularityMin <- 0.20
nACEpgaMin <- 0.05
lAutosomesBinSizes <- c(500000,100000,50000)
lRegionsChrX <- c('fullX','ARpos')
lMeasureTypes <- c('_copynumber','_Segment_Mean','_Copies','_TPs0','_TPs1','_TPs2')
lParsimonyOpts <- c(TRUE)
bRunACEanalyses <- FALSE
bRunSeedRefresh <- FALSE
bForceRT <- TRUE
bRunPatientAnalyses <- TRUE
bRunArchivalAnalyses <- TRUE
bAddExtraDist <- FALSE


###################################################################################################
### Read Input Param
###################################################################################################
option_list = list(
  make_option(c("-a", "--archival"), type="integer", default=1, help="Run Archival Analysis", metavar="integer"),
  make_option(c("-b", "--binsize"), type="integer", default=500000, help="Set bin size", metavar="integer"),
  make_option(c("-s", "--sampleset"), type="character", default="metsplsm,metsonly,allsampl",help="Set of samples", metavar="character"),
  make_option(c("-x", "--chrXregion"), type="character", default="fullX",help="ChrX region", metavar="character"),
  make_option(c("-m", "--measure"), type="character", default="_TPs2,_Copies",help="Measure to use", metavar="character")
)
oOpt <- parse_args(
  OptionParser(option_list=option_list),
  commandArgs(trailingOnly = TRUE)
)
if (oOpt[['archival']] %in% c(0,1)) {
  bRunArchivalAnalyses <- oOpt[['archival']]
}
if (oOpt[['binsize']] >= 1000) {
  nAutosomesBinSize <- oOpt[['binsize']]
  lAutosomesBinSizes <- intersect(lAutosomesBinSizes,nAutosomesBinSize)
}
if (oOpt[['chrXregion']] != 'ALL') {
  sChrXRegion <- oOpt[['chrXregion']]
  lRegionsChrX <- intersect(lRegionsChrX,sChrXRegion)
}
if (oOpt[['measure']] != 'ALL') {
  lMeasureTypes <- intersect(lMeasureTypes,strsplit(oOpt[['measure']],',')[[1]])
  sMeasureType <- strsplit(oOpt[['measure']],',')[[1]][1]
}
sSampleSel <-  oOpt[['sampleset']]


###################################################################################################
### Define region of interest
###################################################################################################
oRegions <- data.frame(
  ID = c('AR','ARenh','centromere'),
  Chromosome = c('X','X','X'),
  Start = c(66763874,66117800,58632012),
  End = c(66950461,66128800,61632012),
  color = c('red','blue','black'),
  stringsAsFactors = FALSE
)
lRegionsCol <- oRegions[,'color']
names(lRegionsCol) <- oRegions[,'ID']
lRegionsCol <- c(c(other='gray'),lRegionsCol)


###################################################################################################
### Define control genes
###################################################################################################
oControlGenes <- data.frame(
  ID = c('NSUN3','NSUN3ext','AP3B1','ZXDB','HCN','AR','ARenh'),
  Chromosome = c('3','3','5','X','11','X','X'),
  Start = c(93781760,90000001,77296349,57618269,65265233,66763874,66117800),
  End = c(93847389,95000000,77590579,57623906,65273940,66950461,66128800),
  stringsAsFactors = FALSE
)
rownames(oControlGenes) <- oControlGenes[,'ID']


###################################################################################################
### Load samples information
###################################################################################################
oDataStats <- read.delim(
  paste0(sPathTXT,'/samples_stats_pga.txt'),
  header=TRUE,
  sep='\t',
  stringsAsFactors=FALSE
)
oDataStats <- oDataStats[
  !is.na(oDataStats[,'Trial.ID']) &
    oDataStats[,'Trial.ID'] == "CA63" &
    !is.na(oDataStats[,'Run_ID']) &
    !is.na(oDataStats[,'Core_ID']) &
    !is.na(oDataStats[,'Organ']) &
    !is.na(oDataStats[,'ACE_cellularity_all']) &
    !is.na(oDataStats[,'ACE_ploidy_all']),
]
oDataStats[,'ARcnLABlog2'] <- log2(oDataStats[,'ARcnLAB'])
oDataStats <- oDataStats[
  oDataStats[,'Trial.ID'] %in% names(table(unique(oDataStats[!is.na(oDataStats[,'Run_ID']),c('Trial.ID','Sample_ID')])[,'Trial.ID']))[
    table(unique(oDataStats[!is.na(oDataStats[,'Run_ID']),c('Trial.ID','Sample_ID')])[,'Trial.ID'])>1],
]


###################################################################################################
### Define sample sets for the analyses
###################################################################################################
lSampleSets <- list(
  allsampl = oDataStats[,'Core_ID'],
  metsonly = oDataStats[
    regexpr('archival',oDataStats[,'Tissue.Site']) == -1 &
      !oDataStats[,'Organ'] %in% c('Plasma'),
    'Core_ID'],
  metsplsm = oDataStats[
    regexpr('archival',oDataStats[,'Tissue.Site']) == -1,
    'Core_ID'],
  metsnoln = oDataStats[
    regexpr('archival',oDataStats[,'Tissue.Site']) == -1 &
      !oDataStats[,'Organ'] %in% c('Plasma','Lymph nodes'),
    'Core_ID']
)
if (oOpt[['sampleset']] != 'ALL') {
  if (oOpt[['sampleset']] == 'MAIN') {
    lSampleSets <- lSampleSets[c('metsonly','allsampl')]
  } else {
    lSampleSets <- lSampleSets[strsplit(sSampleSel,',')[[1]]]
  }
  sSampleSel <- names(lSampleSets)[1]
}


###################################################################################################
### List sample BAMs
###################################################################################################
lBAMs <- setdiff(
  list.files(sPathBAM,'.RMDUP.SORTED.bam'),
  list.files(sPathBAM,'.RMDUP.SORTED.bam.bai')
)
names(lBAMs) <- gsub('.RMDUP.SORTED.bam','',lBAMs,fixed=TRUE)
names(lBAMs) <- sapply(names(lBAMs),function(s) strsplit(s,'_WGSsf_')[[1]][1])
lSamplesBAMs <- lBAMs[names(lBAMs) %in% oDataStats[,'Core_ID']]
print('Duplicated samples by run:')
print(table(oDataStats[,c('Sample_ID','Run_ID')])[rowSums(table(oDataStats[,c('Sample_ID','Run_ID')]))>1,])
print('Unmatched BAMs:')
setdiff(names(lBAMs),names(lSamplesBAMs))
print('Unmatched SamplesIDs:')
setdiff(unique(na.omit(oDataStats[,'Core_ID'])),names(lSamplesBAMs))


##################################################################################################
### Load RoseTree TC
###################################################################################################
oDataRSTtc <- read.delim(
  file.path(sPathRST,'patient_rostree_estimate_tc.txt'),
  header=TRUE,
  sep='\t',
  stringsAsFactors=FALSE
)


###################################################################################################
### Load SV files from HC
###################################################################################################
if (!file.exists(file.path(sPathTXT,'SV_archive.txt'))) {
  oResSV <- NULL
  for (sSVfile in list.files(file.path(sPathTXT,'../SV'))) {
    sSampleID <- gsub('.txt','',gsub('all_SV_','',sSVfile))
    oDataSV <- read.delim(
      file.path(sPathTXT,'../SV',sSVfile),
      header = FALSE,
      sep = '\t',
      stringsAsFactors = FALSE
    )
    colnames(oDataSV) <- c('chr_start','pos_start','type','chr_end','pos_end')
    if (paste0(sSampleID,'hc') %in% oDataStats[,'Core_ID']) {
      oDataSV[,'Core_ID'] <- paste0(sSampleID,'hc')
    } else {
      oDataSV[,'Core_ID'] <- sSampleID
    }
    oDataSV <- merge(
      oDataStats[,c('Trial.ID','Sample_ID','Core_ID')],
      oDataSV,
      all.y=TRUE
    )
    oResSV <- rbind(oResSV,oDataSV)
  }
  oResSV <- oResSV[complete.cases(oResSV),]
  oResSV <- unique(oResSV)
  write.table(
    oResSV,
    file = file.path(sPathTXT,'SV_archive.txt'),
    col.names = TRUE,
    row.names = FALSE,
    sep= '\t',
    quote= FALSE
  )
} else {
  oResSV <- read.delim(
    file.path(sPathTXT,'SV_archive.txt'),
    header = TRUE,
    sep = '\t',
    stringsAsFactors = FALSE
  )
}


##################################################################################################
### Load distances from HC
###################################################################################################
oDataHCdist <- read.delim(
  file.path(sPathTXT,'sample_distance_by_clonal_mutation.txt'),
  header=TRUE,
  sep='\t',
  stringsAsFactors=FALSE
)


###################################################################################################
### Loop through differnt bin sizes to created sample data
###################################################################################################
dir.create(file.path(sPathCASCADE,'details_bypatient'),recursive = TRUE)
if (bRunACEanalyses) {
  nAutosomesBinSize <- 500000
  for (nAutosomesBinSize in lAutosomesBinSizes) {
    if (!file.exists(paste0(sPathCASCADE,'/CN__',as.integer(nAutosomesBinSize),'__byPatient.txt.gz'))) {
      ###################################################################################################
      ### Load bin info
      ###################################################################################################
      oBins <- getBinAnnotations(as.integer(nAutosomesBinSize/1000), genome="hg19", type="SR50")
      
      
      ##################################################################################################
      ### Load RoseTree AF data
      ###################################################################################################
      oDataRSTcn <- read.delim(
        file.path(sPathRST,'patient_rostree_estimate_cn.txt'),
        header=TRUE,
        sep='\t',
        stringsAsFactors=FALSE
      )
      oDataRSTcn[,'chr_ori'] <- oDataRSTcn[,'chr']
      oDataRSTcn[,'start_ori'] <- oDataRSTcn[,'start']
      oDataRSTcn[,'end_ori'] <- oDataRSTcn[,'end']
      oDataRSTcn[,'chr'] <- gsub('chr','',oDataRSTcn[,'chr_ori'])
      oDataRSTcn[,'start'] <- as.integer(ceiling(oDataRSTcn[,'start_ori'] / nAutosomesBinSize) * nAutosomesBinSize)
      oDataRSTcn[,'end'] <- as.integer(floor(oDataRSTcn[,'end_ori'] / nAutosomesBinSize) * nAutosomesBinSize)
      oDataRSTcn <- merge(oDataRSTtc,oDataRSTcn,by='Trial.ID')
      
      
      ###################################################################################################
      ### Process all samples
      ###################################################################################################
      sTrialID <- '5017'
      for(sTrialID in unique(oDataStats[,'Trial.ID'])) {try({
        ###################################################################################################
        ### Set outpit files
        ###################################################################################################
        sFileOut <- paste0(
          sPathCASCADE,'/details_bypatient/',
          sTrialID,'__',
          as.integer(nAutosomesBinSize))
        
        
        ###################################################################################################
        ### Run Sample analyses
        ###################################################################################################
        if (!file.exists(paste0(sFileOut,'.txt')) & (
          file.exists(paste0(sFileOut,'.RData')) |
          length(intersect(
            names(lSamplesBAMs),
            oDataStats[oDataStats[,'Trial.ID']==sTrialID,'Core_ID']
          )) > 0 
        )) {
          write.table(NULL,file=paste0(sFileOut,'.txt'))
          
          
          ###################################################################################################
          ### Run QDNAseq
          ###################################################################################################
          if (file.exists(paste0(sFileOut,'.RData'))) {
            load(paste0(sFileOut,'.RData'))
          } else {
            lCores <- intersect(
              names(lSamplesBAMs),
              oDataStats[oDataStats[,'Trial.ID']==sTrialID,'Core_ID']
            )
            oQM <- list()
            oQM[['ReadCounts']] <- binReadCounts(
              oBins,
              bamfiles=file.path(sPathBAM,lSamplesBAMs[lCores]),
              chunkSize=1000000,
              cache=TRUE,
              bamnames=lCores)
            
          }
          if (!any(names(oQM)=='seeds') | bRunSeedRefresh) {
            oQM[['seeds']] <- .GlobalEnv$.Random.seed
          } else {
            .GlobalEnv$.Random.seed <- oQM[['seeds']]
          }
          if (!any(names(oQM)=='SegmentsNorm') | bRunSeedRefresh) {
            oQM[['ReadFiltered']] <- applyFilters(oQM[['ReadCounts']], residual=TRUE, blacklist=TRUE)
            oQM[['ReadCorrected']] <- estimateCorrection(oQM[['ReadFiltered']])
            oQM[['Bins']]  <- correctBins(
              applyFilters(oQM[['ReadCorrected']] ,
                           residual=TRUE,
                           blacklist=TRUE,
                           chromosomes=lACEchrExlude),
              method="ratio",
              adjustIncompletes=TRUE,
              chromosomes=lACEchrExlude)
            oQM[['BinsNorm']] <- normalizeBins(oQM[['Bins']])
            oQM[['BinsSmooth']] <- smoothOutlierBins(oQM[['BinsNorm']])
            #~~~~~~~~~~~~~~~~~~~~~~``
            oQM[['Segments']] <- segmentBins(oQM[['BinsSmooth']],transformFun="log2")
            oQM[['SegmentsNorm']] <- normalizeSegmentedBins(oQM[['Segments']])
            oQM[['SegmentsCalls']]  <- callBins(
              oQM[['SegmentsNorm']] ,
              method='CGHcall',
              organism='human',
              prior = 'all',
              cellularity=1)
            oQM[['CGH']]  <- makeCgh(oQM[['SegmentsCalls']] ,filter=FALSE)
          }
          
          
          ###################################################################################################
          ### Run ACE
          ###################################################################################################
          oSM <- list(data = oQM[['SegmentsNorm']])
          sSampleID <- sampleNames(oQM[['SegmentsNorm']])[1]
          for (sSampleID in sampleNames(oQM[['SegmentsNorm']])) {
            nSampleID <- which(sampleNames(oQM[['SegmentsNorm']])==sSampleID)
            oSM[['model']][[sSampleID]] <- squaremodel(
              oQM[['SegmentsNorm']],
              QDNAseqobjectsample = nSampleID,
              prows=100,
              ptop=5,
              pbottom=1,
              method = 'RMSE',
              exclude = lACEchrExlude,
              penalty = nPenCell,
              penploidy = nPenPloi,
              highlightminima = TRUE
            )
            oSM[['minimadf']][[sSampleID]] <- oSM[['model']][[sSampleID]][['minimadf']]
            oRTtc <-  oDataRSTtc[
              oDataRSTtc[,'Sample_ID'] == oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'] |
                paste0(oDataRSTtc[,'Sample_ID'],'hc') == oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'] |
                paste0(oDataRSTtc[,'Sample_ID'],'hcfull') == oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'],
            ]
            if (nrow(oRTtc)>0) {
              if (any(
                oSM[['minimadf']][[sSampleID]][,'ploidy'] >= oRTtc[,'RT_ploidy_min'] &
                oSM[['minimadf']][[sSampleID]][,'ploidy'] <= oRTtc[,'RT_ploidy_max'] &
                oSM[['minimadf']][[sSampleID]][,'cellularity'] >= oRTtc[,'RT_cellularity_min'] &
                oSM[['minimadf']][[sSampleID]][,'cellularity'] <= oRTtc[,'RT_cellularity_max']
              )) {
                oSM[['minimadf']][[sSampleID]] <- oSM[['minimadf']][[sSampleID]][
                  oSM[['minimadf']][[sSampleID]][,'ploidy'] >= oRTtc[,'RT_ploidy_min'] &
                    oSM[['minimadf']][[sSampleID]][,'ploidy'] <= oRTtc[,'RT_ploidy_max'] &
                    oSM[['minimadf']][[sSampleID]][,'cellularity'] >= oRTtc[,'RT_cellularity_min'] &
                    oSM[['minimadf']][[sSampleID]][,'cellularity'] <= oRTtc[,'RT_cellularity_max'],
                ]
              } else {
                if (bForceRT) {
                  oSM[['minimadf']][[sSampleID]] <- rbind(
                    data.frame(
                      ploidy = (oRTtc[,'RT_ploidy_min'] + oRTtc[,'RT_ploidy_max']) / 2,
                      cellularity = (oRTtc[,'RT_cellularity_min'] + oRTtc[,'RT_cellularity_max']) / 2,
                      error = -1,
                      minimum = TRUE
                    ),
                    oSM[['minimadf']][[sSampleID]]
                  )
                }
              }
            }
            oSM[['bestparam']][[sSampleID]] <- oSM[['minimadf']][[sSampleID]][oSM[['minimadf']][[sSampleID]][,'error']==min(oSM[['minimadf']][[sSampleID]][,'error']),]
            oSM[['bestparam']][[sSampleID]][,'source'] <- 'ACE'
            if (any(oSM[['bestparam']][[sSampleID]][,'error']==-1)) oSM[['bestparam']][[sSampleID]][oSM[['bestparam']][[sSampleID]][,'error']==-1,'source'] <- 'RT'
            if (nrow(oSM[['bestparam']][[sSampleID]])>1) {
              oSM[['bestparam']][[sSampleID]] <- oSM[['bestparam']][[sSampleID]][order(
                oSM[['bestparam']][[sSampleID]][,'cellularity'],
                oSM[['bestparam']][[sSampleID]][,'ploidy'],
                method='radix',
                decreasing=c(TRUE,FALSE)
              ),]
            }
            oSM[['CN']][[sSampleID]] <- getadjustedsegments(
              objectsampletotemplate(oSM[['data']],nSampleID),
              cellularity = oSM[['bestparam']][[sSampleID]][1,'cellularity'],
              ploidy = oSM[['bestparam']][[sSampleID]][1,'ploidy'],
              standard = 1,
              log = FALSE
            )
          }
          
          
          ###################################################################################################
          ### Select QDNAseq data
          ###################################################################################################
          oRes <- NULL
          for (sType in c('copynumber','segmented','calls')) {
            if (sType == 'copynumber') oTmp <- data.frame(copynumber(oQM[['CGH']]))
            if (sType == 'segmented') oTmp <- data.frame(segmented(oQM[['CGH']]))
            if (sType == 'calls') oTmp <- data.frame(calls(oQM[['CGH']]))
            colnames(oTmp) <- paste(colnames(oTmp),sType,sep='_')
            oTmp[,'region'] <- as.character(rownames(oTmp))
            if (is.null(oRes)) {
              oRes <- oTmp
            } else {
              oRes <- merge(oRes,oTmp,by='region',all=TRUE)
            }
          }
          oRes[,'Chromosome'] <- sapply(oRes[,'region'],function(s) strsplit(s,':')[[1]][1])
          oRes[,'Start'] <- sapply(oRes[,'region'],function(s) as.integer(strsplit(s,'-')[[1]][2])-nAutosomesBinSize+1)
          oRes[,'End'] <- as.integer(oRes[,'Start'] + nAutosomesBinSize - 1)
          
          
          ###################################################################################################
          ### Select ACE data
          ###################################################################################################
          for (sSampleID in sampleNames(oQM[['SegmentsNorm']])) {try({
            print(sSampleID)
            oTmp <- oSM[['CN']][[sSampleID]]
            oTmp[!is.na(oTmp[,'P_log10']) & oTmp[,'P_log10'] == -Inf,'P_log10'] <- -100
            colnames(oTmp)[colnames(oTmp)=='Chromosome'] <- 'segment_chr'
            colnames(oTmp)[colnames(oTmp)=='Start'] <- 'segment_start'
            colnames(oTmp)[colnames(oTmp)=='End'] <- 'segment_end'
            oTmp[,'segment_start'] <- as.integer(oTmp[,'segment_start'])
            oTmp[,'segment_end'] <- as.integer(oTmp[,'segment_end'])
            oTmp[,'segment'] <- paste0(
              oTmp[,'segment_chr'],':',
              formatC(oTmp[,'segment_start']),'-',
              formatC(oTmp[,'segment_end'])
            )
            oTmp <- oTmp[order(
              as.integer(gsub('Y','24', gsub('X','23',oTmp[,'segment_chr']))),
              oTmp[,'segment_start']),]
            oTmp[,'cn_type_start'] <- 0
            oTmp[,'cn_type_end'] <- 0
            for (nRow in 1:nrow(oTmp)) {
              if (nRow>1) {
                if (oTmp[nRow,'segment_chr']==oTmp[nRow-1,'segment_chr']) {
                  oTmp[nRow,'cn_type_start'] <- sign(oTmp[nRow,'Copies'] - oTmp[nRow-1,'Copies'])
                }
              }
              if (nRow<nrow(oTmp)) {
                if (oTmp[nRow,'segment_chr']==oTmp[nRow+1,'segment_chr']) {
                  oTmp[nRow,'cn_type_end'] <- sign(oTmp[nRow,'Copies'] - oTmp[nRow+1,'Copies'])
                }
              }
            }
            colnames(oTmp)[!colnames(oTmp) %in% c('segment_chr','segment_start','segment_end')] <- paste(
              sSampleID,
              colnames(oTmp)[!colnames(oTmp) %in% c('segment_chr','segment_start','segment_end')],
              sep = '_'
            )
            oTmp[,'queryHits'] <- 1:nrow(oTmp)
            
            
            ###################################################################################################
            ### Fix ACE segments calls
            ###################################################################################################
            if (any(oTmp[,'segment_start']>oTmp[,'segment_end'])) {
              for (sSegToFix in oTmp[oTmp[,'segment_start']>oTmp[,'segment_end'],'segment']) {
                print(paste('Fixing segment:',sSegToFix))
                if (any(oTmp[,'segment_chr']==oTmp[oTmp[,'segment']==sSegToFix,'segment_chr'] &
                        oTmp[,'segment_start'] > oTmp[oTmp[,'segment']==sSegToFix,'segment_start'])) {
                  nNewEnd <- min(oTmp[
                    oTmp[,'segment_chr']==oTmp[oTmp[,'segment']==sSegToFix,'segment_chr'] &
                      oTmp[,'segment_start'] > oTmp[oTmp[,'segment']==sSegToFix,'segment_start'],
                    'segment_start'])-1
                  print(paste('Set new End to next segment:',nNewEnd))
                } else {
                  nNewEnd <- max(oTmp[
                    !is.na(oTmp[,'calls']) &
                      oTmp[,'segment_chr']==oTmp[oTmp[,'segment']==sSegToFix,'segment_chr'],
                    'segment_end'])
                  print(paste('Set new End chromosome stop:',nNewEnd))
                }
                oTmp[oTmp[,'segment']==sSegToFix,'segment_end'] <- nNewEnd
              }
              oTmp[,'segment'] <- paste0(
                oTmp[,'segment_chr'],':',
                formatC(oTmp[,'segment_start']),'-',
                formatC(oTmp[,'segment_end'])
              )
            }
            
            
            ###################################################################################################
            ### Merge QDNAseq and ACE results
            ###################################################################################################
            oTmp[,'queryHits'] <- 1:nrow(oTmp)
            oRes[,'subjectHits'] <- 1:nrow(oRes)
            oRes <- merge(
              oRes,
              as.data.frame(findOverlaps(
                makeGRangesFromDataFrame(
                  oTmp,
                  ignore.strand = TRUE,
                  keep.extra.columns = TRUE,
                  seqnames.field = 'segment_chr',
                  start.field = 'segment_start',
                  end.field = 'segment_end'
                ),
                makeGRangesFromDataFrame(
                  oRes,
                  ignore.strand = TRUE,
                  keep.extra.columns = TRUE,
                  seqnames.field = 'Chromosome',
                  start.field = 'Start',
                  end.field = 'End'
                ),
                ignore.strand = TRUE
              )),
              by = 'subjectHits',
              all = TRUE
            )
            oRes <- merge(oRes,oTmp,by = 'queryHits',all.x=TRUE)
            oRes <- oRes[,colnames(oRes)!='queryHits']
            oRes <- oRes[,colnames(oRes)!='subjectHits']
            if (any(!is.na(oRes[,'segment_chr']) & !is.na(oRes[,'Chromosome']) & oRes[,'segment_chr'] != oRes[,'Chromosome'])) print(sSampleID)
            colnames(oRes)[colnames(oRes) %in% c('segment','segment_chr','segment_start','segment_end')] <- paste(
              sSampleID,
              colnames(oRes)[colnames(oRes) %in% c('segment','segment_chr','segment_start','segment_end')],
              sep = '_'
            )
            
          })}
          oRes <- oRes[order(
            as.integer(gsub('Y','24', gsub('X','23',oRes[,'Chromosome']))),
            oRes[,'Start']),]
          
          
          ###################################################################################################
          ### Export data
          ###################################################################################################
          save(oQM,oSM,oRes,file=paste0(sFileOut,'.RData'))
          write.table(
            oRes,
            file=paste0(sFileOut,'.txt'),
            sep='\t',
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE
          )
        }
        
        
        ###################################################################################################
        ### Create sample report
        ###################################################################################################
        if (!file.exists(paste0(sFileOut,'.models.pdf'))) {try({
          if (file.exists(paste0(sFileOut,'.RData'))) {
            ###################################################################################################
            ### Load data objects
            ###################################################################################################
            load(paste0(sFileOut,'.RData'))
            
            
            ###################################################################################################
            ### Create plots for each sample
            ###################################################################################################
            pdf(paste0(sFileOut,'.models.pdf'),paper='special',width=24,height=12)
            sSampleID <- sampleNames(oSM[['data']])[1]
            for (sSampleID in sampleNames(oSM[['data']])) {
              nSampleID <- (1:length(sampleNames(oQM[['SegmentsNorm']])))[
                sampleNames(oQM[['SegmentsNorm']])==sSampleID]
              
              
              ###################################################################################################
              ### Get RT data
              ###################################################################################################
              oRTtc <-  oDataRSTtc[
                oDataRSTtc[,'Sample_ID'] == oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'] |
                  paste0(oDataRSTtc[,'Sample_ID'],'hc') == oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'] |
                  paste0(oDataRSTtc[,'Sample_ID'],'hcfull') == oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'],
              ]
              
              
              ###################################################################################################
              ### Get Sample CN profile from ACE
              ###################################################################################################
              oPlotCN <- oRes[,c('Chromosome','Start','End',paste0(sSampleID,c('_Segment_Mean','_Copies')))]
              oPlotCN <- oPlotCN[order(
                as.integer(gsub('Y','24',gsub('X','23',oPlotCN[,'Chromosome']))),
                oPlotCN[,'Start']),]
              oPlotCN[,'bin'] <- 1:nrow(oPlotCN)
              oPlotCN[sapply(c(1:22,'X','Y'),function(s) min(oPlotCN[oPlotCN[,'Chromosome']==s,'bin'])),'chrstart'] <- c(1:22,'X','Y')
              oPlotCN[sapply(c(1:22,'X','Y'),function(s) floor(mean(oPlotCN[oPlotCN[,'Chromosome']==s,'bin']))),'chrmid'] <- c(1:22,'X','Y')
              
              
              ###################################################################################################
              ### Start grid plot 
              ###################################################################################################
              grid.arrange(
                ###################################################################################################
                ### Plot CN without t.c. adj.
                ###################################################################################################
                singleplot(
                  objectsampletotemplate(oSM[['data']],nSampleID),
                  cellularity = 1,
                  ploidy = 2,
                  onlyautosomes = FALSE,
                  title = paste(c(sSampleID,'- no adjustment'))
                ),
                
                
                ###################################################################################################
                ### Plot CN with t.c. adj.
                ###################################################################################################
                singleplot(
                  objectsampletotemplate(oSM[['data']],nSampleID),
                  cellularity = oSM[['bestparam']][[sSampleID]][1,'cellularity'],
                  ploidy = oSM[['bestparam']][[sSampleID]][1,'ploidy'],
                  error = round(oSM[['bestparam']][[sSampleID]][1,'error']),
                  onlyautosomes = FALSE,
                  standard = 1,
                  cap = max(oPlotCN[,paste0(sSampleID,'_Copies')],na.rm=TRUE)+2,
                  title = paste(c(
                    sSampleID,
                    'cellularity',oSM[['bestparam']][[sSampleID]][1,'cellularity'],
                    'ploidy',oSM[['bestparam']][[sSampleID]][1,'ploidy'],
                    'error',round(oSM[['bestparam']][[sSampleID]][1,'error'])
                  ))),
                
                
                ###################################################################################################
                ### Plot ACE optimal solutions
                ###################################################################################################
                oSM[['model']][[sSampleID]]$matrixplot +
                  geom_point(
                    data = oSM[['model']][[sSampleID]][['minimadf']][
                      oSM[['model']][[sSampleID]][['minimadf']][,'error']==min(
                        oSM[['model']][[sSampleID]][['minimadf']][,'error']),
                    ],
                    aes(x=cellularity,y=ploidy),
                    shape=4,
                    size=4,
                    stroke=2,
                    color='black'
                  ) +
                  geom_point(
                    data = oSM[['bestparam']][[sSampleID]],
                    aes(x=cellularity,y=ploidy,shape=source),
                    size=4,
                    stroke=2,
                    color='darkgreen'
                  ) +
                  geom_rect(
                    data = oRTtc,
                    aes(
                      xmin=RT_cellularity_min,
                      xmax=RT_cellularity_max,
                      ymin=RT_ploidy_min,
                      ymax=RT_ploidy_max
                    ),
                    alpha = 0.5
                  ) +
                  scale_shape_manual(values=c('ACE'=4,'RT'=13)) +
                  theme_classic() +
                  ggtitle(paste(
                    sSampleID,'(',
                    oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'],'):',
                    'cellularity',oSM[['bestparam']][[sSampleID]][1,'cellularity'],
                    'ploidy',oSM[['bestparam']][[sSampleID]][1,'ploidy'],
                    'error',round(oSM[['bestparam']][[sSampleID]][1,'error'],2)
                  )),
                
                
                ###################################################################################################
                ### Plot ACE optimal solutions
                ###################################################################################################
                ggplot2::ggplot() +
                  geom_vline(
                    data=oPlotCN[!is.na(oPlotCN[,'chrstart']),],
                    aes(xintercept = bin), 
                    color = "#666666", 
                    linetype = "dashed") +
                  geom_hline(
                    yintercept = seq(0,max(oPlotCN[,paste0(sSampleID,'_Copies')],na.rm=TRUE)+2), 
                    color = 'lightgray', 
                    size = 0.5) +
                  geom_point(
                    data=oPlotCN,
                    aes(x=bin,y=get(paste0(sSampleID,'_Segment_Mean'))),
                    size = 0.5, 
                    color = 'gray', 
                    shape = 24
                  ) +
                  geom_point(
                    data=oPlotCN,
                    aes(x=bin,y=get(paste0(sSampleID,'_Copies'))),
                    size = 1, 
                    color = 'darkorange', 
                    shape = 24
                  ) +
                  scale_y_continuous(
                    name = "copies", 
                    limits = c(0,max(oPlotCN[,paste0(sSampleID,'_Copies')],na.rm=TRUE)+2), 
                    breaks = seq(0,max(oPlotCN[,paste0(sSampleID,'_Copies')],na.rm=TRUE)+2), 
                    expand=c(0,0)) +
                  scale_x_continuous(
                    name = "chromosome", 
                    limits = c(0,max(oPlotCN[,'bin'],na.rm=TRUE)), 
                    breaks = oPlotCN[!is.na(oPlotCN[,'chrmid']),'bin'], labels = oPlotCN[!is.na(oPlotCN[,'chrmid']),'chrmid'], expand = c(0,0)) +
                  theme_classic() + theme(
                    axis.line = element_line(color='black'), 
                    axis.ticks = element_line(color='black'), 
                    axis.text = element_text(color='black')) +
                  ggtitle(paste(
                    sSampleID,'(',
                    oDataStats[oDataStats[,'Core_ID']==sSampleID,'Sample_ID'],'):',
                    'cellularity',oSM[['bestparam']][[sSampleID]][1,'cellularity'],
                    'ploidy',oSM[['bestparam']][[sSampleID]][1,'ploidy'],
                    'error',round(oSM[['bestparam']][[sSampleID]][1,'error'],2)
                  )) +
                  theme(plot.title = element_text(hjust = 0.5)),
                
                
                ###################################################################################################
                ### Define grid arrange
                ###################################################################################################
                ncol=2,nrow=2,
                layout_matrix=matrix(nrow=2,ncol=2,c(1,2,3,4))
              )
              
              
              ###################################################################################################
              ### Remove objects
              ###################################################################################################
              rm(oPlotCN,oRTtc)
            }
            try(dev.off());try(dev.off());try(dev.off());try(dev.off())
          }
        })}
      })}
    }
  }
}


###################################################################################################
### Filter samples to analyse
###################################################################################################
oDataStats[oDataStats[,'Run_ID'] %in% c('HiSeqII','G20','G51','G84','H06','H30','H49','I06','I23'),'IN_PLOT'] <- 1
oDataStats[is.na(oDataStats[,'IN_PLOT']),'IN_PLOT'] <- 0
oDataStats <- oDataStats[oDataStats[,'IN_PLOT']==1,]
oDataStats <- oDataStats[!oDataStats[,'Core_ID'] %in% c('H670009','H670010','H670011','H670012'),]
rownames(oDataStats) <- oDataStats[,'Core_ID']


###################################################################################################
### Loop through differnt bin sizes to created sample data
###################################################################################################
if (bRunPatientAnalyses) {
  lInfoSamples <- list() 
  lDataPlotChrX <- list()
  lDataPlotChrXregionCol <- list()
  lSamplesOrderCluster <- list()
  
  
  ###################################################################################################
  ### Loop through differnt bin sizes and sample selections
  ###################################################################################################
  for (sSampleSel in names(lSampleSets)) {
    for (nAutosomesBinSize in c(lAutosomesBinSizes)) {
      dir.create(file.path(
        sPathCASCADE,
        paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize))
      ))
      nChrXBinSize <- nAutosomesBinSize
      
      
      ###################################################################################################
      ### Loop through differnt patients
      ###################################################################################################
      for (sTrialID in unique(oDataStats[,'Trial.ID'])) {try({
        if (file.exists(file.path(
          sPathCASCADE,
          'details_bypatient',
          paste0(sTrialID,'__',as.integer(nAutosomesBinSize),'.RData')))) {
          
          
          ###################################################################################################
          ### Load patient CN analysis data
          ###################################################################################################
          print(paste(Sys.time(),'Load patient CN analysis data for',sTrialID))
          load(file.path(
            sPathCASCADE,
            'details_bypatient',
            paste0(sTrialID,'__',as.integer(nAutosomesBinSize),'.RData')))
          rownames(oRes) <- oRes[,'region']
          
          
          ###################################################################################################
          ### Get samples metadata
          ###################################################################################################
          print(paste(Sys.time(),'Get samples metadata for',sTrialID))
          lTrees <- list()
          lTrees[['info_all']] <- unique(oDataStats[
            oDataStats[,'Trial.ID'] == sTrialID,
            c('Trial.ID','Core_ID','Sample_ID','Run_ID','Figure_N',
              'Tissue.Site','Organ','Organ.Color','ARcnLAB','ARcnLABlog2',
              'QCAmin37_MEAN_COVERAGE'
            )])
          rownames(lTrees[['info_all']]) <- lTrees[['info_all']][,'Core_ID']
          
          
          ###################################################################################################
          ### Add ACE information for the samples
          ###################################################################################################
          lTrees[['info_all']][,'ACE_rosetree_cellularity'] <- sapply(
            lTrees[['info_all']][,'Core_ID'],
            function(s) 
              oSM[['bestparam']][[s]][1,'cellularity']
          )
          lTrees[['info_all']][,'ACE_rosetree_ploidy'] <- sapply(
            lTrees[['info_all']][,'Core_ID'],
            function(s) 
              oSM[['bestparam']][[s]][1,'ploidy']
          )
          lTrees[['info_all']][,'ACE_rosetree_error'] <- sapply(
            lTrees[['info_all']][,'Core_ID'],
            function(s) 
              oSM[['bestparam']][[s]][1,'error']
          )
          lTrees[['info_all']][,'ACE_pga_autosomes'] <- sapply(
            lTrees[['info_all']][,'Core_ID'],
            function(s) (
              sum(na.omit(oRes[oRes[,'Chromosome']!='X',paste0(s,'_Copies')])!=2) / 
                length(na.omit(oRes[oRes[,'Chromosome']!='X',paste0(s,'_Copies')]))
            )
          )
          
          
          ###################################################################################################
          ### Add Rosetree information for the samples
          ###################################################################################################
          lTrees[['info_all']] <- merge(
            lTrees[['info_all']],
            oDataRSTtc,
            all.x=TRUE
          )
          lTrees[['info_all']][,'sample_removed'] <- ''
          rownames(lTrees[['info_all']]) <- lTrees[['info_all']][,'Core_ID']
          
          
          ###################################################################################################
          ### Filter samples with low t.c.
          ###################################################################################################
          print(paste(Sys.time(),'Remove samples with low t.c. for',sTrialID))
          print('List of removed samples:')
          print(lTrees[['info_all']][
            lTrees[['info_all']][,'ACE_rosetree_cellularity'] < nACEcellularityMin |
              (
                lTrees[['info_all']][,'ACE_rosetree_cellularity'] >= 0.95 &
                  lTrees[['info_all']][,'ACE_pga_autosomes'] < nACEpgaMin
              ),
            'Core_ID'
          ])
          lTrees[['info_all']][
            lTrees[['info_all']][,'ACE_rosetree_cellularity'] < nACEcellularityMin |
              (
                lTrees[['info_all']][,'ACE_rosetree_cellularity'] >= 0.95 &
                  lTrees[['info_all']][,'ACE_pga_autosomes'] < nACEpgaMin
              ),
            'sample_removed'
          ] <- 'due to low t.c.'
          
          
          ###################################################################################################
          ### Remove samples to skip
          ###################################################################################################
          lTrees[['info_arch']] <- lTrees[['info_all']][
            (
              lTrees[['info_all']][,'Organ'] %in% c('Plasma','Prostate') |
                regexpr('archival',lTrees[['info_all']][,'Tissue.Site'])>-1 |
                regexpr('prostate',lTrees[['info_all']][,'Tissue.Site'])>-1 |
                regexpr('Prostate',lTrees[['info_all']][,'Tissue.Site'])>-1
            ) &
              lTrees[['info_all']][,'Run_ID'] == 'CASCADE' &
              lTrees[['info_all']][,'sample_removed']=='',
          ]
          
          
          ###################################################################################################
          ### Filter samples by sample set
          ###################################################################################################
          lTrees[['info_all']][!lTrees[['info_all']][,'Core_ID'] %in% lSampleSets[[sSampleSel]],'sample_removed'] <- paste(
            lTrees[['info_all']][!lTrees[['info_all']][,'Core_ID'] %in% lSampleSets[[sSampleSel]],'sample_removed'],
            'due to organ/site selection',
            sep=';'
          )
          
          
          ###################################################################################################
          ### Export full samples list
          ###################################################################################################
          if (all(lTrees[['info_all']][,'Run_ID']=='CASCADE')) {
            lInfoSamples[[paste0(sSampleSel,'__A',as.integer(nAutosomesBinSize))]] <- rbind(
              lInfoSamples[[paste0(sSampleSel,'__A',as.integer(nAutosomesBinSize))]],
              lTrees[['info_all']]
            )
          } 
          write.table(
            lTrees[['info_all']],
            file=file.path(
              sPathCASCADE,
              paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
              paste0('Models__',sTrialID,'__',sSampleSel,'__A',as.integer(nAutosomesBinSize),'__unfiltered_sample_list.txt')),
            sep='\t',
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE
          )
          
          
          ###################################################################################################
          ### Remove samples to skip
          ###################################################################################################
          lTrees[['info_all']] <- lTrees[['info_all']][lTrees[['info_all']][,'sample_removed']=='',]
          
          
          ###################################################################################################
          ### Set organ colors
          ###################################################################################################
          print(paste(Sys.time(),'Set organ colors for',sTrialID))
          lOrganCols <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ.Color']
          names(lOrganCols) <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ']
          
          
          ###################################################################################################
          ### Run correlation analysis for Autosomes
          ###################################################################################################
          if (nrow(lTrees[['info_all']])>2) {
            ###################################################################################################
            ### Select regions corresponding to Autosomal breakpoints
            ###################################################################################################
            print(paste(Sys.time(),'Select regions corresponding to Autosomal breakpoints for',sTrialID))
            lTrees[['regions_aut']] <- data.frame(
              region=as.character(),
              side=as.character(),
              Chromosome=as.character(),
              Start=as.integer(),
              End=as.integer()
            )
            lTrees[['TPs0']] <- data.frame(
              tp_side=as.character(),
              tp_type=as.integer(),
              tp_chr=as.character(),
              tp_pos_min=as.integer(),
              tp_pos_max=as.integer()
            )
            for (sCurrentSample in lTrees[['info_all']][,'Core_ID']) {try({
              oTmp <- rbind(
                unique(oRes[
                  !is.na(oRes[,paste0(sCurrentSample,'_segment')]) &
                    !is.na(oRes[,paste0(sCurrentSample,'_segmented')]) &
                    oRes[,paste0(sCurrentSample,'_segment_chr')] == oRes[,'Chromosome'] &
                    oRes[,paste0(sCurrentSample,'_segment_start')] == oRes[,'Start'] & 
                    oRes[,paste0(sCurrentSample,'_cn_type_start')] != 0,
                  c('Chromosome','Start','End',
                    paste0(sCurrentSample,'_cn_type_start'))]) %>% mutate(
                      tp_id = paste(
                        'S',
                        ifelse(get(paste0(sCurrentSample,'_cn_type_start'))>0,'U','D'),
                        Chromosome,
                        as.integer(Start),
                        sep='_'
                      ),
                      tp_side = 'S',
                      tp_type = get(paste0(sCurrentSample,'_cn_type_start')),
                      tp_chr = Chromosome,
                      tp_pos_min = Start,
                      tp_pos_max = Start
                    ) %>% dplyr::select(
                      tp_id, tp_side, tp_type, tp_chr, tp_pos_min, tp_pos_max
                    ),
                unique(oRes[
                  !is.na(oRes[,paste0(sCurrentSample,'_segment')]) &
                    !is.na(oRes[,paste0(sCurrentSample,'_segmented')]) &
                    oRes[,paste0(sCurrentSample,'_segment_chr')] == oRes[,'Chromosome'] &
                    oRes[,paste0(sCurrentSample,'_segment_end')] == oRes[,'End'] & 
                    oRes[,paste0(sCurrentSample,'_cn_type_end')] != 0,
                  c('Chromosome','Start','End',
                    paste0(sCurrentSample,'_cn_type_end'))]) %>% mutate(
                      tp_id = paste(
                        'E',
                        ifelse(get(paste0(sCurrentSample,'_cn_type_end'))>0,'U','D'),
                        Chromosome,
                        as.integer(End),
                        sep='_'
                      ),
                      tp_side = 'E',
                      tp_type = get(paste0(sCurrentSample,'_cn_type_end')),
                      tp_chr = Chromosome,
                      tp_pos_min = End,
                      tp_pos_max = End
                    ) %>% dplyr::select(
                      tp_id, tp_side, tp_type, tp_chr,  tp_pos_min, tp_pos_max
                    )
              )
              colnames(oTmp)[colnames(oTmp)=='tp_id'] <- sCurrentSample
              lTrees[['TPs0']] <- merge(
                lTrees[['TPs0']],
                oTmp,
                all = TRUE
              )
              lTrees[['regions_aut']] <- merge(lTrees[['regions_aut']], rbind(
                cbind(
                  unique(oRes[
                    !is.na(oRes[,paste0(sCurrentSample,'_segment')]) &
                      !is.na(oRes[,paste0(sCurrentSample,'_segmented')]) &
                      oRes[,paste0(sCurrentSample,'_segment_chr')] == oRes[,'Chromosome'] &
                      oRes[,paste0(sCurrentSample,'_segment_start')] == oRes[,'Start'] &
                      oRes[,paste0(sCurrentSample,'_cn_type_start')] != 0,
                    c('region','Chromosome','Start','End',paste0(sCurrentSample,'_segment')),
                    drop = FALSE
                  ]),
                  data.frame(side='S', stringsAsFactors=FALSE)
                ),
                cbind(
                  unique(oRes[
                    !is.na(oRes[,paste0(sCurrentSample,'_segment')]) &
                      !is.na(oRes[,paste0(sCurrentSample,'_segmented')]) &
                      oRes[,paste0(sCurrentSample,'_segment_chr')] == oRes[,'Chromosome'] &
                      oRes[,paste0(sCurrentSample,'_segment_end')] == oRes[,'End'] &
                      oRes[,paste0(sCurrentSample,'_cn_type_end')] != 0,
                    c('region','Chromosome','Start','End',paste0(sCurrentSample,'_segment')),
                    drop = FALSE
                  ]),
                  data.frame(side='E', stringsAsFactors=FALSE)
                )),
                all.x = TRUE,
                all.y = TRUE
              )
            })}
            colnames(lTrees[['regions_aut']]) <- gsub('_segment','',colnames(lTrees[['regions_aut']]),fixed=TRUE)
            lTrees[['regions_chrX']] <- lTrees[['regions_aut']][lTrees[['regions_aut']][,'Chromosome']=='X',]
            lTrees[['regions_aut']] <- lTrees[['regions_aut']][lTrees[['regions_aut']][,'Chromosome'] %in% c(1:22),]
            lTrees[['regions_aut']] <- lTrees[['regions_aut']][
              order(as.integer(lTrees[['regions_aut']][,'Chromosome']),
                    as.integer(lTrees[['regions_aut']][,'Start'])),
            ]
            lTrees[['regions_chrX']] <- lTrees[['regions_chrX']][
              order(as.integer(lTrees[['regions_chrX']][,'Start'])),
            ]
            oRes <- oRes[order(
              as.integer(gsub('Y','24', gsub('X','23',oRes[,'Chromosome']))),
              oRes[,'Start']),]
            if (nrow(lTrees[['regions_aut']]) > 0) rownames(lTrees[['regions_aut']]) <- paste0(lTrees[['regions_aut']][,'region'],'_',lTrees[['regions_aut']][,'side'])
            if (nrow(lTrees[['regions_chrX']]) > 0) rownames(lTrees[['regions_chrX']]) <- paste0(lTrees[['regions_chrX']][,'region'],'_',lTrees[['regions_chrX']][,'side'])
            
            
            ###################################################################################################
            ### Identify TPs at different shift
            ###################################################################################################
            lTrees[['TPs0']] <- lTrees[['TPs0']][
              order(
                as.integer(gsub('Y','24', gsub('X','23',lTrees[['TPs0']][,'tp_chr']))),
                as.integer(lTrees[['TPs0']][,'tp_pos_min'])),
            ]
            rownames(lTrees[['TPs0']]) <- paste(
              lTrees[['TPs0']][,'tp_chr'],
              as.integer(lTrees[['TPs0']][,'tp_pos_min']),
              as.integer(lTrees[['TPs0']][,'tp_pos_max']),
              lTrees[['TPs0']][,'tp_side'],
              ifelse(lTrees[['TPs0']][,'tp_type']>0,'U','D'),
              sep='_'
            )
            for (nShift in lShifts[lShifts > 0]) {
              sTPtable <- paste0('TPs',nShift)
              lTrees[[sTPtable]] <- lTrees[['TPs0']]
              for (nRow in 1:nrow(lTrees[[sTPtable]])) {
                lTrees[[sTPtable]][
                  lTrees[[sTPtable]][,'tp_side'] == lTrees[[sTPtable]][nRow,'tp_side'] &
                    lTrees[[sTPtable]][,'tp_type'] == lTrees[[sTPtable]][nRow,'tp_type'] &
                    lTrees[[sTPtable]][,'tp_chr'] == lTrees[[sTPtable]][nRow,'tp_chr'] &
                    lTrees[[sTPtable]][,'tp_pos_min'] >= lTrees[[sTPtable]][nRow,'tp_pos_min'] &
                    lTrees[[sTPtable]][,'tp_pos_max'] <= lTrees[[sTPtable]][nRow,'tp_pos_max'] + (nShift * nAutosomesBinSize),
                  'tp_pos_min'
                ] <- min(lTrees[[sTPtable]][
                  lTrees[[sTPtable]][,'tp_side'] == lTrees[[sTPtable]][nRow,'tp_side'] &
                    lTrees[[sTPtable]][,'tp_type'] == lTrees[[sTPtable]][nRow,'tp_type'] &
                    lTrees[[sTPtable]][,'tp_chr'] == lTrees[[sTPtable]][nRow,'tp_chr'] &
                    lTrees[[sTPtable]][,'tp_pos_min'] >= lTrees[[sTPtable]][nRow,'tp_pos_min'] &
                    lTrees[[sTPtable]][,'tp_pos_max'] <= lTrees[[sTPtable]][nRow,'tp_pos_max'] + (nShift * nAutosomesBinSize),
                  'tp_pos_min'
                ])
                lTrees[[sTPtable]][
                  lTrees[[sTPtable]][,'tp_side'] == lTrees[[sTPtable]][nRow,'tp_side'] &
                    lTrees[[sTPtable]][,'tp_type'] == lTrees[[sTPtable]][nRow,'tp_type'] &
                    lTrees[[sTPtable]][,'tp_chr'] == lTrees[[sTPtable]][nRow,'tp_chr'] &
                    lTrees[[sTPtable]][,'tp_pos_min'] >= lTrees[[sTPtable]][nRow,'tp_pos_min'] &
                    lTrees[[sTPtable]][,'tp_pos_max'] <= lTrees[[sTPtable]][nRow,'tp_pos_max'] + (nShift * nAutosomesBinSize),
                  'tp_pos_max'
                ] <- max(lTrees[[sTPtable]][
                  lTrees[[sTPtable]][,'tp_side'] == lTrees[[sTPtable]][nRow,'tp_side'] &
                    lTrees[[sTPtable]][,'tp_type'] == lTrees[[sTPtable]][nRow,'tp_type'] &
                    lTrees[[sTPtable]][,'tp_chr'] == lTrees[[sTPtable]][nRow,'tp_chr'] &
                    lTrees[[sTPtable]][,'tp_pos_min'] >= lTrees[[sTPtable]][nRow,'tp_pos_min'] &
                    lTrees[[sTPtable]][,'tp_pos_max'] <= lTrees[[sTPtable]][nRow,'tp_pos_max'] + (nShift * nAutosomesBinSize),
                  'tp_pos_max'
                ])
              }
              lTrees[[sTPtable]] <- ddply(
                lTrees[[sTPtable]],
                c('tp_side','tp_type','tp_chr','tp_pos_min','tp_pos_max'),
                function(o) {
                  oTmp <- o[,!colnames(o) %in% c('tp_side','tp_type','tp_chr','tp_pos_min','tp_pos_max')]
                  for (sCol in colnames(oTmp)) {
                    oTmp[,sCol] <- paste(unique(na.omit(oTmp[,sCol])),sep=',',collapse=',')
                  }
                  return(unique(oTmp))
                }
              )
              lTrees[[sTPtable]] <- lTrees[[sTPtable]][
                order(
                  as.integer(gsub('Y','24', gsub('X','23',lTrees[[sTPtable]][,'tp_chr']))),
                  as.integer(lTrees[[sTPtable]][,'tp_pos_min'])),
              ]
              lTrees[[sTPtable]][,'shift_warning'] <- ifelse(
                lTrees[[sTPtable]][,'tp_pos_max'] - lTrees[[sTPtable]][,'tp_pos_min'] > nShift * nAutosomesBinSize,
                1,
                0
              )
              lTrees[[sTPtable]][,'shift_error'] <- ifelse(
                lTrees[[sTPtable]][,'tp_pos_max'] - lTrees[[sTPtable]][,'tp_pos_min'] > 2 * nShift * nAutosomesBinSize,
                lTrees[[sTPtable]][,'tp_pos_max'] - lTrees[[sTPtable]][,'tp_pos_min'],
                0
              )
              rownames(lTrees[[sTPtable]]) <- paste(
                lTrees[[sTPtable]][,'tp_chr'],
                as.integer(lTrees[[sTPtable]][,'tp_pos_min']),
                as.integer(lTrees[[sTPtable]][,'tp_pos_max']),
                lTrees[[sTPtable]][,'tp_side'],
                ifelse(lTrees[[sTPtable]][,'tp_type']>0,'U','D'),
                sep='_'
              )
            }
            
            
            ###################################################################################################
            ### Export TPs data
            ###################################################################################################
            for (nShift in lShifts) {
              sTPtable <- paste0('TPs',nShift)
              lTrees[[sTPtable]][is.na(lTrees[[sTPtable]])] <- ''
              write.table(
                lTrees[[sTPtable]],
                file=file.path(
                  sPathCASCADE,
                  paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                  paste0('Models__',sTrialID,'__',sSampleSel,'__A',as.integer(nAutosomesBinSize),'__',sTPtable,'.txt')),
                sep='\t',
                quote=FALSE,
                col.names=TRUE,
                row.names=TRUE
              )
            }
            
            
            ###################################################################################################
            ### Created data form hamming dist. based on transition points
            ###################################################################################################
            for (nShift in lShifts) {
              lTrees[[paste0('TPs',nShift,'_aut')]] <- lTrees[[paste0('TPs',nShift)]][,!colnames(lTrees[[paste0('TPs',nShift)]]) %in% c(
                'tp_side','tp_type','tp_chr','tp_pos_min','tp_pos_max','shift_warning','shift_error')]
              lTrees[[paste0('TPs',nShift,'_aut')]][lTrees[[paste0('TPs',nShift,'_aut')]] != ''] <- 1
              lTrees[[paste0('TPs',nShift,'_aut')]][lTrees[[paste0('TPs',nShift,'_aut')]] == ''] <- 0
              lTrees[[paste0('TPs',nShift,'_aut')]] <- as.matrix(lTrees[[paste0('TPs',nShift,'_aut')]])
              class(lTrees[[paste0('TPs',nShift,'_aut')]]) <- "numeric"
              lTrees[[paste0('TPs',nShift,'_chrX')]] <- lTrees[[paste0('TPs',nShift,'_aut')]][regexpr('^X_',rownames(lTrees[[paste0('TPs',nShift,'_aut')]]))>-1,]
              lTrees[[paste0('TPs',nShift,'_aut')]] <- lTrees[[paste0('TPs',nShift,'_aut')]][regexpr('^X_',rownames(lTrees[[paste0('TPs',nShift,'_aut')]]))==-1,]
            }
            
            
            ###################################################################################################
            ### If Autosomal breakpoints are detected
            ###################################################################################################
            if (length(lTrees[['regions_aut']])>0) {try({
              lTrees[['data_aut_all']] <- oRes[oRes[,'Chromosome'] %in% c(1:22),]
              lTrees[['data_chrX_all']] <- oRes[oRes[,'Chromosome']=='X',]
              
              
              ###################################################################################################
              ### Add annotation for the region of interest
              ###################################################################################################
              print(paste(Sys.time(),'Add annotation for the region of interest for',sTrialID))
              lTrees[['data_chrX_all']][,'annotation'] <- 'other'
              for (nRegion in 1:nrow(oRegions)) {
                lTrees[['data_chrX_all']][
                  lTrees[['data_chrX_all']][,'Chromosome'] == oRegions[nRegion,'Chromosome'] & (
                    (lTrees[['data_chrX_all']][,'Start']<=oRegions[nRegion,'Start'] & lTrees[['data_chrX_all']][,'End']>=oRegions[nRegion,'End']) |
                      (lTrees[['data_chrX_all']][,'Start']>=oRegions[nRegion,'Start'] & lTrees[['data_chrX_all']][,'End']<=oRegions[nRegion,'End']) |
                      (lTrees[['data_chrX_all']][,'Start']>=oRegions[nRegion,'Start'] & lTrees[['data_chrX_all']][,'Start']<=oRegions[nRegion,'End']) |
                      (lTrees[['data_chrX_all']][,'End']>=oRegions[nRegion,'Start'] & lTrees[['data_chrX_all']][,'End']<=oRegions[nRegion,'End'])
                  ),
                  'annotation'] <- oRegions[nRegion,'ID']
              }
              rownames(lTrees[['data_chrX_all']]) <- lTrees[['data_chrX_all']][,'region']
              
              
              ###################################################################################################
              ### Add annotation for SVs
              ###################################################################################################
              print(paste(Sys.time(),'Add annotation for SVs for',sTrialID))
              for (sSampleID in unique(oResSV[oResSV[,'Trial.ID']==sTrialID,'Sample_ID'])) {
                print(paste(Sys.time(),'Start SV:',sSampleID))
                lTrees[['data_chrX_all']] <- lTrees[['data_chrX_all']][,colnames(lTrees[['data_chrX_all']]) != paste0('SV_',sSampleID)]
                oResSVtmpS <- oResSV[
                  complete.cases(oResSV) &
                    oResSV[,'Sample_ID']==sSampleID,
                  c('chr_start','pos_start','type')]
                oResSVtmpS[,'region'] <- paste0(
                  oResSVtmpS[,'chr_start'],':',
                  as.integer(floor(oResSVtmpS[,'pos_start']/nChrXBinSize)*nChrXBinSize+1),'-',
                  as.integer(ceiling(oResSVtmpS[,'pos_start']/nChrXBinSize)*nChrXBinSize)
                )
                oResSVtmpE <- oResSV[
                  complete.cases(oResSV) &
                    oResSV[,'Sample_ID']==sSampleID,
                  c('chr_end','pos_end','type')]
                oResSVtmpE[,'region'] <- paste0(
                  oResSVtmpE[,'chr_end'],':',
                  as.integer(floor(oResSVtmpE[,'pos_end']/nChrXBinSize)*nChrXBinSize+1),'-',
                  as.integer(ceiling(oResSVtmpE[,'pos_end']/nChrXBinSize)*nChrXBinSize)
                )
                oResSVtmp <- rbind(
                  oResSVtmpS[,c('region','type')],
                  oResSVtmpE[,c('region','type')]
                )
                oResSVtmp <- ddply(oResSVtmp,.(region),function(o) 
                  data.frame(tmp=paste0(sort(unique(o[,'type'])),collapse='_'))
                )
                rownames(oResSVtmp) <- oResSVtmp[,'region']
                colnames(oResSVtmp) <- c('region',paste0('SV_',sSampleID))
                lTrees[['data_chrX_all']] <- merge(lTrees[['data_chrX_all']],oResSVtmp,all.x=TRUE)
                print(paste(Sys.time(),'End'))
              }
              if (length(grep(paste0('SV_',sTrialID),colnames(lTrees[['data_chrX_all']]),value=TRUE))==1) {
                lTrees[['data_chrX_all']][,paste0('SV_ALL_',sTrialID)] <- ifelse(
                  !is.na(lTrees[['data_chrX_all']][,grep(paste0('SV_',sTrialID),colnames(lTrees[['data_chrX_all']]),value=TRUE)]),
                  'SVs',
                  NA)
              }
              if (length(grep(paste0('SV_',sTrialID),colnames(lTrees[['data_chrX_all']]),value=TRUE))>1) {
                lTrees[['data_chrX_all']][,paste0('SV_ALL_',sTrialID)] <- ifelse(
                  rowSums(!is.na(lTrees[['data_chrX_all']][,grep(paste0('SV_',sTrialID),colnames(lTrees[['data_chrX_all']]),value=TRUE)]))>0,
                  'SVs',
                  NA)
              }
              rownames(lTrees[['data_chrX_all']]) <- lTrees[['data_chrX_all']][,'region']
              gc();gc();gc();gc()
              
              
              ###################################################################################################
              ### Define SV colors 
              ###################################################################################################
              print(paste(Sys.time(),'Define SV colors for',sTrialID))
              lSVtypes <- unique(c('SVs',unlist(lTrees[['data_chrX_all']][,grep('SV_',colnames(lTrees[['data_chrX_all']]),value=TRUE)])))
              names(lSVtypes) <- lSVtypes
              lSVtypes[!is.na(lSVtypes) & regexpr('DEL',lSVtypes)> -1 & regexpr('DUP',lSVtypes)> -1] <- 'green'
              lSVtypes[!is.na(lSVtypes) & regexpr('DEL',lSVtypes)> -1] <- 'blue'
              lSVtypes[!is.na(lSVtypes) & regexpr('DUP',lSVtypes)> -1] <- 'red'
              lSVtypes[!is.na(lSVtypes) & !lSVtypes %in% c('green','red','blue')] <- 'black'
              
              
              ###################################################################################################
              ### Sort data by chr loc
              ###################################################################################################
              lTrees[['data_aut_all']] <- lTrees[['data_aut_all']][order(
                as.integer(lTrees[['data_aut_all']][,'Chromosome']),
                as.integer(lTrees[['data_aut_all']][,'Start'])
              ),]
              lTrees[['data_chrX_all']] <- lTrees[['data_chrX_all']][order(
                as.integer(lTrees[['data_chrX_all']][,'Start'])
              ),]
              
              
              ###################################################################################################
              ### Run the analysis based on differnt metrics / regions
              ###################################################################################################
              print(paste(Sys.time(),'Run the analysis based on differnt metrics / regions for',sTrialID))
              for (sMeasureType in lMeasureTypes) {
                for (sChrXRegion in lRegionsChrX) {try({
                  lTrees[['info']] <- lTrees[['info_all']]
                  
                  
                  ###################################################################################################
                  ### Select the region to analyze
                  ###################################################################################################
                  print(paste(Sys.time(),'Run the analysis based on differnt metrics / regions for',sTrialID,sMeasureType,sChrXRegion))
                  if (sChrXRegion=='fullX') {
                    lChrXRegions <- lTrees[['data_chrX_all']][lTrees[['data_chrX_all']][,'Chromosome']=='X','region']
                  }
                  if (sChrXRegion=='ARpos') {
                    lChrXRegions <- lTrees[['data_chrX_all']][
                      lTrees[['data_chrX_all']][,'Chromosome']=='X' & (
                        (lTrees[['data_chrX_all']][,'Start']<=nARregionStart & lTrees[['data_chrX_all']][,'End']>=nARregionEnd) |
                          (lTrees[['data_chrX_all']][,'Start']>=nARregionStart & lTrees[['data_chrX_all']][,'End']<=nARregionEnd) |
                          (lTrees[['data_chrX_all']][,'Start']>=nARregionStart & lTrees[['data_chrX_all']][,'Start']<=nARregionEnd) |
                          (lTrees[['data_chrX_all']][,'End']>nARregionStart & lTrees[['data_chrX_all']][,'End']<=nARregionEnd)
                      ),
                      'region']
                  }
                  lTrees[['data_chrX_all']] <- lTrees[['data_chrX_all']]
                  oRegionCol <- data.frame(
                    Region=lTrees[['data_chrX_all']][lChrXRegions,'annotation'],
                    stringsAsFactors=FALSE
                  )
                  rownames(oRegionCol) <- lChrXRegions
                  lRegionCol <- list(Region=lRegionsCol)
                  
                  
                  ###################################################################################################
                  ### Select the measure for the clustering
                  ###################################################################################################
                  print(paste(Sys.time(),'Select the measure for the clustering for',sTrialID,sMeasureType,sChrXRegion))
                  if (regexpr('TPs',sMeasureType)>-1) {
                    lTrees[['data_aut']] <- lTrees[[paste0(substr(sMeasureType,2,10),'_aut')]]
                    lTrees[['data_chrX']] <- lTrees[[paste0(substr(sMeasureType,2,10),'_chrX')]]
                  } else {
                    lTrees[['data_aut']] <- lTrees[['data_aut_all']][unique(lTrees[['regions_aut']][,'region']),]
                    lTrees[['data_aut']] <- lTrees[['data_aut']][,grep(paste0(sMeasureType,'$'),colnames(lTrees[['data_aut']]),value=TRUE)]
                    colnames(lTrees[['data_aut']]) <- gsub(sMeasureType,'',colnames(lTrees[['data_aut']]))
                    lTrees[['data_chrX']] <- lTrees[['data_chrX_all']][lChrXRegions,]
                    lTrees[['data_chrX']] <- lTrees[['data_chrX']][,grep(paste0(sMeasureType,'$'),colnames(lTrees[['data_chrX']]),value=TRUE)]
                    colnames(lTrees[['data_chrX']]) <- gsub(sMeasureType,'',colnames(lTrees[['data_chrX']]))
                  }
                  
                  
                  ###################################################################################################
                  ### Remove samples missing one data type
                  ###################################################################################################
                  print(paste(Sys.time(),'Remove samples missing one data type for',sTrialID,sMeasureType,sChrXRegion))
                  print("List of samples missing one data type")
                  print(lTrees[['info']][
                    !lTrees[['info']][,'Core_ID'] %in% colnames(lTrees[['data_aut']]) |
                      !lTrees[['info']][,'Core_ID'] %in% colnames(lTrees[['data_aut']]),
                    'Core_ID'
                  ])
                  lTrees[['info']] <- lTrees[['info']][intersect(
                    lTrees[['info']][,'Core_ID'],
                    colnames(lTrees[['data_chrX']])
                  ),]
                  lTrees[['info']] <- lTrees[['info']][intersect(
                    lTrees[['info']][,'Core_ID'],
                    colnames(lTrees[['data_aut']])
                  ),]
                  lTrees[['data_aut']] <- lTrees[['data_aut']][,intersect(
                    lTrees[['info']][,'Core_ID'],
                    colnames(lTrees[['data_aut']])
                  )]
                  lTrees[['data_chrX']] <- lTrees[['data_chrX']][,intersect(
                    lTrees[['info']][,'Core_ID'],
                    colnames(lTrees[['data_chrX']])
                  )]
                  print(dim(lTrees[['info']]))
                  print(dim(lTrees[['data_aut']]))
                  print(dim(lTrees[['data_chrX']]))
                  
                  
                  ###################################################################################################
                  ### Aggreate chrX regions for the heatmap
                  ###################################################################################################
                  print(paste(Sys.time(),'Aggreate chrX regions for the heatmap for',sTrialID,sMeasureType,sChrXRegion))
                  lChrXRegionsGroupKB <- paste0(
                    lTrees[['data_chrX_all']][lChrXRegions,'Chromosome'],':',
                    formatC(as.integer(((floor(lTrees[['data_chrX_all']][lChrXRegions,'Start']/nHMrowgroups)*nHMrowgroups))+1),width=10,flag='0'),'-',
                    formatC(as.integer(((ceiling(lTrees[['data_chrX_all']][lChrXRegions,'End']/nHMrowgroups)+1)*nHMrowgroups)),width=10,flag='0')
                  )
                  names(lChrXRegionsGroupKB) <- lChrXRegions
                  lDataPlotAnno <- lTrees[['data_chrX_all']][lChrXRegions,'annotation']
                  names(lDataPlotAnno) <- lChrXRegions
                  
                  
                  ###################################################################################################
                  ### Set heatmap annotations
                  ###################################################################################################
                  print(paste(Sys.time(),'Set heatmap annotations for',sTrialID,sMeasureType,sChrXRegion))
                  lHMcolAR <- colorRamp2(c(0,2,3,5,10), c("green","green", "yellow", "orange","red"))
                  lHMcolCor <- colorRamp2(c(0,0.25,0.5,1,2), c("red","orange","white","white",'blue'))
                  lHMcolCorP <- colorRamp2(c(0,5,20,100), c("white","yellow","orange","red"))
                  lHMcolOrgan <- unique(lTrees[['info']][,c('Organ','Organ.Color')])[,'Organ.Color']
                  names(lHMcolOrgan) <- unique(lTrees[['info']][,c('Organ','Organ.Color')])[,'Organ']
                  if (sMeasureType == '_Copies') {
                    lHMcol <- lHMcolAR
                  } 
                  if (sMeasureType == '_Segment_Mean') {
                    lHMcol <- colorRamp2(c(
                      min(lTrees[['data_chrX']],na.rm=TRUE),
                      quantile(lTrees[['data_chrX']],0.25,na.rm=TRUE),
                      quantile(lTrees[['data_chrX']],0.50,na.rm=TRUE),
                      quantile(lTrees[['data_chrX']],0.75,na.rm=TRUE), 
                      min(quantile(lTrees[['data_chrX']],0.99,na.rm=TRUE),10)), 
                      c("#ffffff","#ffffff", "#fff1cb", "#e52722","#e52722"))
                  }
                  if (sMeasureType == '_copynumber') {
                    lHMcol <- colorRamp2(c(
                      min(c(-2,min(lTrees[['data_chrX']],na.rm=TRUE),na.rm=TRUE)),
                      -1,
                      0,
                      quantile(lTrees[['data_chrX']],0.75,na.rm=TRUE), 
                      quantile(lTrees[['data_chrX']],0.99,na.rm=TRUE)), 
                      c("#ffffff","#ffffff", "#fff1cb", "#e52722","#e52722"))
                  }
                  
                  
                  ###################################################################################################
                  ### Remove centromeres and not significant segments in chrX
                  ###################################################################################################
                  print(paste(Sys.time(),'Remove centromeres and not significant segments in chrX for',sTrialID,sMeasureType,sChrXRegion))
                  if (regexpr('TPs',sMeasureType)==-1) {
                    lTrees[['data_chrX']][
                      rownames(lTrees[['data_chrX']]) %in% lTrees[['data_chrX_all']][
                        lTrees[['data_chrX_all']][,'annotation'] == 'centromere','region'],
                    ] <- NA
                  }
                  
                  
                  ###################################################################################################
                  ### Remove NA samples
                  ###################################################################################################
                  print(paste(Sys.time(),'Remove NA samples for',sTrialID,sMeasureType,sChrXRegion))
                  lTrees[['data_chrX']] <- lTrees[['data_chrX']][,c(apply(
                    lTrees[['data_chrX']],
                    2,
                    function(l) !all(is.na(l))
                  ))]
                  lTrees[['data_aut']] <- lTrees[['data_aut']][,
                                                               colnames(lTrees[['data_chrX']])
                  ]
                  lTrees[['info']] <- lTrees[['info']][
                    lTrees[['info']][,'Core_ID'] %in% colnames(lTrees[['data_chrX']]),
                  ]
                  print(dim(lTrees[['info']]))
                  print(dim(lTrees[['data_aut']]))
                  print(dim(lTrees[['data_chrX']]))
                  
                  
                  ###################################################################################################
                  ### Estimate distances and related p val.
                  ###################################################################################################
                  print(paste(Sys.time(),'Estimate distances and related p val. for',sTrialID,sMeasureType,sChrXRegion))
                  if (regexpr('TPs',sMeasureType)>-1) {
                    lTrees[['dist_aut']] <- lTrees[['data_aut']]
                    lTrees[['dist_aut']][] <- lTrees[['data_aut']][]
                    class(lTrees[['dist_aut']]) <- "character"
                    lTrees[['dist_aut']] <- phyDat(
                      data = t(lTrees[['dist_aut']]),
                      type = 'USER',
                      contrast=matrix(data = c(1,0,0,1,1,1),ncol = 2, byrow = TRUE, dimnames = list(c('0','1','?'),c('0','1'))),
                      ambiguity = c("?")
                    )
                    lTrees[['dist_aut']] <- dist.hamming(lTrees[['dist_aut']], ratio = TRUE, exclude = "none")
                    lTrees[['cor.p_aut']] <- as.matrix(lTrees[['dist_aut']])
                    lTrees[['cor.p_aut']][] <- 0
                    lTrees[['dist_chrX']] <- lTrees[['data_chrX']]
                    lTrees[['dist_chrX']][] <- lTrees[['data_chrX']][]
                    class(lTrees[['dist_chrX']]) <- "character"
                    lTrees[['dist_chrX']] <- phyDat(
                      data = t(lTrees[['dist_chrX']]),
                      type = 'USER',
                      contrast=matrix(data = c(1,0,0,1,1,1),ncol = 2, byrow = TRUE, dimnames = list(c('0','1','?'),c('0','1'))),
                      ambiguity = c("?")
                    )
                    lTrees[['dist_chrX']] <- dist.hamming(lTrees[['dist_chrX']], ratio = TRUE, exclude = "none")
                    lTrees[['cor.p_chrX']] <- as.matrix(lTrees[['dist_chrX']])
                    lTrees[['cor.p_chrX']][] <- 0
                  } else {
                    lTrees[['dist_aut']] <- fClustering(
                      lTrees[['data_aut']]
                    )
                    lTrees[['cor.p_aut']] <- fClustering(
                      lTrees[['data_aut']],
                      dist_type='cor.p'
                    )
                    lTrees[['dist_chrX']] <- fClustering(
                      lTrees[['data_chrX']]
                    )
                    lTrees[['cor.p_chrX']] <- fClustering(
                      lTrees[['data_chrX']],
                      dist_type='cor.p'
                    )
                  }
                  
                  
                  ###################################################################################################
                  ### Create the PDF for clustering results
                  ###################################################################################################
                  print(paste(Sys.time(),'Create the PDF for clustering results for',sTrialID,sMeasureType,sChrXRegion))
                  pdf(paste0(
                    sPathCASCADE,
                    paste0('/S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                    '/Models__',sTrialID,'__',sSampleSel,'_',sMeasureType,'__',sChrXRegion,'__A',as.integer(nAutosomesBinSize),'__X',as.integer(nChrXBinSize),'.pdf'),
                    paper='special',width=16,height=8)
                  
                  
                  ###################################################################################################
                  ### Run cluster analysis based on the different data sources
                  ###################################################################################################
                  for (sChrSource in c('aut','chrX')) {try({
                    ###################################################################################################
                    ### Define best clustering algorithm
                    ###################################################################################################
                    print(paste(Sys.time(),'Define best clustering algorithm for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    lTrees[[paste0('trees_comp_',sChrSource)]] <- sapply(
                      c( average='average',single='single',complete='complete',ward='ward',diana='diana'),
                      fTreeEst,
                      data = lTrees[[paste0('dist_',sChrSource)]]
                    )
                    
                    
                    ###################################################################################################
                    ### Define best number of clusters
                    ###################################################################################################
                    print(paste(Sys.time(),'Define best number of clusters for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    lTrees[[paste0('clust_comp_',sChrSource)]] <- fviz_nbclust(
                      x = as.matrix(lTrees[[paste0('dist_',sChrSource)]]), 
                      FUNcluster  = hcut, 
                      k.max = min(ncol(lTrees[[paste0('data_',sChrSource)]])/2,10),
                      method = "silhouette"
                    )
                    
                    
                    ###################################################################################################
                    ### Create manual tree if all distances are 0s
                    ###################################################################################################
                    print(paste(Sys.time(),'Create manual tree if all distances are 0s for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    if (all(
                      rowSums(as.matrix(lTrees[[paste0('dist_',sChrSource)]])==0)==ncol(
                        lTrees[[paste0('data_',sChrSource)]])
                    )) {
                      ###################################################################################################
                      ### Crete manual tree with 0s dist
                      ###################################################################################################
                      lTrees[[paste0('tree_',sChrSource)]] <- fCreateManualTree(
                        colnames(lTrees[[paste0('data_',sChrSource)]])
                      )
                      
                      
                      ###################################################################################################
                      ### Set all samples in the same cluster
                      ###################################################################################################
                      lTrees[[paste0('clusters_',sChrSource)]] <- cutree(
                        lTrees[[paste0('tree_',sChrSource)]], 
                        k = 1
                      )
                      lTrees[[paste0('clusters_str_',sChrSource)]] <- as.character(cutree(
                        lTrees[[paste0('tree_',sChrSource)]], 
                        k = 1
                      ))
                      names(lTrees[[paste0('clusters_str_',sChrSource)]]) <- names(cutree(
                        lTrees[[paste0('tree_',sChrSource)]], 
                        k = 1
                      ))
                    } else {
                      ###################################################################################################
                      ### Crete the tree
                      ###################################################################################################
                      lTrees[[paste0('tree_',sChrSource)]] <- fTreeEst(names(lTrees[[paste0('trees_comp_',sChrSource)]])[
                        !is.na(lTrees[[paste0('trees_comp_',sChrSource)]]) &
                          lTrees[[paste0('trees_comp_',sChrSource)]]==max(
                            lTrees[[paste0('trees_comp_',sChrSource)]],
                            na.rm=TRUE
                          )][1],
                        data = lTrees[[paste0('dist_',sChrSource)]],
                        'tree'
                      )
                      
                      
                      ###################################################################################################
                      ### Cut the tree
                      ###################################################################################################
                      lTrees[[paste0('clust_comp_',sChrSource)]][['data']][,'y'] <- round(
                        lTrees[[paste0('clust_comp_',sChrSource)]][['data']][,'y'],
                        2
                      )
                      lTrees[[paste0('clust_selection_',sChrSource)]] <- lTrees[[paste0('clust_comp_',sChrSource)]][['data']]
                      lTrees[[paste0('clust_selection_',sChrSource)]][,'clusters'] <- as.numeric(as.character(
                        lTrees[[paste0('clust_selection_',sChrSource)]][,'clusters']
                      ))
                      lTrees[[paste0('clust_selection_',sChrSource)]][,'plateau'] <- c(
                        FALSE,
                        lTrees[[paste0('clust_selection_',sChrSource)]][-1,'y'] >= max(
                          lTrees[[paste0('clust_selection_',sChrSource)]][-1,'y']
                        )
                      )
                      lTrees[[paste0('clust_selection_',sChrSource)]][,'optimal'] <- (
                        lTrees[[paste0('clust_selection_',sChrSource)]][,'clusters'] == min(
                          lTrees[[paste0('clust_selection_',sChrSource)]][
                            lTrees[[paste0('clust_selection_',sChrSource)]][,'plateau'],
                            'clusters'])
                      )
                      lTrees[[paste0('clust_n_optimal_',sChrSource)]] <- (
                        lTrees[[paste0('clust_selection_',sChrSource)]][
                          lTrees[[paste0('clust_selection_',sChrSource)]][,'optimal'],'clusters']
                      )
                      lTrees[[paste0('clusters_',sChrSource)]] <- cutree(
                        lTrees[[paste0('tree_',sChrSource)]], 
                        k = lTrees[[paste0('clust_n_optimal_',sChrSource)]]
                      )
                      lTrees[[paste0('clusters_str_',sChrSource)]] <- as.character(cutree(
                        lTrees[[paste0('tree_',sChrSource)]], 
                        k = lTrees[[paste0('clust_n_optimal_',sChrSource)]]
                      ))
                      names(lTrees[[paste0('clusters_str_',sChrSource)]]) <- names(cutree(
                        lTrees[[paste0('tree_',sChrSource)]], 
                        k = lTrees[[paste0('clust_n_optimal_',sChrSource)]]
                      ))
                    }
                    lTrees[[paste0('clusters_str_',sChrSource)]] <- lTrees[[paste0('clusters_str_',sChrSource)]][
                      lTrees[[paste0('tree_',sChrSource)]][['order']]]
                    lTrees[['info']][,paste0('clusters_str_',sChrSource)] <- lTrees[[paste0('clusters_str_',sChrSource)]][rownames(lTrees[['info']])]
                    
                    
                    ###################################################################################################
                    ### Create the silhouette
                    ###################################################################################################
                    print(paste(Sys.time(),'Create the silhouette for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    lTrees[[paste0('silhouette_',sChrSource)]] <- silhouette(
                      x=lTrees[[paste0('clusters_',sChrSource)]], 
                      dist=lTrees[[paste0('dist_',sChrSource)]]
                    )
                    
                    
                    ###################################################################################################
                    ### Create the PCA
                    ###################################################################################################
                    print(paste(Sys.time(),'Create the PCA for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    lTrees[[paste0('PCA_',sChrSource)]] <- try(fviz_cluster(list(
                      data = lTrees[[paste0('dist_',sChrSource)]], 
                      cluster = lTrees[[paste0('clusters_',sChrSource)]] ))
                    )
                    
                    
                    ###################################################################################################
                    ### Set informative labels
                    ###################################################################################################
                    print(paste(Sys.time(),'Set tree labels for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    lTrees[[paste0('tree_',sChrSource)]][['order.lab.core_id']] <- c(
                      lTrees[[paste0('tree_',sChrSource)]][['order.lab']]
                    )
                    lTrees[[paste0('tree_',sChrSource)]][['order.lab.organ']] <- c(
                      lTrees[['info']][lTrees[[paste0('tree_',sChrSource)]][['order.lab']],'Organ']
                    )   
                    lTrees[[paste0('tree_',sChrSource)]][['order.lab.ARcnLAB']] <- paste(
                      'N.',lTrees[['info']][lTrees[[paste0('tree_',sChrSource)]][['order.lab']],'Figure_N'],'(',
                      lTrees[['info']][lTrees[[paste0('tree_',sChrSource)]][['order.lab']],'Sample_ID'],') - AR',
                      lTrees[['info']][lTrees[[paste0('tree_',sChrSource)]][['order.lab']],'ARcnLAB']
                    )
                    
                    
                    ###################################################################################################
                    ### Add tree dist
                    ###################################################################################################
                    try(oDataHCdist <- fSetLPdist(
                      data = oDataHCdist,
                      dist = lTrees[[paste0('dist_',sChrSource)]],
                      info = lTrees[['info']],
                      source = sChrSource,
                      color = lOrganCols,
                      title = paste(
                        'dist',
                        sSampleSel,
                        as.integer(nAutosomesBinSize),
                        sMeasureType,
                        sChrXRegion,
                        sChrSource,
                        sep='_'))
                    )
                    if (bAddExtraDist) {
                      try(oDataHCdist <- fSetLPdist(
                        oDataHCdist,
                        cophenetic(lTrees[[paste0('tree_',sChrSource)]]),
                        lTrees[['info']],
                        sChrSource,
                        lOrganCols,
                        paste(
                          'dist_cophenetic_',
                          sSampleSel,
                          as.integer(nAutosomesBinSize),
                          sMeasureType,
                          sChrXRegion,
                          sChrSource,
                          sep='_'))
                      )
                    }
                    
                    
                    ###################################################################################################
                    ### Set Title
                    ###################################################################################################
                    print(paste(Sys.time(),'Set title for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    lTrees[[paste0('title_',sChrSource)]] <- paste(
                      sTrialID,
                      sSampleSel,
                      sChrXRegion,
                      sMeasureType,
                      sChrSource,
                      '(',ifelse(sChrSource=='aut',as.integer(nAutosomesBinSize),as.integer(nChrXBinSize))/1000,'kb )'
                    )
                    
                    
                    ###################################################################################################
                    ### Plot the QC and the trees
                    ###################################################################################################
                    print(paste(Sys.time(),'Plot QC amd trees for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    try(plot(lTrees[[paste0('clust_comp_',sChrSource)]]+ ggtitle(lTrees[[paste0('title_',sChrSource)]])))
                    try(plot(lTrees[[paste0('silhouette_',sChrSource)]],main=lTrees[[paste0('title_',sChrSource)]]))
                    try(plot(lTrees[[paste0('PCA_',sChrSource)]] + ggtitle(lTrees[[paste0('title_',sChrSource)]])))
                    marOld <-  par('mar')
                    marNew <-  c(2,0,2,0)
                    par(mfrow=c(3,1),mar=marNew)
                    for (sLabelType in c('order.lab','order.lab.organ','order.lab.ARcnLAB','order.lab.core_id')) {try({
                      lTrees[[paste0('tree_',sChrSource)]][['order.lab']] <- c(
                        lTrees[[paste0('tree_',sChrSource)]][[sLabelType]]
                      )
                      if (sLabelType != 'order.lab.core_id') {
                        pltree(
                          lTrees[[paste0('tree_',sChrSource)]],
                          main=lTrees[[paste0('title_',sChrSource)]],
                          xlab=gsub('order.lab.','',sLabelType))
                        rect.hclust(
                          lTrees[[paste0('tree_',sChrSource)]], 
                          k = max(lTrees[[paste0('clusters_',sChrSource)]]), 
                          border = 2:(max(lTrees[[paste0('clusters_',sChrSource)]])+1))
                      }
                    })}
                    par(mfrow=c(1,1),mar=marOld)
                    
                    
                    ###################################################################################################
                    ### Plot data heatmap for chrX CNs
                    ###################################################################################################
                    print(paste(Sys.time(),'Plot chrX HM for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    if (sChrSource=='chrX') {
                      ###################################################################################################
                      ### Define colors and clusters
                      ###################################################################################################
                      if(max(lTrees[[paste0('clusters_',sChrSource)]])==1) {
                        oClustCol <- FALSE
                      } else {
                        oClustCol <- as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])
                      }
                      oColsRows <- data.frame(Region = lDataPlotAnno[
                        rownames(lTrees[[paste0('data_',sChrSource)]])],
                        stringsAsFactors=FALSE)
                      rownames(oColsRows) <- rownames(lTrees[[paste0('data_',sChrSource)]])
                      lColsRows <- list(
                        Region=lRegionsCol
                      )
                      
                      
                      ###################################################################################################
                      ### Define colors for all samples heatmap
                      ###################################################################################################
                      for (sColSV in c(
                        paste0('SV_ALL_',sTrialID),
                        grep(paste0('SV_',sTrialID),colnames(lTrees[['data_chrX_all']]),value=TRUE)
                      )) {
                        if (!all(is.na(lTrees[['data_chrX_all']][rownames(lTrees[[paste0('data_',sChrSource)]]),sColSV]))) {
                          oColsRows[,sColSV] <- lTrees[['data_chrX_all']][rownames(lTrees[[paste0('data_',sChrSource)]]),sColSV]
                          lColsRows <- c(lColsRows,list(tmp=lSVtypes[names(lSVtypes) %in% unique(na.omit(oColsRows[,sColSV]))]))
                          names(lColsRows)[names(lColsRows)=='tmp'] <- sColSV
                          if (sColSV==paste0('SV_ALL_',sTrialID)) {
                            oRegionCol[,sColSV] <- oColsRows[rownames(oRegionCol),sColSV]
                            lRegionCol <- c(lRegionCol,lColsRows[sColSV])
                          }
                        }
                      }
                      
                      
                      ###################################################################################################
                      ### Merege data for all samples heatmap
                      ###################################################################################################
                      if(all(lTrees[['info']][,'Run_ID']=='CASCADE')) {
                        sHMplotName <- paste0(
                          sSampleSel,'__',
                          sMeasureType,'_',
                          sChrXRegion,
                          '__A',as.integer(nAutosomesBinSize),
                          '__X',as.integer(nChrXBinSize)
                        )
                        print(paste('HM ROW HEADER',sHMplotName))
                        
                        
                        ###################################################################################################
                        ### Set info samples for final Heatmap
                        ###################################################################################################
                        lInfoSamples[[sHMplotName]] <- rbind(lInfoSamples[[sHMplotName]],lTrees[['info']]) 
                        
                        
                        ###################################################################################################
                        ### Set order cluster for final Heatmap
                        ###################################################################################################
                        print(paste(Sys.time(),'Set ordered labels for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                        lSamplesOrderCluster[[sHMplotName]] <- c(
                          lSamplesOrderCluster[[sHMplotName]],
                          lTrees[[paste0('clusters_str_',sChrSource)]]
                        )
                        
                        
                        ###################################################################################################
                        ### Set row anno. for final Heatmap
                        ###################################################################################################
                        if(is.null(lDataPlotChrX[[sHMplotName]])) {
                          lDataPlotChrXregionCol[[sHMplotName]] <- oRegionCol
                          lDataPlotChrX[[sHMplotName]] <- lTrees[[paste0('data_',sChrSource)]]
                        } else {
                          lDataPlotChrXregionCol[[sHMplotName]] <- merge(
                            lDataPlotChrXregionCol[[sHMplotName]],
                            oRegionCol,
                            by = 'row.names'
                          )
                          rownames(lDataPlotChrXregionCol[[sHMplotName]]) <- lDataPlotChrXregionCol[[sHMplotName]][,'Row.names']
                          lDataPlotChrXregionCol[[sHMplotName]] <- lDataPlotChrXregionCol[[sHMplotName]][,!colnames(lDataPlotChrXregionCol[[sHMplotName]]) %in% c('row.names','Row.names','Region.y'),drop=FALSE]  
                          colnames(lDataPlotChrXregionCol[[sHMplotName]])[colnames(lDataPlotChrXregionCol[[sHMplotName]])=='Region.x'] <- 'Region'
                          lDataPlotChrX[[sHMplotName]] <- merge(lDataPlotChrX[[sHMplotName]],
                                                                lTrees[[paste0('data_',sChrSource)]][,names(lTrees[[paste0('clusters_str_',sChrSource)]])],
                                                                by='row.names'
                          )
                          rownames(lDataPlotChrX[[sHMplotName]]) <- lDataPlotChrX[[sHMplotName]][,'Row.names']
                          lDataPlotChrX[[sHMplotName]] <- lDataPlotChrX[[sHMplotName]][,colnames(lDataPlotChrX[[sHMplotName]])!='Row.names']  
                        }
                      }
                      
                      
                      ###################################################################################################
                      ### Plot chrX heatmap
                      ###################################################################################################
                      lTrees[[paste0('data_',sChrSource)]] <-  lTrees[[paste0('data_',sChrSource)]][order(sapply(rownames( lTrees[[paste0('data_',sChrSource)]]),function(s) as.numeric(strsplit(s,'-')[[1]][2]))),]
                      oColsRows <- oColsRows[rownames(lTrees[[paste0('data_',sChrSource)]]),]
                      for (oRowSplit in list(
                        a=NULL
                        #b=lChrXRegionsGroupKB[rownames(lTrees[[paste0('data_',sChrSource)]])]
                      )) {
                        try(draw(Heatmap(
                          lTrees[[paste0('data_',sChrSource)]],
                          row_title = paste(
                            sTrialID,sSampleSel,
                            'CN profile of',sChrXRegion,
                            'by',sChrSource,sMeasureType,'(',ifelse(sChrSource=='aut',as.integer(nAutosomesBinSize),as.integer(nChrXBinSize))/1000,'kb)'),
                          name = paste0('CN',sMeasureType),
                          row_split = oRowSplit,
                          cluster_rows = FALSE,
                          cluster_columns = oClustCol,
                          show_row_names = FALSE,
                          show_column_names = TRUE,
                          col = lHMcol,
                          top_annotation = HeatmapAnnotation(
                            Organ = lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'Organ'],
                            ARcnLAB = lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ARcnLAB'],
                            clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                            ACEtc = ifelse(
                              lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ACE_cellularity_all']==1,
                              0,
                              lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ACE_cellularity_all']
                            ),
                            foo = anno_block(gp = gpar(fill = c('white','gray')),
                                             labels_gp = gpar(col = "black", fontsize = 16)),
                            col = list(
                              Organ = lHMcolOrgan,
                              ARcnLAB = lHMcolAR,
                              ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green')))),
                          left_annotation = rowAnnotation(
                            df=oColsRows,
                            col=lColsRows)
                        )))
                      }
                    }
                    
                    
                    ###################################################################################################
                    ### Plot data heatmap for Cor. dist
                    ###################################################################################################
                    print(paste(Sys.time(),'Plot cor dist HM for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    try(draw(Heatmap(
                      as.matrix(lTrees[[paste0('dist_',sChrSource)]])[
                        labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),
                        labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]]))
                      ],
                      row_title = paste(
                        sTrialID,sSampleSel,
                        'CN profile of',sChrXRegion,
                        'by',sChrSource,sMeasureType,'(',ifelse(sChrSource=='aut',as.integer(nAutosomesBinSize),as.integer(nChrXBinSize))/1000,'kb)'),
                      name = paste0('Cor. dist. for',sMeasureType),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      col = lHMcolCor,
                      top_annotation = HeatmapAnnotation(
                        Organ = lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'Organ'],
                        ARcnLAB = lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'ARcnLAB'],
                        clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                        ACEtc = ifelse(
                          lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'ACE_cellularity_all']==1,
                          0,
                          lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'ACE_cellularity_all']
                        ),
                        foo = anno_block(gp = gpar(fill = c('white','gray')),
                                         labels_gp = gpar(col = "black", fontsize = 16)),
                        col = list(
                          Organ = lHMcolOrgan,
                          ARcnLAB = lHMcolAR,
                          ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green')))),
                      left_annotation = rowAnnotation(
                        Organ = lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'Organ'],
                        ARcnLAB = lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'ARcnLAB'],
                        clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                        ACEtc = ifelse(
                          lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'ACE_cellularity_all']==1,
                          0,
                          lTrees[['info']][labels(as.dendrogram(lTrees[[paste0('tree_',sChrSource)]])),'ACE_cellularity_all']
                        ),
                        foo = anno_block(gp = gpar(fill = c('white','gray')),
                                         labels_gp = gpar(col = "black", fontsize = 16)),
                        col = list(
                          Organ = lHMcolOrgan,
                          ARcnLAB = lHMcolAR,
                          ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green'))))
                    )))
                    
                    
                    ###################################################################################################
                    ### Plot data heatmap for Cor. p
                    ###################################################################################################
                    print(paste(Sys.time(),'Plot p. val. HM for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    try(draw(Heatmap(
                      -log10(as.matrix(lTrees[[paste0('cor.p_',sChrSource)]])),
                      row_title = paste(
                        sTrialID,sSampleSel,
                        'CN profile of',sChrXRegion,
                        'by',sChrSource,sMeasureType,'(',ifelse(sChrSource=='aut',as.integer(nAutosomesBinSize),as.integer(nChrXBinSize))/1000,'kb)'),
                      name = paste0('-log10(Cor. p val.) for',sMeasureType),
                      cluster_rows = as.dendrogram(lTrees[[paste0('tree_',sChrSource)]]),
                      cluster_columns = as.dendrogram(lTrees[[paste0('tree_',sChrSource)]]),
                      col = lHMcolCorP,
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      top_annotation = HeatmapAnnotation(
                        Organ = lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'Organ'],
                        ARcnLAB = lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ARcnLAB'],
                        clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                        ACEtc = ifelse(
                          lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ACE_cellularity_all']==1,
                          0,
                          lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ACE_cellularity_all']
                        ),
                        foo = anno_block(gp = gpar(fill = c('white','gray')),
                                         labels_gp = gpar(col = "black", fontsize = 16)),
                        col = list(
                          Organ = lHMcolOrgan,
                          ARcnLAB = lHMcolAR,
                          ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green')))),
                      left_annotation = rowAnnotation(
                        Organ = lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'Organ'],
                        ARcnLAB = lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ARcnLAB'],
                        clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                        ACEtc = ifelse(
                          lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ACE_cellularity_all']==1,
                          0,
                          lTrees[['info']][colnames(lTrees[[paste0('data_',sChrSource)]]),'ACE_cellularity_all']
                        ),
                        foo = anno_block(gp = gpar(fill = c('white','gray')),
                                         labels_gp = gpar(col = "black", fontsize = 16)),
                        col = list(
                          Organ = lHMcolOrgan,
                          ARcnLAB = lHMcolAR,
                          ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green'))))
                    )))
                    
                    
                    ###################################################################################################
                    ### Create Phylo Trees
                    ###################################################################################################
                    lTrees[[paste0('trees_phylo_',sChrSource)]] <- try(fTreePhylo(lTrees[[paste0('dist_',sChrSource)]]))
                    if (!class(lTrees[[paste0('trees_phylo_',sChrSource)]])=='try-error' & bAddExtraDist) {
                      try(oDataHCdist <- fSetLPdist(
                        oDataHCdist,
                        cophenetic(lTrees[[paste0('trees_phylo_',sChrSource)]]),
                        lTrees[['info']],
                        sChrSource,
                        lOrganCols,
                        paste(
                          'dist_phylo',
                          sSampleSel,
                          as.integer(nAutosomesBinSize),
                          sMeasureType,
                          sChrXRegion,
                          sChrSource,
                          sep='_'))
                      )
                      fCompareTrees(
                        tree1 = lTrees[[paste0('tree_',sChrSource)]],
                        tree2 = lTrees[[paste0('trees_phylo_',sChrSource)]],
                        treeName1 = 'Correlation Tree',
                        treeName2 = 'NJ based on correlation dist'
                      )
                    }
                    
                    
                    ###################################################################################################
                    ### Create Phylo Trees based on TP hamming dist at different shift
                    ###################################################################################################
                    for (nShift in lShifts) {
                      lTrees[[paste0('dist_hamming_TPs',nShift,'_',sChrSource)]] <- try(fTreePhylo(
                        lTrees[[paste0('TPs',nShift,'_',sChrSource)]],
                        hamming.dist=TRUE,
                        add.root=FALSE,
                        return.dist=TRUE
                      ))
                      if (!class(lTrees[[paste0('dist_hamming_TPs',nShift,'_',sChrSource)]])=='try-error') {
                        try(oDataHCdist <- fSetLPdist(
                          oDataHCdist,
                          lTrees[[paste0('dist_hamming_TPs',nShift,'_',sChrSource)]],
                          lTrees[['info']],
                          sChrSource,
                          lOrganCols,
                          paste(
                            'dist_hamming_TPs',nShift,
                            sSampleSel,
                            as.integer(nAutosomesBinSize),
                            sMeasureType,
                            sChrXRegion,
                            sChrSource,
                            sep='_'))
                        )
                      }
                      
                      
                      ###################################################################################################
                      ### Create Phylo Trees based on TP hamming dist at different shift with / without parsimony
                      ###################################################################################################
                      for (bParsimony in lParsimonyOpts) {
                        lTrees[[paste0('tree_hamming_TPs',nShift,'_',bParsimony,'_',sChrSource)]] <- try(fTreePhylo(
                          lTrees[[paste0('TPs',nShift,'_',sChrSource)]],
                          hamming.dist=TRUE,
                          add.root=TRUE,
                          parsimony=bParsimony,
                          return.dist=FALSE
                        ))
                        if (!class(lTrees[[paste0('tree_hamming_TPs',nShift,'_',bParsimony,'_',sChrSource)]])=='try-error') {
                          fCompareTrees(
                            tree1 = lTrees[[paste0('tree_',sChrSource)]],
                            tree2 = lTrees[[paste0('tree_hamming_TPs',nShift,'_',bParsimony,'_',sChrSource)]],
                            treeName1 = 'Correlation Tree',
                            treeName2 = paste0('tree_hamming_TPs',nShift,'_',ifelse(bParsimony,'parsimony',''),'_',sChrSource),
                            hamming.data = lTrees[[paste0('TPs',nShift,'_',sChrSource)]],
                            info = lTrees[['info']]
                          )
                          if (bAddExtraDist) {
                            try(oDataHCdist <- fSetLPdist(
                              oDataHCdist,
                              cophenetic(lTrees[[paste0('tree_hamming_TPs',nShift,'_',bParsimony,'_',sChrSource)]]),
                              lTrees[['info']],
                              sChrSource,
                              lOrganCols,
                              paste0(
                                'dist_hamming_cophenetic_TPs',nShift,'_',
                                ifelse(bParsimony,'parsimony',''),'_',
                                sSampleSel,'_',
                                as.integer(nAutosomesBinSize),'_',
                                sMeasureType,'_',
                                sChrXRegion,'_',
                                sChrSource
                              ))
                            )
                          }
                        }
                      }
                    }
                    
                    
                    ###################################################################################################
                    ### Plot data heatmap for Hamming Dist TPs2
                    ###################################################################################################
                    print(paste(Sys.time(),'Plot cor dist Hamming Dist TPs2 for',sTrialID,sMeasureType,sChrXRegion,sChrSource))
                    for (bParsimony in lParsimonyOpts) {
                      oMclustTree <- dendextend::prune(as.dendrogram(lTrees[[paste0('tree_hamming_TPs2_',bParsimony,'_',sChrSource)]]),'root')
                      try(draw(Heatmap(
                        as.matrix(lTrees[[paste0('dist_hamming_TPs2_',sChrSource)]])[
                          labels(oMclustTree),labels(oMclustTree)  
                        ],
                        row_title = paste(
                          sTrialID,sSampleSel,
                          'CN profile of',sChrXRegion,
                          'by',sChrSource,sMeasureType,'(',ifelse(sChrSource=='aut',as.integer(nAutosomesBinSize),as.integer(nChrXBinSize))/1000,'kb)'),
                        name = paste0('Hamming ',bParsimony,'TPs2 for',sMeasureType),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        col = lHMcolCor,
                        top_annotation = HeatmapAnnotation(
                          Organ = lTrees[['info']][labels(oMclustTree),'Organ'],
                          ARcnLAB = lTrees[['info']][labels(oMclustTree),'ARcnLAB'],
                          clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                          foo = anno_block(gp = gpar(fill = c('white','gray')),
                                           labels_gp = gpar(col = "black", fontsize = 16)),
                          col = list(
                            Organ = lHMcolOrgan,
                            ARcnLAB = lHMcolAR,
                            ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green')))),
                        left_annotation = rowAnnotation(
                          Organ = lTrees[['info']][labels(oMclustTree),'Organ'],
                          ARcnLAB = lTrees[['info']][labels(oMclustTree),'ARcnLAB'],
                          clusters = lTrees[[paste0('clusters_str_',sChrSource)]],
                          foo = anno_block(gp = gpar(fill = c('white','gray')),
                                           labels_gp = gpar(col = "black", fontsize = 16)),
                          col = list(
                            Organ = lHMcolOrgan,
                            ARcnLAB = lHMcolAR,
                            ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green'))))
                      )))
                    }
                  })}
                  
                  
                  ###################################################################################################
                  ### Plot the tanglegram to compare chrX and aut
                  ###################################################################################################
                  print(paste(Sys.time(),'Plot the tanglegram to compare chrX and aut for',sTrialID,sMeasureType,sChrXRegion))
                  par(mfrow=c(1,2))
                  try(pltree(lTrees[['tree_chrX']],"Dendrogram by chrX"))
                  try(pltree(lTrees[['tree_aut']],"Dendrogram of auosomes"))
                  par(mfrow=c(1,1))
                  try(tanglegram(
                    untangle(
                      as.dendrogram(lTrees[['tree_chrX']]),
                      as.dendrogram(lTrees[['tree_aut']]),
                      method="step2side", 
                    )[[1]],
                    untangle(
                      as.dendrogram(lTrees[['tree_chrX']]),
                      as.dendrogram(lTrees[['tree_aut']]),
                      method="step2side", 
                    )[[2]],
                    highlight_distinct_edges = FALSE, # Turn-off dashed lines
                    common_subtrees_color_lines = FALSE, # Turn-off line colors
                    common_subtrees_color_branches = TRUE, # Color common branches 
                    margin_inner = 10,
                    main_left = 'chrX',
                    main_right = 'autosomes',
                    main = paste("entanglement =", round(entanglement(dendlist(
                      as.dendrogram(lTrees[['tree_chrX']]), 
                      as.dendrogram(lTrees[['tree_aut']])
                    )), 2))
                  ))
                  
                  
                  ###################################################################################################
                  ### Plot the tanglegram to compare chrX and aut by hamming dist
                  ###################################################################################################
                  try(tanglegram(
                    untangle(
                      as.dendrogram(lTrees[['tree_hamming_TPs2_TRUE_chrX']]),
                      as.dendrogram(lTrees[['tree_hamming_TPs2_TRUE_aut']]),
                      method="step2side", 
                    )[[1]],
                    untangle(
                      as.dendrogram(lTrees[['tree_hamming_TPs2_TRUE_chrX']]),
                      as.dendrogram(lTrees[['tree_hamming_TPs2_TRUE_aut']]),
                      method="step2side", 
                    )[[2]],
                    highlight_distinct_edges = FALSE, # Turn-off dashed lines
                    common_subtrees_color_lines = FALSE, # Turn-off line colors
                    common_subtrees_color_branches = TRUE, # Color common branches 
                    margin_inner = 10,
                    main_left = 'chrX hamming TPs2',
                    main_right = 'autosomes hamming TPs2',
                    main = paste("entanglement =", round(entanglement(dendlist(
                      as.dendrogram(lTrees[['tree_hamming_TPs2_TRUE_chrX']]), 
                      as.dendrogram(lTrees[['tree_hamming_TPs2_TRUE_aut']])
                    )), 2))
                  ))
                  
                  
                  ###################################################################################################
                  ### Plot all aut PCA with chr clusters
                  ###################################################################################################
                  print(paste(Sys.time(),'Plot all aut PCA with chr clusters for',sTrialID,sMeasureType,sChrXRegion))
                  try({
                    oPlot <- ggplot(
                      merge(
                        lTrees[['info']][,c('Core_ID','Sample_ID','Organ','Organ.Color','ARcnLAB')],
                        merge(
                          lTrees[['PCA_aut']][['data']],
                          data.frame(
                            name = names(lTrees[['clusters_str_chrX']]),
                            cluster=lTrees[['clusters_str_chrX']],
                            stringsAsFactors = FALSE
                          ),
                          all.x=TRUE,
                          by='name'
                        ),
                        by.x='Core_ID',
                        by.y='name'
                      )
                    )
                    oPlot <- oPlot + geom_point(aes(
                      x=x,
                      y=y,
                      size = ARcnLAB,
                      fill=cluster.y
                    ),
                    shape = 21
                    )
                    oPlot <- oPlot + geom_text_repel(aes(
                      label = Organ,
                      x=x,
                      y=y
                    ))
                    oPlot <- oPlot + ggtitle(lTrees[['title_aut']])
                    oPlot <- oPlot + theme_bw()
                    plot(oPlot)
                  })
                  
                  
                  ###################################################################################################
                  ### Export Cluster data
                  ###################################################################################################
                  print(paste(Sys.time(),'Export Cluster data for',sTrialID,sMeasureType,sChrXRegion))
                  lTrees[['clusters_merge']] <- data.frame(
                    data.frame(clusters_aut= lTrees[['clusters_aut']])
                  )
                  lTrees[['clusters_merge']][,'Row.names'] <- rownames(
                    lTrees[['clusters_merge']]
                  )
                  lTrees[['clusters_merge']] <- merge(
                    lTrees[['clusters_merge']],
                    merge(
                      data.frame(clusters_aut=lTrees[['clusters_aut']]),
                      data.frame(clusters_chrX=lTrees[['clusters_chrX']]),
                      by = 'row.names',
                      all=TRUE
                    ),
                    all=TRUE
                  )
                  lTrees[['clusters_merge']] <- merge(
                    lTrees[['info']][,c('Core_ID','Organ')],
                    lTrees[['clusters_merge']],
                    all.y=TRUE,
                    by.x='Core_ID',
                    by.y='Row.names'
                  )
                  write.table(
                    lTrees[['info']],
                    file=paste0(
                      sPathCASCADE,
                      paste0('/S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                      '/Models__',sTrialID,'__',sSampleSel,'_',sMeasureType,'__',sChrXRegion,'__A',as.integer(nAutosomesBinSize),'__X',as.integer(nChrXBinSize),'.txt'), 
                    sep="\t", 
                    append=FALSE , 
                    row.names=FALSE, 
                    col.names=TRUE
                  )
                  
                  
                  ###################################################################################################
                  ### Save Cluster data objexts
                  ###################################################################################################
                  print(paste(Sys.time(),'Save Cluster data objexts for',sTrialID,sMeasureType,sChrXRegion))
                  save(lTrees,file=paste0(
                    sPathCASCADE,
                    paste0('/S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                    '/Models__',sTrialID,'__',sSampleSel,'_',sMeasureType,'__',sChrXRegion,'__A',as.integer(nAutosomesBinSize),'__X',as.integer(nChrXBinSize),'.Rdata'))
                  
                  
                  ###################################################################################################
                  ### Close PDF
                  ###################################################################################################
                  try(dev.off());try(dev.off());try(dev.off());try(dev.off())
                  
                  
                  ###################################################################################################
                  ### Get TPs distributions data
                  ###################################################################################################
                  oDataPlot <- NULL
                  lRegionGrop <- list()
                  for (sGroup in c(
                    'all',
                    unique(lTrees[['info']][,'Organ']),
                    unique(lTrees[['info']][,'clusters_str_aut'])
                  )) {try({
                    if (sGroup == 'all') {
                      oDataPlot <- rbind(oDataPlot,data.frame(
                        group_type = 'all',
                        group = sGroup,
                        region = rownames(lTrees[['TPs2_aut']]),
                        chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
                        sum = rowSums(lTrees[['TPs2_aut']]),
                        samples = ncol(lTrees[['TPs2_aut']]),
                        sample_count = rowSums(lTrees[['TPs2_aut']]),
                        freq = rowSums(lTrees[['TPs2_aut']]) / ncol(lTrees[['TPs2_aut']]),
                        count = 1,
                        cumsum_dec = 0,
                        cumsum_inc = 0,
                        stringsAsFactors = FALSE        
                      ))
                      oDataPlot <- rbind(oDataPlot,data.frame(
                        group_type = 'all',
                        group = sGroup,
                        region = rownames(lTrees[['TPs2_chrX']]),
                        chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
                        sum = rowSums(lTrees[['TPs2_chrX']]),
                        samples = ncol(lTrees[['TPs2_chrX']]),
                        sample_count = rowSums(lTrees[['TPs2_chrX']]),
                        freq = rowSums(lTrees[['TPs2_chrX']]) / ncol(lTrees[['TPs2_chrX']]),
                        count = 1,
                        cumsum_dec = 0,
                        cumsum_inc = 0,
                        stringsAsFactors = FALSE        
                      ))
                    }
                    if (sGroup %in% unique(lTrees[['info']][,'Organ'])) {
                      oDataPlot <- rbind(oDataPlot,data.frame(
                        group_type = 'Organ',
                        group = sGroup,
                        region = rownames(lTrees[['TPs2_aut']]),
                        chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
                        sum = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]),
                        samples = sum(colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID']),
                        sample_count = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]),
                        freq = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]) / 
                          sum(colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID']),
                        count = 1,
                        cumsum_dec = 0,
                        cumsum_inc = 0,
                        stringsAsFactors = FALSE        
                      ))
                      oDataPlot <- rbind(oDataPlot,data.frame(
                        group_type = 'Organ',
                        group = sGroup,
                        region = rownames(lTrees[['TPs2_chrX']]),
                        chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
                        sum = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]),
                        samples = sum(colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID']),
                        sample_count = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]),
                        freq = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]) / 
                          sum(colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID']),
                        count = 1,
                        cumsum_dec = 0,
                        cumsum_inc = 0,
                        stringsAsFactors = FALSE        
                      ))
                    }
                    if (sGroup %in% unique(lTrees[['info']][,'clusters_str_aut'])) {
                      oDataPlot <- rbind(oDataPlot,data.frame(
                        group_type = 'Cluster',
                        group = sGroup,
                        region = rownames(lTrees[['TPs2_aut']]),
                        chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
                        sum = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]),
                        samples = sum(colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']),
                        sample_count = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]),
                        freq = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]) / 
                          sum(colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']),
                        count = 1,
                        cumsum_dec = 0,
                        cumsum_inc = 0,
                        stringsAsFactors = FALSE        
                      ))
                      oDataPlot <- rbind(oDataPlot,data.frame(
                        group_type = 'Cluster',
                        group = sGroup,
                        region = rownames(lTrees[['TPs2_chrX']]),
                        chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
                        sum = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]),
                        samples = sum(colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']),
                        sample_count = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]),
                        freq = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]) / 
                          sum(colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']),
                        count = 1,
                        cumsum_dec = 0,
                        cumsum_inc = 0,
                        stringsAsFactors = FALSE        
                      ))
                    }
                    lRegionGrop <- c(lRegionGrop,list(tmp=unique(oDataPlot[
                      oDataPlot[,'group']==sGroup &
                        oDataPlot[,'freq'] >= 0.8,
                      'region'
                    ])))
                    names(lRegionGrop)[names(lRegionGrop)=='tmp'] <- sGroup
                    oDataPlot <- oDataPlot[order(oDataPlot[,'sum'],decreasing=TRUE),]
                    oDataPlot[oDataPlot[,'group']==sGroup,'cumsum_dec'] <- cumsum(oDataPlot[oDataPlot[,'group']==sGroup,'count'])
                    oDataPlot <- oDataPlot[order(oDataPlot[,'sum'],decreasing=FALSE),]
                    oDataPlot[oDataPlot[,'group']==sGroup,'cumsum_inc'] <- cumsum(oDataPlot[oDataPlot[,'group']==sGroup,'count'])
                  })}
                  
                  
                  ###################################################################################################
                  ### Define colors
                  ###################################################################################################
                  lOrganCols <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ.Color']
                  names(lOrganCols) <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ']
                  
                  
                  ###################################################################################################
                  ### Get TPs distributions data
                  ###################################################################################################
                  pdf(paste0(
                    sPathCASCADE,
                    paste0('/S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                    '/TrasitionPoints__',sTrialID,'__',sSampleSel,'_',sMeasureType,'__',sChrXRegion,'__A',as.integer(nAutosomesBinSize),'__X',as.integer(nChrXBinSize),'.pdf'),
                    paper='special',width=16,height=8
                  )
                  
                  
                  ###################################################################################################
                  ### Plot densities by Organ and Cluster
                  ###################################################################################################
                  try(grid.arrange(
                    ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('all','Organ'),]) + 
                      geom_density_ridges(
                        aes(x=freq,y=interaction(group,group_type),fill=group)
                      ) + 
                      theme_bw() + 
                      theme(axis.text.x=element_text(angle=45, hjust=1)) + 
                      ggtitle(sTrialID) + 
                      xlab('Fraction of samples') + 
                      ylab('Samples by Organ'),
                    ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('all','Organ'),]) + 
                      geom_density_ridges(
                        aes(x=sample_count,y=interaction(group,group_type),fill=group)
                      ) + 
                      theme_bw() + 
                      theme(axis.text.x=element_text(angle=45, hjust=1)) + 
                      ggtitle(sTrialID) + 
                      xlab('Number of samples') + 
                      ylab('Samples by Organ'),
                    ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('all','Cluster'),]) + 
                      geom_density_ridges(
                        aes(x=freq,y=interaction(group,group_type),fill=group)
                      ) + 
                      theme_bw() + 
                      theme(axis.text.x=element_text(angle=45, hjust=1)) + 
                      ggtitle(sTrialID) + 
                      xlab('Fraction of samples') + 
                      ylab('Samples by Cluster'),
                    ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('all','Cluster'),]) + 
                      geom_density_ridges(
                        aes(x=sample_count,y=interaction(group,group_type),fill=group)
                      ) + 
                      theme_bw() + 
                      theme(axis.text.x=element_text(angle=45, hjust=1)) + 
                      ggtitle(sTrialID) + 
                      xlab('Number of samples') + 
                      ylab('Samples by Cluster'),
                    ncol=2,
                    nrow=2,
                    layout_matrix=matrix(nrow=2,ncol=2,c(1,2,3,4))
                  ))
                  
                  
                  ###################################################################################################
                  ### Plot frequencies by Organ and Cluster
                  ###################################################################################################
                  try(grid.arrange(
                    ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('all','Organ'),]) + 
                      stat_smooth(aes(x=sample_count,y=cumsum_dec,color=group,fill=group,group=group)) + 
                      scale_x_reverse() + 
                      theme_bw() + 
                      ggtitle(sTrialID) + 
                      xlab('Number of samples') + 
                      ylab('Cumulative number of transition points'),
                    ggplot(oDataPlot[oDataPlot[,'freq'] > 0 & oDataPlot[,'group_type'] %in% c('all','Organ'),]) + 
                      stat_smooth(aes(x= 1 / freq,y=cumsum_dec,color=group,fill=group,group=group)) + 
                      theme_bw() + 
                      ggtitle(sTrialID) + 
                      xlab('1 / Fraction of samples') + 
                      ylab('Cumulative number of transition points'),
                    ncol=1,
                    nrow=2,
                    layout_matrix=matrix(nrow=2,ncol=1,c(1,2))
                  ))
                  
                  
                  ###################################################################################################
                  ### Plot Transition points upset
                  ###################################################################################################
                  try(oPlot <- upset(fromList(lRegionGrop[names(lRegionGrop) %in% unique(lTrees[['info']][,'Organ'])]), nsets = sum(names(lRegionGrop) %in% unique(lTrees[['info']][,'Organ'])), order.by = "freq"))
                  try(print(oPlot))
                  
                  
                  ###################################################################################################
                  ### Close PDF
                  ###################################################################################################
                  try(dev.off());try(dev.off());try(dev.off());try(dev.off())
                })}
                
                
                ###################################################################################################
                ### Close pending PDFs
                ###################################################################################################
                try(dev.off());try(dev.off());try(dev.off());try(dev.off())
              }
            })}
          }
          
          
          ###################################################################################################
          ### If Archival samaple is available run archival analysis
          ###################################################################################################
          if (bRunArchivalAnalyses & nrow(lTrees[['info_arch']])>0) {
            ###################################################################################################
            ### Loop through all archival samples
            ###################################################################################################
            sSampleArchival <- lTrees[['info_arch']][,'Core_ID'][2]
            for (sSampleArchival in lTrees[['info_arch']][,'Core_ID']) {
              oResComparison <- list()
              oSampleComarison <- list()
              
              
              ###################################################################################################
              ### Identify Archival sample
              ###################################################################################################
              nSampleArchival <- (1:length(sampleNames(oSM[['data']])))[
                sampleNames(oSM[['data']])==sSampleArchival]
              sPrefix <- 'Other'
              if (lTrees[['info_arch']][sSampleArchival,'Organ']=='Plasma') sPrefix <- 'Plasma'
              if (
                lTrees[['info_arch']][sSampleArchival,'Organ']=='Prostate' |
                regexpr('prostate',lTrees[['info_arch']][sSampleArchival,'Tissue.Site'])>-1 |
                regexpr('Prostate',lTrees[['info_arch']][sSampleArchival,'Tissue.Site'])>-1
              ) {
                sPrefix <- 'Prostate'
              }
              if (regexpr('archival',lTrees[['info_arch']][sSampleArchival,'Tissue.Site'])>-1) sPrefix <- 'Archival'
              
              
              ###################################################################################################
              ### Identify Archival sample segments
              ###################################################################################################
              for (sSampleType in c('p','c')) {
                for (nShift in lShifts) {
                  oResComparison[[paste0('S',nShift,sSampleType)]] <- rbind(
                    cbind(
                      unique(oRes[
                        !is.na(oRes[,paste0(sSampleArchival,'_segment')]) &
                          !is.na(oRes[,paste0(sSampleArchival,'_segmented')]) &
                          oRes[,paste0(sSampleArchival,'_cn_type_start')] != 0,
                        paste(sSampleArchival,
                              c(
                                'segment','segment_chr','segment_start','segment_end',
                                'segmented','Segment_Mean','Copies',
                                'cn_type_start','cn_type_end'
                              ),
                              sep='_'
                        )
                      ]),
                      data.frame(archival_sample_id = sSampleArchival, archival_segment_side='S', stringsAsFactors=FALSE)
                    ),
                    cbind(
                      unique(oRes[
                        !is.na(oRes[,paste0(sSampleArchival,'_segment')]) &
                          !is.na(oRes[,paste0(sSampleArchival,'_segmented')]) &
                          oRes[,paste0(sSampleArchival,'_cn_type_end')] != 0,
                        paste(sSampleArchival,
                              c(
                                'segment','segment_chr','segment_start','segment_end',
                                'segmented','Segment_Mean','Copies',
                                'cn_type_start','cn_type_end'
                              ),
                              sep='_'
                        )
                      ]),
                      data.frame(archival_sample_id = sSampleArchival, archival_segment_side='E', stringsAsFactors=FALSE)
                    )
                  )
                  oResComparison[[paste0('S',nShift,sSampleType)]][,
                                                                   'archival_segment'] <-  oResComparison[[paste0('S',nShift,sSampleType)]][,
                                                                                                                                            paste0(sSampleArchival,'_segment')
                                                                   ]
                  oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_type'] <- '-'
                  oResComparison[[paste0('S',nShift,sSampleType)]][
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side']=='S',
                    'archival_segment_type'] <- ifelse(
                      oResComparison[[paste0('S',nShift,sSampleType)]][
                        oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side']=='S',
                        paste0(sSampleArchival,'_cn_type_start')]>0,
                      'U','D')
                  oResComparison[[paste0('S',nShift,sSampleType)]][
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side']=='E',
                    'archival_segment_type'] <- ifelse(
                      oResComparison[[paste0('S',nShift,sSampleType)]][
                        oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side']=='E',
                        paste0(sSampleArchival,'_cn_type_end')]>0,
                      'U','D')
                  rownames(oResComparison[[paste0('S',nShift,sSampleType)]]) <- paste(
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment'],
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'],
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_type'],
                    sep='_')
                }
              }
              
              
              ###################################################################################################
              ### Create comparisons PDF
              ###################################################################################################
              pdf(file.path(
                sPathCASCADE,
                paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                paste0(sPrefix,'__',sTrialID,'__',sSampleSel,'__',as.integer(nAutosomesBinSize),'__vs__',sSampleArchival,'.pdf')),
                paper='special',
                width=16,
                height=8
              )
              
              
              ###################################################################################################
              ### Identify commont segments with mets
              ###################################################################################################
              for (sSampleCurrent in lTrees[['info_all']][,'Core_ID']) {
                nSampleCurrent <- (1:length(sampleNames(oSM[['data']])))[
                  sampleNames(oSM[['data']])==sSampleCurrent]
                print(paste(sSampleCurrent,nSampleCurrent))
                
                
                ###################################################################################################
                ### Plot common segments 
                ###################################################################################################
                oSampleComarison[[sSampleCurrent]] <- twosamplecompare(
                  oSM[['data']], 
                  index1 = nSampleArchival, 
                  ploidy1 = oSM[['bestparam']][[sSampleArchival]][1,'ploidy'],
                  cellularity1 = oSM[['bestparam']][[sSampleArchival]][1,'cellularity'], 
                  standard1 = 1, 
                  name1 = paste(
                    sSampleArchival,
                    oDataStats[oDataStats[,'Core_ID']==sSampleArchival,'Sample_ID'],
                    oDataStats[oDataStats[,'Core_ID']==sSampleArchival,'Tissue.Site']
                  ), 
                  index2 = nSampleCurrent,
                  ploidy2 = oSM[['bestparam']][[sSampleCurrent]][1,'ploidy'], 
                  cellularity2 = oSM[['bestparam']][[sSampleCurrent]][1,'cellularity'], 
                  standard2 = 1, 
                  name2 = paste(
                    sSampleCurrent,
                    oDataStats[oDataStats[,'Core_ID']==sSampleCurrent,'Sample_ID'],
                    oDataStats[oDataStats[,'Core_ID']==sSampleCurrent,'Tissue.Site']
                  ),
                  equalsegments = FALSE, 
                  altmethod = FALSE, 
                  onlyautosomes = TRUE, 
                  showcorrelation = TRUE)
                plot(oSampleComarison[[sSampleCurrent]][['compareplot']])
                
                
                ###################################################################################################
                ### Identify common segments with mets at differnt shifts
                ###################################################################################################
                sSampleType <- 'p'
                for (nShift in lShifts) {
                  ###################################################################################################
                  ### Identify segments in the met sample
                  ###################################################################################################
                  oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]] <- unique(oRes[
                    !is.na(oRes[,paste0(sSampleCurrent,'_segment')]) &
                      !is.na(oRes[,paste0(sSampleCurrent,'_segmented')]) &
                      (
                        oRes[,paste0(sSampleCurrent,'_cn_type_start')] != 0 |
                          oRes[,paste0(sSampleCurrent,'_cn_type_end')] != 0
                      ),
                    paste(sSampleCurrent,
                          c(
                            'segment','segment_chr','segment_start','segment_end',
                            'segmented','Segment_Mean','Copies',
                            'cn_type_start','cn_type_end'
                          ),
                          sep='_'
                    )
                  ])
                  
                  
                  ###################################################################################################
                  ### If there are testable TPS
                  ###################################################################################################
                  if (nrow(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]])>0) {
                    oTmp <- NULL
                    
                    
                    ###################################################################################################
                    ### Check each TP
                    ###################################################################################################
                    lArchivalSegments <- unique(oResComparison[[paste0('S',nShift,sSampleType)]][order(
                      as.integer(gsub('X','23',gsub('Y','24',oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]))),
                      as.integer(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_start')])
                    ),'archival_segment'])
                    for (sArchivalSegment in lArchivalSegments) {
                      ###################################################################################################
                      ### Create dummy df
                      ###################################################################################################
                      oTmpEmpty <- data.frame(
                        archival_segment = sArchivalSegment,
                        archival_segment_side = c('S','S','E','E'),
                        archival_segment_type = c('U','D','U','D'),
                        segment = NA,
                        segment_chr = NA,
                        segment_start = NA,
                        segment_end = NA,
                        segmented = NA,
                        Segment_Mean = NA,
                        Copies = NA,
                        cn_type_start = NA,
                        cn_type_end = NA,
                        stringsAsFactors = FALSE
                      )
                      colnames(oTmpEmpty)[colnames(oTmpEmpty) %in%  c(
                        'segment','segment_chr','segment_start','segment_end',
                        'segmented','Segment_Mean','Copies',
                        'cn_type_start','cn_type_end'
                      )] <- paste(sSampleCurrent,colnames(oTmpEmpty)[colnames(oTmpEmpty) %in%  c(
                        'segment','segment_chr','segment_start','segment_end',
                        'segmented','Segment_Mean','Copies',
                        'cn_type_start','cn_type_end'
                      )],
                      sep='_'
                      )
                      
                      
                      ###################################################################################################
                      ### Add start matached TP
                      ###################################################################################################
                      for (sArchivalSide in unique(oResComparison[[paste0('S',nShift,sSampleType)]][
                        oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment']==sArchivalSegment,
                        'archival_segment_side'])) {
                        for (sArchivalType in unique(oResComparison[[paste0('S',nShift,sSampleType)]][
                          oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment']==sArchivalSegment &
                          oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side']==sArchivalSide,
                          'archival_segment_type'])) {
                          nRowArch <- min(which(
                            !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]) &
                              !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))]) &
                              !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))]) &
                              oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment'] == sArchivalSegment &
                              oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'] == sArchivalSide &
                              oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_type'] == sArchivalType
                          ))
                          lRowSel <- which(
                            !is.na(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segment_chr')]) &
                              !is.na(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))]) &
                              !is.na(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segment_',ifelse(sArchivalSide=='S','start','end'))]) &
                              !is.na(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segmented')]) &
                              !is.na(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_Segment_Mean')]) &
                              !is.na(oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_Copies')]) &
                              oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segment_chr')] == oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_chr')] &
                              oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))] == oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))] &
                              oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segment_',ifelse(sArchivalSide=='S','start','end'))] >= oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))] - (nShift * nAutosomesBinSize) &
                              oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segment_',ifelse(sArchivalSide=='S','start','end'))] <= oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))] + (nShift * nAutosomesBinSize)
                          )
                          if (length(lRowSel)>0) {
                            if (length(lRowSel)>1) {
                              lRowSel <- lRowSel[abs(
                                oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][lRowSel,paste0(sSampleCurrent,'_segment_',ifelse(sArchivalSide=='S','start','end'))] - 
                                  oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))]
                              ) == min(abs(
                                oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][lRowSel,paste0(sSampleCurrent,'_segment_',ifelse(sArchivalSide=='S','start','end'))] - 
                                  oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))]
                              ))][1]
                            }
                            oTmp <- rbind(
                              oTmp,
                              cbind(
                                data.frame(
                                  archival_segment = sArchivalSegment,
                                  archival_segment_side = sArchivalSide,
                                  archival_segment_type = sArchivalType,
                                  stringsAsFactors = FALSE
                                ),
                                oSampleComarison[[sSampleCurrent]][['segments_current']][[paste0('S',nShift,sSampleType)]][lRowSel,]
                              )
                            )
                          } else {
                            oTmp <- rbind(
                              oTmp,
                              oTmpEmpty[
                                oTmpEmpty[,'archival_segment_side']==sArchivalSide &
                                  oTmpEmpty[,'archival_segment_type']==sArchivalType,
                              ]
                            )
                          }
                        }
                      }
                    }
                    if (!is.null(oTmp)) {
                      oResComparison[[paste0('S',nShift,sSampleType)]] <- merge(
                        oResComparison[[paste0('S',nShift,sSampleType)]],
                        oTmp,
                        all.x=TRUE)
                    }
                    print(paste(nShift,dim(oTmp)))
                  }
                  print(paste(nShift,dim(oResComparison[[paste0('S',nShift,sSampleType)]])))
                  oResComparison[[paste0('S',nShift,sSampleType)]] <- oResComparison[[paste0('S',nShift,sSampleType)]][order(
                    as.integer(gsub('X','23',gsub('Y','24',oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]))),
                    as.integer(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_start')])
                  ),]
                  rownames(oResComparison[[paste0('S',nShift,sSampleType)]]) <- paste(
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment'],
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'],
                    oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_type'],
                    sep='_')
                }
              }
              
              
              ###################################################################################################
              ### Closet PDF
              ###################################################################################################
              try(dev.off());try(dev.off());try(dev.off());try(dev.off())
              
              
              ###################################################################################################
              ### Validate archival breakpoinsts
              ###################################################################################################
              sSampleType <- 'c'
              for (sPatientCheck in setdiff(unique(oDataStats[oDataStats[,'Run_ID']=='CASCADE','Trial.ID']),sTrialID)) {
                print(paste(Sys.time(),'Load Extra patient',sPatientCheck,'CN data for ',sPrefix,' validation in',sTrialID))
                if (file.exists(file.path(
                  sPathCASCADE,
                  'details_bypatient',
                  paste0(sPatientCheck,'__',as.integer(nAutosomesBinSize),'.txt')))) {
                  
                  
                  ###################################################################################################
                  ### Read sample data
                  ###################################################################################################
                  oResCheck <- read.delim(
                    file.path(
                      sPathCASCADE,
                      'details_bypatient',
                      paste0(sPatientCheck,'__',as.integer(nAutosomesBinSize),'.txt')),
                    sep='\t',
                    header=TRUE,
                    stringsAsFactors=FALSE
                  )
                  rownames(oResCheck) <- oResCheck[,'region']
                  
                  
                  ###################################################################################################
                  ### Identify common breakpoints in external samples
                  ###################################################################################################
                  for (sSampleCheck in gsub('_Copies','',grep('_Copies',colnames(oResCheck),value=TRUE),fixed=TRUE)) {
                    for (nShift in lShifts) {
                      oTmp <- NULL
                      lArchivalSegments <- unique(oResComparison[[paste0('S',nShift,sSampleType)]][order(
                        as.integer(gsub('X','23',gsub('Y','24',oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]))),
                        as.integer(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_start')])
                      ),'archival_segment'])
                      for (sArchivalSegment in lArchivalSegments) {
                        ###################################################################################################
                        ### Create dummy df
                        ###################################################################################################
                        oTmpEmpty <- data.frame(
                          archival_segment = sArchivalSegment,
                          archival_segment_side = c('S','S','E','E'),
                          archival_segment_type = c('U','D','U','D'),
                          segment = NA,
                          segment_chr = NA,
                          segment_start = NA,
                          segment_end = NA,
                          segmented = NA,
                          Segment_Mean = NA,
                          Copies = NA,
                          cn_type_start = NA,
                          cn_type_end = NA,
                          stringsAsFactors = FALSE
                        )
                        colnames(oTmpEmpty)[colnames(oTmpEmpty) %in%  c(
                          'segment','segment_chr','segment_start','segment_end',
                          'segmented','Segment_Mean','Copies',
                          'cn_type_start','cn_type_end'
                        )] <- paste(sSampleCheck,colnames(oTmpEmpty)[colnames(oTmpEmpty) %in%  c(
                          'segment','segment_chr','segment_start','segment_end',
                          'segmented','Segment_Mean','Copies',
                          'cn_type_start','cn_type_end'
                        )],
                        sep='_')
                        
                        
                        ###################################################################################################
                        ### Add start matached TP
                        ###################################################################################################
                        for (sArchivalSide in unique(oResComparison[[paste0('S',nShift,sSampleType)]][
                          oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment']==sArchivalSegment,
                          'archival_segment_side'])) {
                          for (sArchivalType in unique(oResComparison[[paste0('S',nShift,sSampleType)]][
                            oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment']==sArchivalSegment &
                            oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side']==sArchivalSide,
                            'archival_segment_type'])) {
                            nRowArch <- min(which(
                              !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]) &
                                !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))]) &
                                !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))]) &
                                oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment'] == sArchivalSegment &
                                oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'] == sArchivalSide &
                                oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_type'] == sArchivalType
                            ))
                            lRowSel <- min(which(
                              !is.na(oResCheck[,paste0(sSampleCheck,'_segment_chr')]) &
                                !is.na(oResCheck[,paste0(sSampleCheck,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))]) &
                                !is.na(oResCheck[,paste0(sSampleCheck,'_segment_',ifelse(sArchivalSide=='S','start','end'))]) &
                                !is.na(oResCheck[,paste0(sSampleCheck,'_segmented')]) &
                                !is.na(oResCheck[,paste0(sSampleCheck,'_Segment_Mean')]) &
                                !is.na(oResCheck[,paste0(sSampleCheck,'_Copies')]) &
                                oResCheck[,paste0(sSampleCheck,'_segment_chr')] == oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_chr')] &
                                oResCheck[,paste0(sSampleCheck,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))] == oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_cn_type_',ifelse(sArchivalSide=='S','start','end'))] &
                                oResCheck[,paste0(sSampleCheck,'_segment_',ifelse(sArchivalSide=='S','start','end'))] >= oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))] - (nShift * nAutosomesBinSize) &
                                oResCheck[,paste0(sSampleCheck,'_segment_',ifelse(sArchivalSide=='S','start','end'))] <= oResComparison[[paste0('S',nShift,sSampleType)]][nRowArch,paste0(sSampleArchival,'_segment_',ifelse(sArchivalSide=='S','start','end'))] + (nShift * nAutosomesBinSize)
                            ))
                            if (length(lRowSel)>0) {
                              oTmp <- rbind(
                                oTmp,
                                cbind(
                                  data.frame(
                                    archival_segment = sArchivalSegment,
                                    archival_segment_side = sArchivalSide,
                                    archival_segment_type = sArchivalType,
                                    stringsAsFactors = FALSE
                                  ),
                                  unique(oResCheck[lRowSel,paste(sSampleCheck,
                                                                 c(
                                                                   'segment','segment_chr','segment_start','segment_end',
                                                                   'segmented','Segment_Mean','Copies',
                                                                   'cn_type_start','cn_type_end'
                                                                 ),
                                                                 sep='_'
                                  )])
                                )
                              )
                            } else {
                              oTmp <- rbind(
                                oTmp,
                                oTmpEmpty[
                                  oTmpEmpty[,'archival_segment_side']==sArchivalSide &
                                    oTmpEmpty[,'archival_segment_type']==sArchivalType,
                                ]
                              )
                            }
                          }
                        }
                      }
                      if (!is.null(oTmp)) {
                        if (all(table(paste(oTmp[,'archival_segment'],oTmp[,'archival_segment_side'],oTmp[,'archival_segment_type'],sep='_'))==1)) {try({
                          oTmp <- merge(
                            oResComparison[[paste0('S',nShift,sSampleType)]],
                            oTmp,
                            all.x=TRUE)
                          print(paste(nShift,dim(oTmp)))
                          rownames(oTmp) <- paste(
                            oTmp[,'archival_segment'],
                            oTmp[,'archival_segment_side'],
                            oTmp[,'archival_segment_type'],
                            sep='_')
                          oTmp <- oTmp[order(
                            as.integer(gsub('X','23',gsub('Y','24',oTmp[,paste0(sSampleArchival,'_segment_chr')]))),
                            as.integer(oTmp[,paste0(sSampleArchival,'_segment_start')])
                          ),]
                          oResComparison[[paste0('S',nShift,sSampleType)]] <- oTmp
                          print(paste(nShift,dim(oResComparison[[paste0('S',nShift,sSampleType)]])))
                        })}
                      }
                      rm(oTmp)
                    }
                  }
                  
                  
                  ###################################################################################################
                  ### Remove unused objects
                  ###################################################################################################
                  rm(oResCheck)
                  gc();gc();gc();gc()
                }
              }
              
              
              ###################################################################################################
              ### Export comparisons
              ###################################################################################################
              for (sSampleType in c('p','c')) {
                for (nShift in lShifts) {
                  oResComparison[[paste0('S',nShift,sSampleType)]] <- oResComparison[[paste0('S',nShift,sSampleType)]][order(
                    as.integer(gsub('Y','24',gsub('X','23',oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]))),
                    oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_start')]),]
                  oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_summary_matched'] <- rowSums(!is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,
                                                                                                                                                                   grep('Copies',colnames(oResComparison[[paste0('S',nShift,sSampleType)]]),value=TRUE)
                  ]))-1
                  oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_summary_tested'] <- length(
                    grep('Copies',colnames(oResComparison[[paste0('S',nShift,sSampleType)]]),value=TRUE)
                  )-1
                  write.table(
                    oResComparison[[paste0('S',nShift,sSampleType)]],
                    file=file.path(
                      sPathCASCADE,
                      paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                      paste0(sPrefix,'__',sTrialID,'__',sSampleSel,'__',as.integer(nAutosomesBinSize),'__vs__',sSampleArchival,'__',paste0('S',nShift,sSampleType),'.txt')),
                    sep='\t',
                    quote=FALSE,
                    col.names=TRUE,
                    row.names=FALSE
                  )
                }
              }
              save(
                oResComparison,
                oSampleComarison,
                file=file.path(
                  sPathCASCADE,
                  paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                  paste0(sPrefix,'__',sTrialID,'__',sSampleSel,'__',as.integer(nAutosomesBinSize),'__vs__',sSampleArchival,'.RData'))
              )
              
              
              ###################################################################################################
              ### Create Archival report PDF
              ###################################################################################################
              pdf(file.path(
                sPathCASCADE,
                paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
                paste0(sPrefix,'__',sTrialID,'__',sSampleSel,'__',as.integer(nAutosomesBinSize),'__vs__',sSampleArchival,'__chr_details.pdf')),
                paper='special',
                width=16,
                height=12
              )
              
              
              ###################################################################################################
              ### Create summary report for archival comparisons
              ###################################################################################################
              for (nShift in lShifts) {
                ###################################################################################################
                ### Noise plot at different shift
                ###################################################################################################
                oDataPlot <- merge(
                  oResComparison[[paste0('S',nShift,'p')]][,c(
                    'archival_segment','archival_segment_side',
                    'archival_summary_matched','archival_summary_tested')],
                  oResComparison[[paste0('S',nShift,'c')]][,c(
                    'archival_segment','archival_segment_side',
                    'archival_summary_matched','archival_summary_tested')],
                  by=c('archival_segment','archival_segment_side')
                )
                colnames(oDataPlot) <- gsub('.x','_patient',colnames(oDataPlot),fixed=TRUE)
                colnames(oDataPlot) <- gsub('.y','_control',colnames(oDataPlot),fixed=TRUE)
                oDataPlot <- merge(
                  oDataPlot,
                  ddply(
                    oDataPlot,
                    c('archival_segment','archival_segment_side'),
                    function(o) data.frame(p.value = prop.test(
                      c(o[,'archival_summary_matched_patient'],o[,'archival_summary_matched_control']),
                      c(o[,'archival_summary_tested_patient'],o[,'archival_summary_tested_control']),
                      alternative = "two.sided",
                      correct = TRUE
                    )[['p.value']]
                    ))
                )
                oDataPlot[is.na(oDataPlot[,'p.value']),'p.value'] <- 1
                oDataPlot[,'Chromosome'] <- sapply(oDataPlot[,'archival_segment'],function(s) strsplit(s,':')[[1]][1])
                oPlot <- ggplot(oDataPlot)
                oPlot <- oPlot + geom_point(aes(
                  x = archival_summary_matched_patient / archival_summary_tested_patient,
                  y = archival_summary_matched_control / archival_summary_tested_control,
                  color=-log10(p.value)
                ),
                size = 8
                )
                oPlot <- oPlot + geom_text_repel(aes(
                  x = archival_summary_matched_patient / archival_summary_tested_patient,
                  y = archival_summary_matched_control / archival_summary_tested_control,
                  label = Chromosome
                ),
                size = 4
                )
                oPlot <- oPlot + theme_bw()
                oPlot <- oPlot + ggtitle(paste(
                  sTrialID,sPrefix,'matches with',
                  nShift,'Shift'))
                oPlot <- oPlot + theme(axis.text.x=element_text(angle=45, hjust=1))
                oPlot <- oPlot + scale_x_continuous(name = 'Fraction matched in patients mets')
                oPlot <- oPlot + scale_y_continuous(name = 'Fraction matched in controls')
                oPlot <- oPlot + scale_colour_gradient2(low='red',mid='yellow',high='green',midpoint=3)
                plot(oPlot)
                
                
                ###################################################################################################
                ### Match with tree at different shift
                ###################################################################################################
                par(mfrow=c(2,1))
                try(plot(as.dendrogram(lTrees[['tree_aut']])))
                try(plot(
                  colSums(!is.na(oResComparison[[paste0('S',nShift,'p')]][,paste0(lTrees[['tree_aut']][['order.lab']],'_Copies')])),
                  type='h',
                  xaxt = 'n',
                  xlab='Sample ID',
                  ylab='count',
                  ylim = c(
                    0,
                    max(colSums(!is.na(oResComparison[[paste0('S',nShift,'p')]][,paste0(lTrees[['tree_aut']][['order.lab']],'_Copies')])))
                  ),
                  main = paste(
                    'Comparison with',
                    sSampleArchival,sPrefix,
                    'at shift',nShift
                  )))
                try(text(
                  x = c(1:length(lTrees[['tree_aut']][['order.lab']])),
                  y = par("usr")[3] - 1,
                  labels = paste(
                    lTrees[['tree_aut']][['order.lab']],
                    lTrees[['tree_aut']][['order.lab.organ']],
                    sep='\n'),
                  xpd = NA,
                  srt = 35,
                  cex = 1.2))
                par(mfrow=c(2,1))
                
                
                ###################################################################################################
                ### Transition points heatmap at different shift
                ###################################################################################################
                try({
                  oDataPlot <- oResComparison[[paste0('S',nShift,'p')]][,
                                                                        grep('_segment$',colnames(oResComparison[[paste0('S',nShift,'p')]]),value=TRUE)
                  ]
                  colnames(oDataPlot) <- gsub('_segment','',colnames(oDataPlot),fixed=TRUE)
                  oDataPlot <- oDataPlot[,colnames(oDataPlot) %in% lTrees[['info']][,'Core_ID']]
                  oDataPlot <- as.data.frame(!is.na(oDataPlot))
                  oDataPlot %<>% mutate_if(is.logical,as.numeric) 
                  rownames(oDataPlot) <- rownames( oResComparison[[paste0('S',nShift,'p')]])
                  lHMcolOrgan <- unique(oDataStats[colnames(oDataPlot),c('Organ','Organ.Color')])[,'Organ.Color']
                  names(lHMcolOrgan) <- unique(oDataStats[colnames(oDataPlot),c('Organ','Organ.Color')])[,'Organ']
                  lHMcolAR <- colorRamp2(c(0,2,3,5,10), c("green","green", "yellow", "orange","red"))
                  draw(Heatmap(
                    as.matrix(oDataPlot),
                    row_title = 'Transition Points',
                    name = paste('Transition Points\nat Shift',nShift),
                    cluster_rows = TRUE,
                    cluster_columns = as.dendrogram(lTrees[['tree_aut']]),
                    col = c('0'='white','1'='black'),
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    row_split = sapply(rownames(oDataPlot),function(s) strsplit(s,':')[[1]][1]),
                    column_names_gp = gpar(fontsize = 4),
                    top_annotation = HeatmapAnnotation(
                      Organ = lTrees[['info']][colnames(oDataPlot),'Organ'],
                      ARcnLAB = oDataStats[colnames(oDataPlot),'ARcnLAB'],
                      foo = anno_block(gp = gpar(fill = c('white','gray')),
                                       labels_gp = gpar(col = "black", fontsize = 16)),
                      col = list(
                        Organ = lHMcolOrgan,
                        ARcnLAB = lHMcolAR))
                  ) 
                  )
                })
              }
              
              
              ###################################################################################################
              ### Create detailed report for archival samples
              ###################################################################################################
              for (sSampleCurrent in names(oSampleComarison)) {
                for (sChr in unique(oResComparison[[paste0('S',max(lShifts),'p')]][,
                                                                                   paste0(sSampleArchival,'_segment_chr')])[order(as.integer(gsub('Y','24',gsub('X','23',
                                                                                                                                                                unique(oResComparison[[paste0('S',max(lShifts),'p')]][,paste0(sSampleArchival,'_segment_chr')])
                                                                                   ))))]) {try({
                                                                                     print(paste0(sSampleCurrent,' chr',sChr))
                                                                                     
                                                                                     
                                                                                     ###################################################################################################
                                                                                     ### Check is the current sample has valid CN segments
                                                                                     ###################################################################################################
                                                                                     if (paste0(sSampleCurrent,'_Copies') %in% colnames(oResComparison[[paste0('S',max(lShifts),'p')]])) {
                                                                                       sSampleType <- 'p'
                                                                                       
                                                                                       
                                                                                       ###################################################################################################
                                                                                       ### Create plot by Copies
                                                                                       ###################################################################################################
                                                                                       oPlotC <- ggplot()
                                                                                       oPlotC <- oPlotC + geom_segment(
                                                                                         data = oRes[
                                                                                           !is.na(oRes[,paste0(sSampleArchival,'_segment_chr')]) &
                                                                                             oRes[,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           c('Start','End',paste0(sSampleArchival,'_Copies'))
                                                                                         ],
                                                                                         aes(
                                                                                           x = Start,
                                                                                           xend = End,
                                                                                           y = get(paste0(sSampleArchival,'_Copies')) + 0.02,
                                                                                           yend = get(paste0(sSampleArchival,'_Copies')) + 0.02
                                                                                         ),
                                                                                         color = 'red'
                                                                                       )
                                                                                       oPlotC <- oPlotC + geom_segment(
                                                                                         data = oRes[
                                                                                           !is.na(oRes[,paste0(sSampleCurrent,'_segment_chr')]) &
                                                                                             oRes[,paste0(sSampleCurrent,'_segment_chr')]==sChr,
                                                                                           c('Start','End',paste0(sSampleCurrent,'_Copies'))
                                                                                         ],
                                                                                         aes(
                                                                                           x = Start,
                                                                                           xend = End,
                                                                                           y = get(paste0(sSampleCurrent,'_Copies')) - 0.02,
                                                                                           yend = get(paste0(sSampleCurrent,'_Copies')) - 0.02 
                                                                                         ),
                                                                                         color = 'blue'
                                                                                       )
                                                                                       oPlotC <- oPlotC + geom_point(
                                                                                         data = oResComparison[[paste0('S',max(lShifts),'p')]][
                                                                                           oResComparison[[paste0('S',max(lShifts),'p')]][,'archival_segment_side'] == 'S' &
                                                                                             oResComparison[[paste0('S',max(lShifts),'p')]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           paste(sSampleArchival,c('segment_start','Copies'),sep='_')
                                                                                         ],
                                                                                         aes(
                                                                                           x = get(paste0(sSampleArchival,'_segment_start')),
                                                                                           y = get(paste0(sSampleArchival,'_Copies')) + 0.02
                                                                                         ),
                                                                                         size=2,
                                                                                         color = 'black'
                                                                                       )
                                                                                       for (nShift in sort(lShifts,decreasing=TRUE)) {
                                                                                         oPlotC <- oPlotC + geom_segment(
                                                                                           data = oResComparison[[paste0('S',nShift,sSampleType)]][
                                                                                             oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'] == 'S' &
                                                                                               !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_Copies')]) &
                                                                                               oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                             c(
                                                                                               paste(sSampleArchival,c('segment_start','Copies'),sep='_'),
                                                                                               paste0(sSampleCurrent,'_Copies')
                                                                                             )
                                                                                           ],
                                                                                           aes(
                                                                                             x = get(paste0(sSampleArchival,'_segment_start')),
                                                                                             xend = get(paste0(sSampleArchival,'_segment_start')),
                                                                                             y = ifelse(
                                                                                               get(paste0(sSampleArchival,'_Copies')) > get(paste0(sSampleCurrent,'_Copies')),
                                                                                               get(paste0(sSampleArchival,'_Copies')) - 0.1,
                                                                                               get(paste0(sSampleArchival,'_Copies')) + 0.1
                                                                                             ),
                                                                                             yend = ifelse(
                                                                                               get(paste0(sSampleArchival,'_Copies')) > get(paste0(sSampleCurrent,'_Copies')),
                                                                                               get(paste0(sSampleCurrent,'_Copies')) + 0.1,
                                                                                               get(paste0(sSampleCurrent,'_Copies')) - 0.1
                                                                                             )
                                                                                           ),
                                                                                           arrow=arrow(),
                                                                                           size=2,
                                                                                           color = ifelse(nShift==0,'green',ifelse(nShift==1,'orange','red'))
                                                                                         )
                                                                                       }
                                                                                       oPlotC <- oPlotC + geom_point(
                                                                                         data = oResComparison[[paste0('S',max(lShifts),'p')]][
                                                                                           oResComparison[[paste0('S',max(lShifts),'p')]][,'archival_segment_side'] == 'E' &
                                                                                             oResComparison[[paste0('S',max(lShifts),'p')]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           paste(sSampleArchival,c('segment_end','Copies'),sep='_')
                                                                                         ],
                                                                                         aes(
                                                                                           x = get(paste0(sSampleArchival,'_segment_end')),
                                                                                           y = get(paste0(sSampleArchival,'_Copies')) + 0.02
                                                                                         ),
                                                                                         size=2,
                                                                                         color = 'black'
                                                                                       )
                                                                                       for (nShift in sort(lShifts,decreasing=TRUE)) {
                                                                                         oPlotC <- oPlotC + geom_segment(
                                                                                           data = oResComparison[[paste0('S',nShift,sSampleType)]][
                                                                                             oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'] == 'E' &
                                                                                               !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_Copies')]) &
                                                                                               oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                             c(
                                                                                               paste(sSampleArchival,c('segment_end','Copies'),sep='_'),
                                                                                               paste0(sSampleCurrent,'_Copies')
                                                                                             )
                                                                                           ],
                                                                                           aes(
                                                                                             x = get(paste0(sSampleArchival,'_segment_end')),
                                                                                             xend = get(paste0(sSampleArchival,'_segment_end')),
                                                                                             y = ifelse(
                                                                                               get(paste0(sSampleArchival,'_Copies')) > get(paste0(sSampleCurrent,'_Copies')),
                                                                                               get(paste0(sSampleArchival,'_Copies')) - 0.1,
                                                                                               get(paste0(sSampleArchival,'_Copies')) + 0.1
                                                                                             ),
                                                                                             yend = ifelse(
                                                                                               get(paste0(sSampleArchival,'_Copies')) > get(paste0(sSampleCurrent,'_Copies')),
                                                                                               get(paste0(sSampleCurrent,'_Copies')) + 0.1,
                                                                                               get(paste0(sSampleCurrent,'_Copies')) - 0.1
                                                                                             )
                                                                                           ),
                                                                                           arrow=arrow(),
                                                                                           size=2,
                                                                                           color = ifelse(nShift==0,'green',ifelse(nShift==1,'orange','red'))
                                                                                         )
                                                                                       }
                                                                                       oPlotC <- oPlotC + theme_bw()
                                                                                       oPlotC <- oPlotC + ggtitle(paste(
                                                                                         sSampleArchival,'(',sPrefix,') vs.',
                                                                                         sSampleCurrent,'by ACE Copies'))
                                                                                       oPlotC <- oPlotC + theme(axis.text.x=element_text(angle=45, hjust=1))
                                                                                       oPlotC <- oPlotC + scale_x_continuous(
                                                                                         name = paste0('Position in chr',sChr),
                                                                                         limits = c(0,max(oRes[oRes[,'Chromosome']==sChr,'End'],na.rm=TRUE)),
                                                                                         n.breaks = 20,
                                                                                         labels=function(n) format(n,scientific=FALSE,digits=0)
                                                                                       )
                                                                                       oPlotC <- oPlotC + scale_y_continuous(
                                                                                         name = 'ACE Copies',
                                                                                         limits = c(-1,max(oRes[,paste0(c(sSampleArchival,sSampleCurrent),'_Copies')],na.rm=TRUE))
                                                                                       )
                                                                                       
                                                                                       
                                                                                       ###################################################################################################
                                                                                       ### Create SV plot 
                                                                                       ###################################################################################################
                                                                                       oPlotSV <- ggplot(oResSV[
                                                                                         oResSV[,'Trial.ID']==sTrialID &
                                                                                           oResSV[,'chr_start']==sChr,
                                                                                       ])
                                                                                       oPlotSV <- oPlotSV + geom_rect(
                                                                                         aes(
                                                                                           xmin = floor(pos_start/nAutosomesBinSize)*nAutosomesBinSize,
                                                                                           xmax = ceiling(pos_start/nAutosomesBinSize)*nAutosomesBinSize,
                                                                                           ymin = c('DUP'=4,'DEL'=1,'INV'=2,'BND'=3)[type] - 0.5,
                                                                                           ymax = c('DUP'=4,'DEL'=1,'INV'=2,'BND'=3)[type] + 0.5,
                                                                                           fill = type
                                                                                         )
                                                                                       )
                                                                                       oPlotSV <- oPlotSV + scale_fill_manual(
                                                                                         values = c('DUP' = 'red','DEL' = 'blue','INV' = 'green','BND' = 'green'),
                                                                                         guide=FALSE
                                                                                       )
                                                                                       oPlotSV <- oPlotSV + theme_bw()
                                                                                       oPlotSV <- oPlotSV + theme(
                                                                                         axis.text.x=element_text(angle=45, hjust=1))
                                                                                       oPlotSV <- oPlotSV + scale_x_continuous(
                                                                                         name = paste0('Position in chr',sChr),
                                                                                         limits = c(0,max(oRes[oRes[,'Chromosome']==sChr,'End'],na.rm=TRUE)),
                                                                                         n.breaks = 20,
                                                                                         labels=function(n) format(n,scientific=FALSE,digits=0)
                                                                                       )
                                                                                       oPlotSV <- oPlotSV + scale_y_continuous(
                                                                                         name = NULL,
                                                                                         breaks  = 1:4,
                                                                                         labels = c('DEL','INV','BND','DUP')
                                                                                       )
                                                                                       
                                                                                       
                                                                                       ###################################################################################################
                                                                                       ### Create plot by segment means
                                                                                       ###################################################################################################
                                                                                       oPlotM <- ggplot()
                                                                                       oPlotM <- oPlotM + geom_segment(
                                                                                         data = oRes[
                                                                                           !is.na(oRes[,paste0(sSampleArchival,'_segment_chr')]) &
                                                                                             oRes[,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           c('Start','End',paste0(sSampleArchival,'_segmented'))
                                                                                         ],
                                                                                         aes(
                                                                                           x = Start,
                                                                                           xend = End,
                                                                                           y = get(paste0(sSampleArchival,'_segmented')),
                                                                                           yend = get(paste0(sSampleArchival,'_segmented'))
                                                                                         ),
                                                                                         color = 'red'
                                                                                       )
                                                                                       oPlotM <- oPlotM + geom_segment(
                                                                                         data = oRes[
                                                                                           !is.na(oRes[,paste0(sSampleArchival,'_segment_chr')]) &
                                                                                             oRes[,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           c('Start','End',paste0(sSampleArchival,'_copynumber'))
                                                                                         ],
                                                                                         aes(
                                                                                           x = Start,
                                                                                           xend = End,
                                                                                           y = get(paste0(sSampleArchival,'_copynumber')),
                                                                                           yend = get(paste0(sSampleArchival,'_copynumber'))
                                                                                         ),
                                                                                         color = 'red'
                                                                                       )
                                                                                       oPlotM <- oPlotM + geom_segment(
                                                                                         data = oRes[
                                                                                           !is.na(oRes[,paste0(sSampleCurrent,'_segment_chr')]) &
                                                                                             oRes[,paste0(sSampleCurrent,'_segment_chr')]==sChr,
                                                                                           c('Start','End',paste0(sSampleCurrent,'_segmented'))
                                                                                         ],
                                                                                         aes(
                                                                                           x = Start,
                                                                                           xend = End,
                                                                                           y = get(paste0(sSampleCurrent,'_segmented')),
                                                                                           yend = get(paste0(sSampleCurrent,'_segmented')) 
                                                                                         ),
                                                                                         color = 'blue'
                                                                                       )
                                                                                       oPlotM <- oPlotM + geom_segment(
                                                                                         data = oRes[
                                                                                           !is.na(oRes[,paste0(sSampleCurrent,'_segment_chr')]) &
                                                                                             oRes[,paste0(sSampleCurrent,'_segment_chr')]==sChr,
                                                                                           c('Start','End',paste0(sSampleCurrent,'_copynumber'))
                                                                                         ],
                                                                                         aes(
                                                                                           x = Start,
                                                                                           xend = End,
                                                                                           y = get(paste0(sSampleCurrent,'_copynumber')),
                                                                                           yend = get(paste0(sSampleCurrent,'_copynumber')) 
                                                                                         ),
                                                                                         color = 'blue'
                                                                                       )
                                                                                       oPlotM <- oPlotM + geom_point(
                                                                                         data = oResComparison[[paste0('S',max(lShifts),'p')]][
                                                                                           oResComparison[[paste0('S',max(lShifts),'p')]][,'archival_segment_side'] == 'S' &
                                                                                             oResComparison[[paste0('S',max(lShifts),'p')]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           paste(sSampleArchival,c('segment_start','segmented'),sep='_')
                                                                                         ],
                                                                                         aes(
                                                                                           x = get(paste0(sSampleArchival,'_segment_start')),
                                                                                           y = get(paste0(sSampleArchival,'_segmented'))
                                                                                         ),
                                                                                         size=2,
                                                                                         color = 'black'
                                                                                       )
                                                                                       for (nShift in sort(lShifts,decreasing=TRUE)) {
                                                                                         oPlotM <- oPlotM + geom_segment(
                                                                                           data = oResComparison[[paste0('S',nShift,sSampleType)]][
                                                                                             oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'] == 'S' &
                                                                                               !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segmented')]) &
                                                                                               oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                             c(
                                                                                               paste(sSampleArchival,c('segment_start','segmented'),sep='_'),
                                                                                               paste0(sSampleCurrent,'_segmented')
                                                                                             )
                                                                                           ],
                                                                                           aes(
                                                                                             x = get(paste0(sSampleArchival,'_segment_start')),
                                                                                             xend = get(paste0(sSampleArchival,'_segment_start')),
                                                                                             y = ifelse(
                                                                                               get(paste0(sSampleArchival,'_segmented')) > get(paste0(sSampleCurrent,'_segmented')),
                                                                                               get(paste0(sSampleArchival,'_segmented')) - 0.1,
                                                                                               get(paste0(sSampleArchival,'_segmented')) + 0.1
                                                                                             ),
                                                                                             yend = ifelse(
                                                                                               get(paste0(sSampleArchival,'_segmented')) > get(paste0(sSampleCurrent,'_segmented')),
                                                                                               get(paste0(sSampleCurrent,'_segmented')) + 0.1,
                                                                                               get(paste0(sSampleCurrent,'_segmented')) - 0.1
                                                                                             )
                                                                                           ),
                                                                                           arrow=arrow(),
                                                                                           size=2,
                                                                                           color = ifelse(nShift==0,'green',ifelse(nShift==1,'orange','red'))
                                                                                         )
                                                                                       }
                                                                                       oPlotM <- oPlotM + geom_point(
                                                                                         data = oResComparison[[paste0('S',max(lShifts),'p')]][
                                                                                           oResComparison[[paste0('S',max(lShifts),'p')]][,'archival_segment_side'] == 'E' &
                                                                                             oResComparison[[paste0('S',max(lShifts),'p')]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                           paste(sSampleArchival,c('segment_end','segmented'),sep='_')
                                                                                         ],
                                                                                         aes(
                                                                                           x = get(paste0(sSampleArchival,'_segment_end')),
                                                                                           y = get(paste0(sSampleArchival,'_segmented'))
                                                                                         ),
                                                                                         size=2,
                                                                                         color = 'black'
                                                                                       )
                                                                                       for (nShift in sort(lShifts,decreasing=TRUE)) {
                                                                                         oPlotM <- oPlotM + geom_segment(
                                                                                           data = oResComparison[[paste0('S',nShift,sSampleType)]][
                                                                                             oResComparison[[paste0('S',nShift,sSampleType)]][,'archival_segment_side'] == 'E' &
                                                                                               !is.na(oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleCurrent,'_segmented')]) &
                                                                                               oResComparison[[paste0('S',nShift,sSampleType)]][,paste0(sSampleArchival,'_segment_chr')]==sChr,
                                                                                             c(
                                                                                               paste(sSampleArchival,c('segment_end','segmented'),sep='_'),
                                                                                               paste0(sSampleCurrent,'_segmented')
                                                                                             )
                                                                                           ],
                                                                                           aes(
                                                                                             x = get(paste0(sSampleArchival,'_segment_end')),
                                                                                             xend = get(paste0(sSampleArchival,'_segment_end')),
                                                                                             y = ifelse(
                                                                                               get(paste0(sSampleArchival,'_segmented')) > get(paste0(sSampleCurrent,'_segmented')),
                                                                                               get(paste0(sSampleArchival,'_segmented')) - 0.1,
                                                                                               get(paste0(sSampleArchival,'_segmented')) + 0.1
                                                                                             ),
                                                                                             yend = ifelse(
                                                                                               get(paste0(sSampleArchival,'_segmented')) > get(paste0(sSampleCurrent,'_segmented')),
                                                                                               get(paste0(sSampleCurrent,'_segmented')) + 0.1,
                                                                                               get(paste0(sSampleCurrent,'_segmented')) - 0.1
                                                                                             )
                                                                                           ),
                                                                                           arrow=arrow(),
                                                                                           size=2,
                                                                                           color = ifelse(nShift==0,'green',ifelse(nShift==1,'orange','red'))
                                                                                         )
                                                                                       }
                                                                                       oPlotM <- oPlotM + theme_bw()
                                                                                       oPlotM <- oPlotM + ggtitle(paste(
                                                                                         sSampleArchival,'(',sPrefix,') vs.',
                                                                                         sSampleCurrent,'by QDNAseq Segment Mean'))
                                                                                       oPlotM <- oPlotM + theme(axis.text.x=element_text(angle=45, hjust=1))
                                                                                       oPlotM <- oPlotM + scale_x_continuous(
                                                                                         name = paste0('Position in chr',sChr),
                                                                                         limits = c(0,max(oRes[oRes[,'Chromosome']==sChr,'End'],na.rm=TRUE)),
                                                                                         n.breaks = 20,
                                                                                         labels=function(n) format(n,scientific=FALSE,digits=0)
                                                                                       )
                                                                                       oPlotM <- oPlotM + scale_y_continuous(
                                                                                         name = 'QDNAseq Segment_Mean',
                                                                                         limits = c(
                                                                                           min(oRes[,paste0(c(sSampleArchival,sSampleCurrent),'_copynumber')],na.rm=TRUE),
                                                                                           max(oRes[,paste0(c(sSampleArchival,sSampleCurrent),'_copynumber')],na.rm=TRUE)
                                                                                         )
                                                                                       )
                                                                                     }
                                                                                     
                                                                                     
                                                                                     ###################################################################################################
                                                                                     ### Subset ACE plot
                                                                                     ###################################################################################################
                                                                                     lTmpBckup <- list()
                                                                                     oPlotA <- oSampleComarison[[sSampleCurrent]][['compareplot']]
                                                                                     for (nLayer in c(1:length(oPlotA[['layers']]))) {
                                                                                       lTmpBckup[[nLayer]] <- duplicate(oPlotA[['layers']][[nLayer]][['data']], shallow = FALSE)
                                                                                       if (any(names(oPlotA[['layers']][[nLayer]][['data']])=='chr')) {
                                                                                         oPlotA[['layers']][[nLayer]][['data']] <- oPlotA[['layers']][[nLayer]][['data']][oPlotA[['layers']][[nLayer]][['data']][,'chr']==sChr,]
                                                                                       }
                                                                                     }
                                                                                     
                                                                                     
                                                                                     ###################################################################################################
                                                                                     ### Plot both graph
                                                                                     ###################################################################################################
                                                                                     if (paste0(sSampleCurrent,'_Copies') %in% colnames(oResComparison[[paste0('S',max(lShifts),'p')]])) {
                                                                                       grid.arrange(
                                                                                         oPlotA,oPlotC,oPlotSV,oPlotM,
                                                                                         layout_matrix = rbind(
                                                                                           c(1),c(1),
                                                                                           c(2),c(2),
                                                                                           c(3),
                                                                                           c(4),c(4)
                                                                                         ))
                                                                                     } else {
                                                                                       plot(oPlotA)
                                                                                     }
                                                                                     
                                                                                     
                                                                                     ###################################################################################################
                                                                                     ### Restore ACE plot data
                                                                                     ###################################################################################################
                                                                                     for (nLayer in c(1:length(oPlotA[['layers']]))) {
                                                                                       oPlotA[['layers']][[nLayer]][['data']] <- duplicate(lTmpBckup[[nLayer]], shallow = FALSE)
                                                                                     }
                                                                                   })}
              }
              
              
              ###################################################################################################
              ### Close PDF
              ###################################################################################################
              try(dev.off());try(dev.off());try(dev.off());try(dev.off())
              
              
              ###################################################################################################
              ### Remove Archival comparisons data
              ###################################################################################################
              rm(oResComparison,oSampleComarison)
              gc();gc();gc();gc()
            }
          }
          
          
          ###################################################################################################
          ### Remove patient analysis objects
          ###################################################################################################
          rm(oQM,oSM,oRes,lTrees)
          gc();gc();gc();gc()
        }
      })}
      
      
      ###################################################################################################
      ### Export full distances table
      ###################################################################################################
      write.table(
        oDataHCdist[,
                    colnames(oDataHCdist)[
                      regexpr('dist_',colnames(oDataHCdist),fixed=TRUE) == -1 | (
                        regexpr('dist_',colnames(oDataHCdist),fixed=TRUE) >-1 &
                          regexpr(sSampleSel,colnames(oDataHCdist),fixed=TRUE) >-1 &
                          regexpr(as.integer(nAutosomesBinSize),colnames(oDataHCdist),fixed=TRUE) >-1
                      )]
        ],
        file=file.path(
          sPathCASCADE,
          paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
          paste0('Distances__',sSampleSel,'__A',as.integer(nAutosomesBinSize),'.txt')),
        sep='\t',
        quote=FALSE,
        col.names=TRUE,
        row.names=FALSE
      )
      
      
      ###################################################################################################
      ### Compare different distances
      ###################################################################################################
      pdf(file.path(
        sPathCASCADE,
        paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
        paste0('Distances__',sSampleSel,'__A',as.integer(nAutosomesBinSize),'.pdf')),
        paper='special',
        width=16,
        height=12)
      for (sColDist in colnames(oDataHCdist)[
        regexpr('dist_',colnames(oDataHCdist),fixed=TRUE) >-1 &
        regexpr(sSampleSel,colnames(oDataHCdist),fixed=TRUE) >-1 &
        regexpr(as.integer(nAutosomesBinSize),colnames(oDataHCdist),fixed=TRUE) >-1
      ]) {
        ###################################################################################################
        ### Plot by Cluster
        ###################################################################################################
        oPlot <- ggplot(oDataHCdist)
        oPlot <- oPlot + geom_point(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist),
          fill = within_same_cluster),
          shape = 21,
          size = 4, 
          color = 'black'
        )
        oPlot <- oPlot + geom_smooth(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist)),
          method=lm,
          color = 'blue'
        )
        oPlot <- oPlot + stat_cor(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist)),
          method = "pearson", 
          label.x = max(oDataHCdist[,'percent_common_nonsilent_mutations'],na.rm=TRUE)*0.8, 
          label.y = max(oDataHCdist[,sColDist],na.rm=TRUE)*0.8)
        oPlot <- oPlot + theme_bw()
        oPlot <- oPlot + xlab('% of common non-silent mutations')
        oPlot <- oPlot + ylab(sColDist)
        plot(oPlot)
        
        
        ###################################################################################################
        ### Plot by Patient
        ###################################################################################################
        oPlot <- ggplot(oDataHCdist)
        oPlot <- oPlot + geom_point(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist),
          fill = oDataStats[Core_ID_1,'Trial.ID']),
          shape = 21,
          size = 4, 
          color = 'black'
        )
        oPlot <- oPlot + geom_text_repel(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist),
          color = oDataStats[Core_ID_1,'Trial.ID'],
          label = paste0(
            oDataStats[Core_ID_1,'Figure_N'],'∩',
            oDataStats[Core_ID_2,'Figure_N'])),
          shape = 21,
          size = 4, 
          color = 'black'
        )
        oPlot <- oPlot + geom_smooth(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist)),
          method=lm,
          color = 'blue'
        )
        oPlot <- oPlot + stat_cor(aes(
          x = percent_common_nonsilent_mutations,
          y = get(sColDist)),
          method = "pearson", 
          label.x = max(oDataHCdist[,'percent_common_nonsilent_mutations'],na.rm=TRUE)*0.8, 
          label.y = max(oDataHCdist[,sColDist],na.rm=TRUE)*0.8)
        oPlot <- oPlot + theme_bw()
        oPlot <- oPlot + xlab('% of common non-silent mutations')
        oPlot <- oPlot + ylab(sColDist)
        plot(oPlot)
      }
      try(dev.off());try(dev.off());try(dev.off());try(dev.off())
      
      
      ###################################################################################################
      ### Save complete HM objects
      ###################################################################################################
      try(save(
        lSamplesOrderCluster,lDataPlotChrX,lDataPlotChrXregionCol,lInfoSamples,
        file=file.path(
          sPathCASCADE,
          paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
          paste0('Heatmaps__',sSampleSel,'__A',as.integer(nAutosomesBinSize),'.RData')
        )
      ))
      
      
      ###################################################################################################
      ### Plot chrX full heatmap
      ###################################################################################################
      for (sHMplotName in names(lInfoSamples)) {
        ###################################################################################################
        ### Export sample list
        ###################################################################################################
        write.table(
          lInfoSamples[[sHMplotName]],
          file=file.path(
            sPathCASCADE,
            paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
            paste0('Heatmaps__',sHMplotName ,'__sample_list.txt')),
          sep='\t',
          quote=FALSE,
          col.names=TRUE,
          row.names=FALSE
        )
        
        
        ###################################################################################################
        ### If the data correspond to a heatmap data table
        ###################################################################################################
        if (sHMplotName %in% names(lDataPlotChrX)) {try({
          ###################################################################################################
          ### Save HM data objects
          ###################################################################################################
          oSamplesOrderCluster <- lSamplesOrderCluster[[sHMplotName]]
          oDataPlotChrX <- lDataPlotChrX[[sHMplotName]][,names(oSamplesOrderCluster)]
          oDataPlotChrXregionCol <- lDataPlotChrXregionCol[[sHMplotName]]
          save(
            oDataPlotChrX,oDataPlotChrXregionCol,oSamplesOrderCluster,
            file=file.path(
              sPathCASCADE,
              paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
              paste0('Heatmaps__',sHMplotName ,'.RData'))
          )
          lHMcolPatient <- c('gray','white')[(1:length(unique(oDataStats[names(oSamplesOrderCluster),'Trial.ID'])) %% 2) + 1]
          names(lHMcolPatient) <- unique(oDataStats[names(oSamplesOrderCluster),'Trial.ID'])
          
          
          ###################################################################################################
          ### Create PDF
          ###################################################################################################
          pdf(file.path(
            sPathCASCADE,
            paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
            paste0('Heatmaps__',sHMplotName ,'.pdf')),
            paper='special',
            width=16,
            height=12)
          
          
          ###################################################################################################
          ### Aggreate chrX regions for the heatmap
          ###################################################################################################
          lChrXRegionsGroupKB <- paste0(
            sapply(rownames(oDataPlotChrX),function(s) strsplit(s,':')[[1]][1]),':',
            formatC(as.integer((floor(as.integer(sapply(rownames(oDataPlotChrX),function(s) strsplit(strsplit(s,':')[[1]][2],'-')[[1]][1]))/nHMrowgroups)*nHMrowgroups)+1),width=10,flag='0'),'-',
            formatC(as.integer((floor((as.integer(sapply(rownames(oDataPlotChrX),function(s) strsplit(strsplit(s,':')[[1]][2],'-')[[1]][2]))/nHMrowgroups)+1)*nHMrowgroups)),width=10,flag='0')
          )
          names(lChrXRegionsGroupKB) <- rownames(oDataPlotChrX)
          lDataPlotAnno <- oDataPlotChrXregionCol[,'Region']
          names(lDataPlotAnno) <- rownames(oDataPlotChrXregionCol)
          lHMcolOrgan <- unique(oDataStats[colnames(oDataPlotChrX),c('Organ','Organ.Color')])[,'Organ.Color']
          names(lHMcolOrgan) <- unique(oDataStats[colnames(oDataPlotChrX),c('Organ','Organ.Color')])[,'Organ']
          lColsRows <- list(Region=c(other='gray',AR='red',ARenh='blue',centromere='black'))
          for (sSVcol in colnames(oDataPlotChrXregionCol)[colnames(oDataPlotChrXregionCol)!='Region']) {
            lColsRows <- c(lColsRows,list(tmp=c(SVs='black'))) 
            names(lColsRows)[names(lColsRows)=='tmp'] <- sSVcol
          }
          lHMcolAR <- colorRamp2(c(0,2,3,5,10), c("green","green", "yellow", "orange","red"))
          
          
          ###################################################################################################
          ### Draw heatmaps
          ###################################################################################################
          oDataPlotChrX <- oDataPlotChrX[order(sapply(rownames(oDataPlotChrX),function(s) as.numeric(strsplit(s,'-')[[1]][2]))),]
          oDataPlotChrXregionCol <- oDataPlotChrXregionCol[rownames(oDataPlotChrX),]
          for (oRowSplit in list(
            a=NULL
            #b=lChrXRegionsGroupKB[rownames(oDataPlotChrX)]
          )) {
            try(draw(Heatmap(
              as.matrix(oDataPlotChrX[,names(oSamplesOrderCluster)]),
              row_title = paste('CN profile of',sHMplotName ),
              name = paste0('CN',sMeasureType),
              column_split = oDataStats[names(oSamplesOrderCluster),'Trial.ID'],
              row_split = oRowSplit,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              clustering_method_rows = 'ward',
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 4),
              col = lHMcolAR,
              top_annotation = HeatmapAnnotation(
                Patient = oDataStats[names(oSamplesOrderCluster),'Trial.ID'],
                Organ = oDataStats[names(oSamplesOrderCluster),'Organ'],
                ARcnLAB = oDataStats[names(oSamplesOrderCluster),'ARcnLAB'],
                ACEtc = ifelse(
                  oDataStats[names(oSamplesOrderCluster),'ACE_cellularity_all']==1,
                  0,
                  oDataStats[names(oSamplesOrderCluster),'ACE_cellularity_all']
                ),
                foo = anno_block(gp = gpar(fill = c('white','gray')),
                                 labels_gp = gpar(col = "black", fontsize = 16)),
                col = list(
                  Patient = lHMcolPatient,
                  Organ = lHMcolOrgan,
                  ARcnLAB = lHMcolAR,
                  ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green')))),
              left_annotation = rowAnnotation(
                df=oDataPlotChrXregionCol[rownames(oDataPlotChrX),],
                col=lColsRows
              )
            )))
            
            ###################################################################################################
            ### Draw extra heatmaps based on extracter table data
            ###################################################################################################
            lHMextras <- list.files(file.path(sPathCASCADE,'details_heatmaps'),pattern='Heatmaps__',full.name=TRUE)
            for (sHMextra in lHMextras) {
              ###################################################################################################
              ### Read extra data table
              ###################################################################################################
              oHMextra <- read.delim(
                sHMextra,
                header=TRUE,
                row.names = 1,
                sep='\t',
                stringsAsFactors=FALSE)
              
              
              ###################################################################################################
              ### Add annotation for the region of interest
              ###################################################################################################
              oHMextraAnn <- data.frame(
                region = rownames(oHMextra),
                Chromosome = sapply(rownames(oHMextra),function(s) strsplit(s,':')[[1]][1]),
                Start = sapply(rownames(oHMextra),function(s) as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][1])),
                End = sapply(rownames(oHMextra),function(s) as.integer(strsplit(s,'-')[[1]][2])),
                Region = 'other',
                stringsAsFactors = FALSE
              )
              for (nRegion in 1:nrow(oRegions)) {
                oHMextraAnn[
                  oHMextraAnn[,'Chromosome'] == oRegions[nRegion,'Chromosome'] & (
                    (oHMextraAnn[,'Start']<=oRegions[nRegion,'Start'] & oHMextraAnn[,'End']>=oRegions[nRegion,'End']) |
                      (oHMextraAnn[,'Start']>=oRegions[nRegion,'Start'] & oHMextraAnn[,'End']<=oRegions[nRegion,'End']) |
                      (oHMextraAnn[,'Start']>=oRegions[nRegion,'Start'] & oHMextraAnn[,'Start']<=oRegions[nRegion,'End']) |
                      (oHMextraAnn[,'End']>=oRegions[nRegion,'Start'] & oHMextraAnn[,'End']<=oRegions[nRegion,'End'])
                  ),
                  'Region'] <- oRegions[nRegion,'ID']
              }
              
              
              
              ###################################################################################################
              ### Add annotation for the region of interest
              ###################################################################################################
              try(draw(Heatmap(
                as.matrix(oHMextra[,names(oSamplesOrderCluster)[
                  names(oSamplesOrderCluster) %in% colnames(oHMextra)
                ]]),
                row_title = paste('CN profile of',gsub('.txt','',basename(sHMextra),fixed=TRUE) ),
                name = paste0('CN',sMeasureType),
                column_split = oDataStats[names(oSamplesOrderCluster)[
                  names(oSamplesOrderCluster) %in% colnames(oHMextra)
                ],'Trial.ID'],
                row_split = oRowSplit,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                clustering_method_rows = 'ward',
                show_row_names = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 4),
                col = lHMcolAR,
                top_annotation = HeatmapAnnotation(
                  Patient = oDataStats[names(oSamplesOrderCluster)[
                    names(oSamplesOrderCluster) %in% colnames(oHMextra)
                  ],'Trial.ID'],
                  Organ = oDataStats[names(oSamplesOrderCluster)[
                    names(oSamplesOrderCluster) %in% colnames(oHMextra)
                  ],'Organ'],
                  ARcnLAB = oDataStats[names(oSamplesOrderCluster)[
                    names(oSamplesOrderCluster) %in% colnames(oHMextra)
                  ],'ARcnLAB'],
                  ACEtc = ifelse(
                    oDataStats[names(oSamplesOrderCluster)[
                      names(oSamplesOrderCluster) %in% colnames(oHMextra)
                    ],'ACE_cellularity_all']==1,
                    0,
                    oDataStats[names(oSamplesOrderCluster)[
                      names(oSamplesOrderCluster) %in% colnames(oHMextra)
                    ],'ACE_cellularity_all']
                  ),
                  foo = anno_block(gp = gpar(fill = c('white','gray')),
                                   labels_gp = gpar(col = "black", fontsize = 16)),
                  col = list(
                    Patient = lHMcolPatient,
                    Organ = lHMcolOrgan,
                    ARcnLAB = lHMcolAR,
                    ACEtc = colorRamp2(c(0,0.2,0.4), c('red', 'orange', 'green')))),
                left_annotation = rowAnnotation(
                  df=oHMextraAnn[,'Region',drop=FALSE],
                  col=lColsRows
                )
              )))
            }
          }
          
          
          ###################################################################################################
          ### Close PDF
          ###################################################################################################
          try(dev.off());try(dev.off());try(dev.off());try(dev.off())
          
          
          ###################################################################################################
          ### Export heatmap data
          ###################################################################################################
          write.table(
            oDataPlotChrX[,names(oSamplesOrderCluster)],
            file=file.path(
              sPathCASCADE,
              paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
              paste0('Heatmaps__',sHMplotName ,'__data.txt')),
            sep='\t',
            quote=FALSE,
            col.names=TRUE,
            row.names=TRUE
          )
          
          
          ###################################################################################################
          ### Remove unused objects
          ###################################################################################################
          rm(oDataPlotChrX,oDataPlotChrXregionCol,oSamplesOrderCluster)
          gc();gc();gc();gc()
        })}
      }
      
      
      ###################################################################################################
      ### Export Control genes table
      ###################################################################################################
      sMeasureType <- '_Copies'
      sChrXRegion <- 'fullX'
      oResControlGenes <- NULL
      # for (sTrialID in c('CA27','CA34','CA35','CA36','CA43','CA63','CA76','CA79','CA83','PEACE001')) {
      for (sTrialID in c('CA63')) {
        load(file.path(
          sPathCASCADE,
          'details_bypatient',
          paste0(sTrialID,'__',as.integer(nAutosomesBinSize),'.RData')))
        rownames(oRes) <- oRes[,'region']
        load(paste0(
          sPathCASCADE,
          paste0('/S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
          '/Models__',sTrialID,'__',sSampleSel,'_',sMeasureType,'__',sChrXRegion,'__A',as.integer(nAutosomesBinSize),'__X',as.integer(nChrXBinSize),'.Rdata'))
        oResTmp <- NULL
        for (sControlGene in rownames(oControlGenes)) {
          oTmp <- as.data.frame(
            rowMeans(t(as.matrix(oRes[
              sapply(rownames(oRes),function(
    s,
    sCC=oControlGenes[sControlGene,'Chromosome'],
    nCS=oControlGenes[sControlGene,'Start'],
    nCE=oControlGenes[sControlGene,'End']
              )
                strsplit(s,':')[[1]][1]==sCC & 
      (
        (
          as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][1])>=nCS &
            as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][1])<=nCE
        ) |
          (
            as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][2])>=nCS &
              as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][2])<=nCE
          ) |
          (
            as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][1])>=nCS &
              as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][2])<=nCE
          ) |
          (
            as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][1])<=nCS &
              as.integer(strsplit(strsplit(s,':')[[1]][2],'-')[[1]][2])>=nCE
          )
      )
              ),
    grep(sMeasureType,colnames(oRes),value=TRUE)
            ])),
    na.rm=TRUE
            )
          )
          colnames(oTmp) <- sControlGene
          rownames(oTmp) <- gsub(sMeasureType,'',rownames(oTmp),fixed=TRUE)
          oTmp[,'Trial.ID'] <- sTrialID
          oTmp[,'Sample_ID'] <- oDataStats[rownames(oTmp),'Sample_ID']
          oTmp[,'Core_ID'] <- rownames(oTmp)
          if (is.null(oResTmp)) {
            oResTmp <- oTmp
          } else {
            oResTmp <- merge(oResTmp,oTmp,all=TRUE)
          }
        }
        oResControlGenes <- rbind(oResControlGenes,oResTmp)
      }
      write.table(
        oResControlGenes,
        file=file.path(
          sPathCASCADE,
          paste0('S',sSampleSel,'__A',as.integer(nAutosomesBinSize)),
          paste0('ControlGenes__data.txt')),
        sep='\t',
        quote=FALSE,
        col.names=TRUE,
        row.names=FALSE
      )
    }
  }
}


# setwd("C:/Data/UCL_stratosphere/reports/CASCADE/Smetsonly__A500000")
# load('Models__CA63__metsonly__TPs2__fullX__A500000__X500000.Rdata')
# lTreesTPs2 <- lTrees
# lOrganCols <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ.Color']
# rm(lTrees)
# load('Models__CA63__metsonly__Copies__fullX__A500000__X500000.Rdata')
# lTreesCopies <- lTrees
# lOrganCols <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ.Color']
# rm(lTrees)
# data <- lTreesTPs2[['data_aut']]
# hamming.dist  <- TRUE
# parsimony <- TRUE
# add.root <- FALSE
# return.dist <- TRUE
# dist.phylo <- lTreesTPs2[['data_aut']]
# dist.phylo <- lTreesTPs2[['data_aut']][]
# class(dist.phylo) <- "character"
# dist.phylo <- phyDat(
#     data = t(dist.phylo),
#     type = 'USER',
#     contrast=matrix(data = c(1,0,0,1,1,1),ncol = 2, byrow = TRUE, dimnames = list(c('0','1','?'),c('0','1'))),
#     ambiguity = c("?")
# )
# dist <- dist.hamming(dist.phylo, ratio = TRUE, exclude = "none")
# data.philo <- data
# data.philo[] <- data[]
# class(data.philo) <- "character"
# data.philo[is.na(data.philo)] <- '?'
# contrast <- matrix(data = c(1,0,0,1,1,1),ncol = 2, byrow = TRUE)
# dimnames(contrast) <- list(c('0','1','?'),c('0','1'))
# data.philo <- phyDat(
#     data = t(data.philo),
#     type = 'USER',
#     contrast=contrast,
#     ambiguity = c("?")
# )
# data.dist <- dist.hamming(data.philo, ratio = TRUE, exclude = "none")
# as.matrix(data.dist)[1:5,1:5]
# if (parsimony) {
#     tree <- optim.parsimony(NJ(data.dist[]), data.philo)
# } else {
#     tree <- NJ(data.dist)
# }
# if (add.root) {        
#     tree <- root(tree,outgroup = "root", resolve.root = TRUE)
# }
# plot(root(optim.parsimony(NJ(dist), dist.phylo),outgroup = "root", resolve.root = TRUE))
# lTreesCopies[['dist_aut']] <- fClustering(lTreesCopies[['data_aut']])
# sChrSource <- 'aut'
# nShift <- 2
# as.matrix(fTreePhylo(lTreesTPs2[[paste0('TPs',nShift,'_',sChrSource)]],hamming.dist=TRUE,add.root=FALSE,return.dist=TRUE))[1:5,1:5]
# as.matrix(lTreesTPs2[['dist_aut']])[1:5,1:5]
# as.matrix(dist)[1:5,1:5]
# sum(lTreesTPs2[['data_aut']][,'CA63_11'] != lTreesTPs2[['data_aut']][,'CA63_14']) / nrow(lTreesTPs2[['data_aut']])
# plot(fTreePhylo(lTreesTPs2[[paste0('TPs',nShift,'_',sChrSource)]],hamming.dist=FALSE,add.root=TRUE,return.dist=FALSE))
# tree <- optim.parsimony(NJ(dist), dist.phylo)
# tree <- root(tree,outgroup = "root", resolve.root = TRUE)
# plot(tree)
# boxplot(
#     as.numeric(as.matrix(lTreesTPs2[['dist_aut']])[lTreesTPs2[['info']][lTreesTPs2[['info']][,'Organ']=='Liver','Core_ID'],lTreesTPs2[['info']][lTreesTPs2[['info']][,'Organ']=='Liver','Core_ID']]),
#     as.numeric(as.matrix(lTreesCopies[['dist_aut']])[lTreesCopies[['info']][lTreesCopies[['info']][,'Organ']=='Liver','Core_ID'],lTreesCopies[['info']][lTreesCopies[['info']][,'Organ']=='Liver','Core_ID']]),
#     as.numeric(as.matrix(dist)[lTreesTPs2[['info']][lTreesTPs2[['info']][,'Organ']=='Liver','Core_ID'],lTreesTPs2[['info']][lTreesTPs2[['info']][,'Organ']=='Liver','Core_ID']])
#     )
# fSetLPdist(
#     data = oDataHCdist,
#     dist = lTreesCopies[['dist_aut']],
#     info = lTreesCopies[['info']],
#     source = sChrSource,
#     color = lOrganCols,
#     title = paste(
#     'dist',
#     sSampleSel,
#     as.integer(nAutosomesBinSize),
#     sMeasureType,
#     sChrXRegion,
#     sChrSource,
#     sep='_')
# )
# data <- oDataHCdist
# dist <- lTreesTPs2[[paste0('dist_',sChrSource)]]
# info <- lTreesTPs2[['info']]
# source <- sChrSource
# color <- lOrganCols
# title <- paste('dist',sSampleSel,as.integer(nAutosomesBinSize),sMeasureType,sChrXRegion,sChrSource,sep='_')

for (sSampleSet in names(lSampleSets)) {
  setwd(paste0('C:/Data/UCL_stratosphere/reports/CASCADE/S',sSampleSet,'__A500000'))
  pdf(paste0('Transition point distributions (all ',sSampleSet,').pdf'),paper='special',width=24,height=12)
  # for (sTrialID in c('CA27','CA34','CA35','CA36','CA43','CA63','CA76','CA79','CA83','PEACE001')) {
  for (sTrialID in c('CA63')) {
    
    for (sFileToLoadMain in intersect(
      list.files(".",pattern=paste0('Models__',sTrialID)),
      list.files(".",pattern=paste0('__Copies__fullX__A500000__X500000.Rdata'))
    )) {try({
      print(sFileToLoadMain)
      load(sFileToLoadMain)
      oDataPlot <- NULL
      lRegionGrop <- list()
      lSampleGrop <- list()
      for (sGroup in c(
        'All',
        'ARgain',
        'ARnormal',
        unique(lTrees[['info']][,'Organ']),
        unique(lTrees[['info']][,'clusters_str_aut']),
        sort(lTrees[['info']][lTrees[['info']][,'Figure_N']=='arch','Core_ID'],decreasing=TRUE)
      )) {
        if (sGroup == 'All') {
          lSampleGrop <- c(lSampleGrop,list(All=colnames(lTrees[['TPs2_aut']])))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'All',
            group = 'All',
            region = rownames(lTrees[['TPs2_aut']]),
            chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = ncol(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch','Core_ID'],drop=FALSE]),
            sample_match = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch','Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'All',
            group = 'All',
            region = rownames(lTrees[['TPs2_chrX']]),
            chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = ncol(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch','Core_ID'],drop=FALSE]),
            sample_match = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch','Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
        }
        if (sGroup == 'ARgain') {
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'AR',
            group = 'ARgain',
            region = rownames(lTrees[['TPs2_aut']]),
            chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = ncol(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']>2,'Core_ID'],drop=FALSE]),
            sample_match = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']>2,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'AR',
            group = 'ARgain',
            region = rownames(lTrees[['TPs2_chrX']]),
            chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = ncol(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']>2,'Core_ID'],drop=FALSE]),
            sample_match = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']>2,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
        }
        if (sGroup == 'ARnormal') {
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'AR',
            group = 'ARnormal',
            region = rownames(lTrees[['TPs2_aut']]),
            chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = ncol(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']<=2,'Core_ID'],drop=FALSE]),
            sample_match = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']<=2,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'AR',
            group = 'ARnormal',
            region = rownames(lTrees[['TPs2_chrX']]),
            chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = ncol(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']<=2,'Core_ID'],drop=FALSE]),
            sample_match = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'ARcnLAB']<=2,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
        }
        if (sGroup %in% unique(lTrees[['info']][,'Organ'])) {
          lSampleGrop <- c(lSampleGrop,list(tmp=colnames(lTrees[['TPs2_aut']])[colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Organ']==sGroup,'Core_ID']]))
          names(lSampleGrop)[names(lSampleGrop)=='tmp'] <- sGroup
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'Organ',
            group = sGroup,
            region = rownames(lTrees[['TPs2_aut']]),
            chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = sum(colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'Organ']==sGroup,'Core_ID']),
            sample_match = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'Organ',
            group = sGroup,
            region = rownames(lTrees[['TPs2_chrX']]),
            chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = sum(colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'Organ']==sGroup,'Core_ID']),
            sample_match = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'Organ']==sGroup,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
        }
        if (sGroup %in% unique(lTrees[['info']][,'clusters_str_aut'])) {
          lSampleGrop <- c(lSampleGrop,list(tmp=colnames(lTrees[['TPs2_aut']])[colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']]))
          names(lSampleGrop)[names(lSampleGrop)=='tmp'] <- sGroup
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'Cluster',
            group = sGroup,
            region = rownames(lTrees[['TPs2_aut']]),
            chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = sum(colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']),
            sample_match = rowSums(lTrees[['TPs2_aut']][,colnames(lTrees[['TPs2_aut']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'Cluster',
            group = sGroup,
            region = rownames(lTrees[['TPs2_chrX']]),
            chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = sum(colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID']),
            sample_match = rowSums(lTrees[['TPs2_chrX']][,colnames(lTrees[['TPs2_chrX']]) %in% lTrees[['info']][lTrees[['info']][,'Figure_N']!='arch' & lTrees[['info']][,'clusters_str_aut']==sGroup,'Core_ID'],drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
        }
        if (sGroup %in% lTrees[['info']][lTrees[['info']][,'Figure_N']=='arch','Core_ID']) {
          lSampleGrop <- c(lSampleGrop,list(tmp=sGroup))
          names(lSampleGrop)[names(lSampleGrop)=='tmp'] <- sGroup
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'Archival',
            group = sGroup,
            region = rownames(lTrees[['TPs2_aut']]),
            chr = sapply(rownames(lTrees[['TPs2_aut']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = sum(colnames(lTrees[['TPs2_aut']]) %in% sGroup),
            sample_match = rowSums(lTrees[['TPs2_aut']][,sGroup,drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
          oDataPlot <- rbind(oDataPlot,data.frame(
            group_type = 'Archival',
            group = sGroup,
            region = rownames(lTrees[['TPs2_chrX']]),
            chr = sapply(rownames(lTrees[['TPs2_chrX']]),function(s) strsplit(s,'_')[[1]][1]),
            sample_count = sum(colnames(lTrees[['TPs2_chrX']]) %in% sGroup),
            sample_match = rowSums(lTrees[['TPs2_chrX']][,sGroup,drop=FALSE]),
            count = 1,
            total = 0,
            cumsum_dec = 0,
            cumsum_inc = 0,
            stringsAsFactors = FALSE        
          ))
        }
        lRegionGrop <- c(lRegionGrop,list(tmp=unique(oDataPlot[
          oDataPlot[,'group']==sGroup &
            (oDataPlot[,'sample_match'] / oDataPlot[,'sample_count']) >= 0.8,
          'region'
        ])))
        names(lRegionGrop)[names(lRegionGrop)=='tmp'] <- sGroup
        oDataPlot <- oDataPlot[oDataPlot[,'sample_match']>0,]
        oDataPlot[oDataPlot[,'group']==sGroup,'total'] <- sum(oDataPlot[oDataPlot[,'group']==sGroup,'count'])
        oDataPlot <- oDataPlot[order(oDataPlot[,'sample_match'],decreasing=TRUE),]
        oDataPlot[oDataPlot[,'group']==sGroup,'cumsum_dec'] <- cumsum(oDataPlot[oDataPlot[,'group']==sGroup,'count'])
        oDataPlot <- oDataPlot[order(oDataPlot[,'sample_match'],decreasing=FALSE),]
        oDataPlot[oDataPlot[,'group']==sGroup,'cumsum_inc'] <- cumsum(oDataPlot[oDataPlot[,'group']==sGroup,'count'])
      }
      lRegionSample <- list()
      for (sSample in colnames(lTrees[['TPs2_aut']])) {
        lRegionSample <- c(lRegionSample,list(tmp=unique(rownames(lTrees[['TPs2_aut']])[lTrees[['TPs2_aut']][,sSample]==1])))
        names(lRegionSample)[names(lRegionSample)=='tmp'] <- sSample
      }
      lOrganCols <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ.Color']
      names(lOrganCols) <- unique(lTrees[['info_all']][,c('Organ','Organ.Color')])[,'Organ']
      grid.arrange(
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Organ'),]) + 
          geom_density_ridges(
            aes(x=sample_match / sample_count,y=interaction(group,group_type),fill=group)
          ) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=45, hjust=1)) + 
          ggtitle(sTrialID) + 
          xlab('Fraction of samples') + 
          ylab('Samples by Organ'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Organ'),]) + 
          geom_density_ridges(
            aes(x=sample_count,y=interaction(group,group_type),fill=group)
          ) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=45, hjust=1)) + 
          ggtitle(sTrialID) + 
          xlab('Number of samples') + 
          ylab('Samples by Organ'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Cluster'),]) + 
          geom_density_ridges(
            aes(x=sample_match / sample_count,y=interaction(group,group_type),fill=group)
          ) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=45, hjust=1)) + 
          ggtitle(sTrialID) + 
          xlab('Fraction of samples') + 
          ylab('Samples by Cluster'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Cluster'),]) + 
          geom_density_ridges(
            aes(x=sample_count,y=interaction(group,group_type),fill=group)
          ) + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=45, hjust=1)) + 
          ggtitle(sTrialID) + 
          xlab('Number of samples') + 
          ylab('Samples by Cluster'),
        ncol=2,
        nrow=2,
        layout_matrix=matrix(nrow=2,ncol=2,c(1,2,3,4))
      )
      grid.arrange(
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Organ'),]) + 
          geom_line(
            aes(x=sample_match,y=cumsum_dec,color=group,fill=group,group=group),
            position=position_dodge(width=0.1),
            size=1
          ) + 
          scale_x_reverse() + 
          scale_color_manual(name = 'Organ', values=c(All='black',lOrganCols), breaks=names(c(All='black',lOrganCols)),labels=names(c(All='black',lOrganCols))) +
          theme_bw() + 
          ggtitle(sTrialID) + 
          xlab('Number of samples') + 
          ylab('Cumulative number of transition points'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Organ'),]) + 
          geom_line(
            aes(x=sample_count/sample_match,y=cumsum_dec,color=group,fill=group,group=group),
            position=position_dodge(width=0.1),
            size=1
          ) + 
          scale_color_manual(name = 'Organ', values=c(All='black',lOrganCols), breaks=names(c(All='black',lOrganCols)),labels=names(c(All='black',lOrganCols))) +
          theme_bw() + 
          ggtitle(sTrialID) + 
          xlab('1 / Fraction of samples') + 
          ylab('Cumulative number of transition points'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','Organ'),]) + 
          geom_line(
            aes(x=((1/sample_match) - (1/sample_count)) / (1-1/sample_count) ,y=cumsum_dec / total,color=group,fill=group,group=group),
            position=position_dodge(width=0.01),
            size=1
          ) + 
          scale_color_manual(name = 'Organ', values=c(All='black',lOrganCols), breaks=names(c(All='black',lOrganCols)),labels=names(c(All='black',lOrganCols))) +
          theme_bw() + 
          ggtitle(sTrialID) + 
          xlab('1 / Fraction of samples (normalized to the 0-1 range)') + 
          ylab('Fraction of transition points'),
        ncol=1,
        nrow=3,
        layout_matrix=matrix(nrow=3,ncol=1,c(1,2,3))
      )
      grid.arrange(
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','AR'),]) + 
          geom_line(
            aes(x=sample_match,y=cumsum_dec,color=group,fill=group,group=group),
            position=position_dodge(width=0.1),
            size=1
          ) + 
          scale_x_reverse() + 
          theme_bw() + 
          ggtitle(sTrialID) + 
          xlab('Number of samples') + 
          ylab('Cumulative number of transition points'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','AR'),]) + 
          geom_line(
            aes(x=sample_count/sample_match,y=cumsum_dec,color=group,fill=group,group=group),
            position=position_dodge(width=0.1),
            size=1
          ) + 
          theme_bw() + 
          ggtitle(sTrialID) + 
          xlab('1 / Fraction of samples') + 
          ylab('Cumulative number of transition points'),
        ggplot(oDataPlot[oDataPlot[,'group_type'] %in% c('All','AR'),]) + 
          geom_line(
            aes(x=((1/sample_match) - (1/sample_count)) / (1-1/sample_count) ,y=cumsum_dec / total,color=group,fill=group,group=group),
            position=position_dodge(width=0.01),
            size=1
          ) + 
          theme_bw() + 
          ggtitle(sTrialID) + 
          xlab('1 / Fraction of samples (normalized to the 0-1 range)') + 
          ylab('Fraction of transition points'),
        ncol=1,
        nrow=3,
        layout_matrix=matrix(nrow=3,ncol=1,c(1,2,3))
      )
      try(oPlot <- upset(fromList(lRegionGrop[names(lRegionGrop) %in% c(unique(lTrees[['info']][,'Organ']),lTrees[['info']][lTrees[['info']][,'Figure_N']=='arch','Core_ID'])]), nsets = sum(names(lRegionGrop) %in% unique(lTrees[['info']][,'Organ'])), order.by = "freq",keep.order=TRUE))
      print(oPlot)
      try(oPlot <- upset(fromList(lRegionGrop[!names(lRegionGrop) %in% c('All','ARnormal','ARgain')]), sets = names(lRegionGrop)[!names(lRegionGrop) %in% c('All','ARnormal','ARgain')],nsets = sum(!names(lRegionGrop) %in% !names(lRegionGrop) %in% c('All','ARnormal','ARnormal')), order.by = "freq",keep.order=TRUE))
      print(oPlot)
      lRegionSample <- list()
      for (sSample in colnames(lTrees[['TPs2_aut']])) {
        lRegionSample <- c(lRegionSample,list(tmp=unique(rownames(lTrees[['TPs2_aut']])[lTrees[['TPs2_aut']][,sSample]==1])))
        names(lRegionSample)[names(lRegionSample)=='tmp'] <- sSample
      }
      for (sTreeName in c('tree_aut','tree_hamming_TPs2_TRUE_aut')) {try({
        oPlot <- upset(
          fromList(lRegionSample), 
          nsets = length(lRegionSample),
          nintersects = NA,
          sets = labels(as.dendrogram(lTrees[[sTreeName]]))[
            labels(as.dendrogram(lTrees[[sTreeName]])) %in% names(lRegionSample)
          ],
          mainbar.y.label = sTreeName,
          sets.bar.color = lTrees[['info']][labels(as.dendrogram(lTrees[[sTreeName]]))[
            labels(as.dendrogram(lTrees[[sTreeName]])) %in% names(lRegionSample)
          ],
          'Organ.Color'
          ],
          set.metadata = list(
            data = lTrees[['info']][
              names(lRegionSample),
              c('Core_ID','Organ','ARcnLAB')],
            plots = list(list(type = "heat",column = "ARcnLAB",assigne=20))
          ),
          order.by = "freq",
          keep.order = TRUE
        )
        print(oPlot)
      })}
      for (sFileToLoadArchival in intersect(
        list.files(".",pattern=paste0('Archival__',sTrialID)),
        list.files(".",pattern='.RData')
      )) {try({
        load(sFileToLoadArchival)
        print(sFileToLoadArchival)
        oDataArch <- oResComparison[['S2p']][,
                                             regexpr('_Copies',colnames(oResComparison[['S2p']]),fixed=TRUE) > -1 
        ]
        oDataArch <- !is.na(as.matrix(oDataArch))
        class(oDataArch) <- "numeric"
        colnames(oDataArch) <- gsub('_Copies','',colnames(oDataArch),fixed=TRUE)
        oDataArch <- melt(oDataArch)
        colnames(oDataArch) <- c('TP','Core_ID','TPcount')
        oDataArch <- merge(oDataArch,oDataStats[,c('Trial.ID','Sample_ID','Core_ID','Figure_N','Tissue.Site','Organ','Organ.Color')])
        lOrganCols <- unique(oDataArch[,c('Organ','Organ.Color')])[,'Organ.Color']
        names(lOrganCols) <- unique(oDataArch[,c('Organ','Organ.Color')])[,'Organ']
        oDataArch <- merge(
          oDataArch,
          data.frame(
            Core_ID = colnames(lTrees[['TPs2_aut']]),
            TP_count_aut = colSums(lTrees[['TPs2_aut']]),
            TP_count_chrX = colSums(lTrees[['TPs2_chrX']]),
            stringsAsFactors = FALSE
          ),
          all.x=TRUE
        )
        oDataArch <- merge(
          oDataArch,
          data.frame(
            Organ = names(lapply(lRegionGrop,length)),
            TP_Organ = unlist(lapply(lRegionGrop,length)),
            stringsAsFactors = FALSE
          ),
          all.x=TRUE
        )
        oDataArch <- oDataArch[order(
          oDataArch[,'Organ'],
          as.integer(gsub('arch','-1',oDataArch[,'Figure_N']))
        ),]
        lSamplesOrder <- unique(paste(
          oDataArch[,'Figure_N'],
          oDataArch[,'Organ'],
          sep = '.'
        ))
        lSamplesOrder <- c(
          lSamplesOrder[regexpr('arch.',lSamplesOrder,fixed=TRUE)>-1],
          lSamplesOrder[regexpr('arch.',lSamplesOrder,fixed=TRUE)==-1]
        )
        oPlot <- ggplot()
        oPlot <- oPlot + geom_bar(
          data=unique(oDataArch[,c('Figure_N','Organ','TP_count_aut','TP_count_chrX')]),aes(
            x = interaction(Figure_N,Organ),
            y = TP_count_aut + TP_count_chrX,
            group= Organ,
            fill = 'Private'),
          color='black',
          stat="identity")
        oPlot <- oPlot + geom_bar(
          data=unique(oDataArch[,c('Figure_N','Organ','TP_Organ')]),aes(
            x = interaction(Figure_N,Organ),
            y = ifelse(Figure_N=='arch',0,TP_Organ),
            group= Organ,
            fill = Organ
          ),
          color='black',
          stat="identity")
        oPlot <- oPlot + geom_bar(
          data=oDataArch,aes(
            x = interaction(Figure_N,Organ),
            y = TPcount,
            group= Organ),
          fill = 'black',
          position="stack", 
          color='black',
          stat="identity")
        oPlot <- oPlot + scale_fill_manual(
          name = 'Organ',
          values=c(Private='white',Archival='black',lOrganCols),
          breaks=names(c(Private='white',Archival='black',lOrganCols)),
          labels=names(c(Private='white',Archival='black',lOrganCols))
        )
        oPlot <- oPlot + scale_x_discrete(limits=lSamplesOrder)
        oPlot <- oPlot + ggtitle(sTrialID)
        oPlot <- oPlot + xlab('Sample / Organ')
        oPlot <- oPlot + ylab('Transition point count')
        oPlot <- oPlot + theme_bw()
        oPlot <- oPlot + theme(axis.text.x=element_text(angle=45, hjust=1))
        plot(oPlot)
        rm(oResComparison)
      })}
      rm(lTrees)
    })}
  }
  try(dev.off());try(dev.off());try(dev.off());try(dev.off())
}



