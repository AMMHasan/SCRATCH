require(QDNAseq)
require(QDNAseq.hg19)
require(ACE)
require(Biobase)
require(BSgenome.Hsapiens.UCSC.hg19)
require(CGHbase)
require(CGHcall)
require(tidyverse)


options(stringsAsFactors = FALSE)



# PCE_list <- list()

#' generate CN df
#'
#' @param BAM_path path of .bam files
#' @param pattern_PEA pattern for the project samples
#' @param bin_size bin size to calculate copy number by QDNAseq
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM/"
#' bin_size=500
#' pattern_PEA = "PEA310"
#' generate_CN_df(BAM_path,pattern_PEA, bin_size)
#' 
#' 
generate_CN_df <- function(BAM_path,pattern_PEA, bin_size){
  
  sampleFileNames <- setdiff(
    list.files(path = BAM_path, 
               pattern = pattern_PEA), 
    list.files(path = BAM_path, 
               pattern = ".bai")
  )
  
  
  bamFiles <- file.path(BAM_path, sampleFileNames)
  
  
  bins <- getBinAnnotations(binSize=bin_size) 
  
  readCounts <- binReadCounts(
    bins = bins,
    bamfiles=bamFiles, 
    path=NULL, # as the bamfiles contain the full path 
    ext="bam", 
    bamnames=NULL, 
    phenofile=NULL,
    chunkSize=10000000, 
    cache=TRUE, 
    isPaired=NA,
    isProperPair=NA, 
    isUnmappedQuery=FALSE, 
    hasUnmappedMate=NA, 
    isMinusStrand=NA,
    isMateMinusStrand=NA, 
    isFirstMateRead=NA, 
    isSecondMateRead=NA, 
    isSecondaryAlignment=NA,
    isNotPassingQualityControls=FALSE, 
    isDuplicate=FALSE, 
    minMapq=37, 
    pairedEnds=NULL
  )
  
  
  readCountsFiltered <- applyFilters(
    readCounts, 
    residual=TRUE, 
    blacklist=TRUE,
    mappability=NA, 
    bases=NA,
    chromosomes=c("Y", "MT")
  )
  
  readCountsFiltered <- estimateCorrection(
    readCountsFiltered, 
    span=0.65, 
    family="symmetric", 
    adjustIncompletes=TRUE,
    maxIter=1, 
    cutoff=4, 
    variables=c("gc", "mappability")
  )
  
  copyNumbers <- correctBins(readCountsFiltered)
  
  copyNumbersNormalized <- normalizeBins(
    copyNumbers,
    method="median", 
    force=FALSE)
  
  
  copyNumbersSmooth <- smoothOutlierBins(
    copyNumbersNormalized,
    logTransform=TRUE, 
    force=FALSE)
  
  copyNumbersSegmented <- segmentBins(
    copyNumbersSmooth, 
    transformFun="log2",
    smoothBy=FALSE, 
    alpha=1e-10, 
    undo.splits="sdundo", 
    undo.SD=1,
    force=FALSE
    # seeds=NULL
  )
  
  
  copyNumbersSegmented <- normalizeSegmentedBins(
    copyNumbersSegmented,
    force = F
  )
  
  
  copyNumbersCalled <- callBins(
    copyNumbersSegmented, 
    organism="human", 
    method="CGHcall",
    prior = "all"
    #cutoffs=log2(c(deletion = 0.5, loss = 1.5, gain = 2.5, amplification = 10)/2)
  )
  
  cgh <- makeCgh(copyNumbersCalled)
  
  CN_df <- CGHbase::segmented(cgh) %>% 
    data.frame %>% 
    mutate(position = rownames(.)) 
  return(CN_df)
}



BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM/"
bin_size=500

#######################
#######################               PEA310
#######################

pattern_PEA = "PEA310"

# CN_df <- generate_CN_df(BAM_path,pattern_PEA, bin_size)
# 
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_10044875) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_10044875)))
# TP_PEA310_10044875 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_10044916) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_10044916)))
# TP_PEA310_10044916 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_10044940) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_10044940)))
# TP_PEA310_10044940 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_10044953) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_10044953)))
# TP_PEA310_10044953 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_10044958) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_10044958)))
# TP_PEA310_10044958 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_2018_10_17) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_2018_10_17)))
# TP_PEA310_2018_10_17 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_2021_01_26) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_2021_01_26)))
# TP_PEA310_2021_01_26 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA310_2022_03_30) %>% 
#   mutate(diff = c(NA,diff(sub_PEA310_2022_03_30)))
# TP_PEA310_2022_03_30 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_list <- list()
# TP_list[["PEA310_10044875"]] <- TP_PEA310_10044875
# TP_list[["PEA310_10044916"]] <- TP_PEA310_10044916
# TP_list[["PEA310_10044940"]] <- TP_PEA310_10044940
# TP_list[["PEA310_10044953"]] <- TP_PEA310_10044953
# TP_list[["PEA310_10044958"]] <- TP_PEA310_10044958
# TP_list[["PEA310_2018_10_17"]] <- TP_PEA310_2018_10_17
# TP_list[["PEA310_2021_01_26"]] <- TP_PEA310_2021_01_26 
# TP_list[["PEA310_2022_03_30"]] <- TP_PEA310_2022_03_30 
# 
# TP_binary_matrix <- t(splitstackshape:::charMat(listOfValues = TP_list, fill = 0L))
# colnames(TP_binary_matrix) <- names(TP_list)
# 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA310_2018_10_17 > 0) %>% 
#   select(starts_with("PEA310_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA310_2021_01_26 > 0) %>% 
#   select(starts_with("PEA310_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA310_2022_03_30 > 0) %>% 
#   select(starts_with("PEA310_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100
#   
#   
# ################################################################################
# 
# #######################
# #######################               PEA236
# #######################
# 
# pattern_PEA = "PEA236"
# 
# CN_df <- generate_CN_df(BAM_path,pattern_PEA, bin_size)
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_10024218) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_10024218)))
# TP_PEA236_10024218 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_10024219) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_10024219)))
# TP_PEA236_10024219 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_10024225) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_10024225)))
# TP_PEA236_10024225 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_10024236) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_10024236)))
# TP_PEA236_10024236 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_2019_11_11) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_2019_11_11)))
# TP_PEA236_2019_11_11 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_2019_12_11) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_2019_12_11)))
# TP_PEA236_2019_12_11 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA236_2020_01_08) %>% 
#   mutate(diff = c(NA,diff(sub_PEA236_2020_01_08)))
# TP_PEA236_2020_01_08 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# 
# TP_list <- list()
# 
# TP_list[["PEA236_10024218"]] <- TP_PEA236_10024218
# TP_list[["PEA236_10024219"]] <- TP_PEA236_10024219
# TP_list[["PEA236_10024225"]] <- TP_PEA236_10024225
# TP_list[["PEA236_10024236"]] <- TP_PEA236_10024236
# TP_list[["PEA236_2019_11_11"]] <- TP_PEA236_2019_11_11
# TP_list[["PEA236_2019_12_11"]] <- TP_PEA236_2019_12_11
# TP_list[["PEA236_2020_01_08"]] <- TP_PEA236_2020_01_08 
# 
# 
# TP_binary_matrix <- t(splitstackshape:::charMat(listOfValues = TP_list, fill = 0L))
# colnames(TP_binary_matrix) <- names(TP_list)
# 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA236_2019_11_11 > 0) %>% 
#   select(starts_with("PEA236_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA236_2019_12_11 > 0) %>% 
#   select(starts_with("PEA236_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA236_2020_01_08 > 0) %>% 
#   select(starts_with("PEA236_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# ##################################################
# 
# #######################
# #######################               PEA307
# #######################
# 
# 
# pattern_PEA = "PEA307"
# 
# CN_df <- generate_CN_df(BAM_path,pattern_PEA, bin_size)
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_10052064) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_10052064)))
# TP_PEA307_10052064 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_10052074) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_10052074)))
# TP_PEA307_10052074 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_10052088) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_10052088)))
# TP_PEA307_10052088 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_10052111) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_10052111)))
# TP_PEA307_10052111 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_10052131) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_10052131)))
# TP_PEA307_10052131 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_10052145) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_10052145)))
# TP_PEA307_10052145 <- TP_df[which(TP_df$diff != 0),]$position
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_2021_04_16) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_2021_04_16)))
# TP_PEA307_2021_04_16 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA307_2020_11_04) %>% 
#   mutate(diff = c(NA,diff(sub_PEA307_2020_11_04)))
# TP_PEA307_2020_11_04 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_list <- list()
# 
# TP_list[["PEA307_10052064"]] <- TP_PEA307_10052064
# TP_list[["PEA307_10052074"]] <- TP_PEA307_10052074
# TP_list[["PEA307_10052088"]] <- TP_PEA307_10052088
# TP_list[["PEA307_10052111"]] <- TP_PEA307_10052111
# TP_list[["PEA307_10052131"]] <- TP_PEA307_10052131
# TP_list[["PEA307_10052145"]] <- TP_PEA307_10052145
# TP_list[["PEA307_2021_04_16"]] <- TP_PEA307_2021_04_16
# TP_list[["PEA307_2020_11_04"]] <- TP_PEA307_2020_11_04 
# 
# 
# TP_binary_matrix <- t(splitstackshape:::charMat(listOfValues = TP_list, fill = 0L))
# colnames(TP_binary_matrix) <- names(TP_list)
# 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA307_2021_04_16 > 0) %>% 
#   select(starts_with("PEA307_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA307_2020_11_04 > 0) %>% 
#   select(starts_with("PEA307_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# 
# 
# #######################
# #######################               PEA349
# #######################
# 
# 
# pattern_PEA = "PEA349"
# 
# CN_df <- generate_CN_df(BAM_path,pattern_PEA, bin_size)
# 
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA349_10045150) %>% 
#   mutate(diff = c(NA,diff(sub_PEA349_10045150)))
# TP_PEA349_10045150 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA349_10045161) %>% 
#   mutate(diff = c(NA,diff(sub_PEA349_10045161)))
# TP_PEA349_10045161 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA349_10045168) %>% 
#   mutate(diff = c(NA,diff(sub_PEA349_10045168)))
# TP_PEA349_10045168 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA349_10045199) %>% 
#   mutate(diff = c(NA,diff(sub_PEA349_10045199)))
# TP_PEA349_10045199 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA349_2021_06_09) %>% 
#   mutate(diff = c(NA,diff(sub_PEA349_2021_06_09)))
# TP_PEA349_2021_06_09 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_df <- CN_df %>% 
#   select(position, sub_PEA349_2022_04_06) %>% 
#   mutate(diff = c(NA,diff(sub_PEA349_2022_04_06)))
# TP_PEA349_2022_04_06 <- TP_df[which(TP_df$diff != 0),]$position
# 
# 
# TP_list <- list()
# 
# TP_list[["PEA349_10045150"]] <- TP_PEA349_10045150
# TP_list[["PEA349_10045161"]] <- TP_PEA349_10045161
# TP_list[["PEA349_10045168"]] <- TP_PEA349_10045168
# TP_list[["PEA349_10045199"]] <- TP_PEA349_10045199
# TP_list[["PEA349_2021_06_09"]] <- TP_PEA349_2021_06_09
# TP_list[["PEA349_2022_04_06"]] <- TP_PEA349_2022_04_06
# 
# TP_binary_matrix <- t(splitstackshape:::charMat(listOfValues = TP_list, fill = 0L))
# colnames(TP_binary_matrix) <- names(TP_list)
# 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA349_2021_06_09 > 0) %>% 
#   select(starts_with("PEA349_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# 
# 
# TP_binary_matrix %>% 
#   as.data.frame() %>% 
#   filter(PEA349_2022_04_06 > 0) %>% 
#   select(starts_with("PEA349_1")) %>% 
#   mutate(row_sum = rowSums(.)) %>% 
#   pull(row_sum) %>% table %>% prop.table()*100 
# 
# 
# 
# 
# 
# 
# ##################################################
# plot_data <- rio::import("~/Documents/Research/prostate_cancer/PEACE/TP/seq_plasma_vs_PM_tumour.xlsx") %>% 
#   filter(`Number of mets` > 0)
# 
# ggplot(plot_data) +
#   geom_line(aes(x=`Number of mets`, y=PEA310_2018_10_17), colour = "red", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA310_2021_01_26), colour = "black", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA310_2022_03_30), colour = "darkgreen", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA236_2019_11_11), colour = "blue", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA236_2019_12_11), colour = "gold", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA236_2020_01_08), colour = "brown", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA307_2021_04_16), colour = "pink", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA307_2020_11_04), colour = "grey", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA349_2021_06_09), colour = "violet", size=1.5) +
#   geom_line(aes(x=`Number of mets`, y=PEA349_2022_04_06), colour = "orange", size=1.5) +
#   ylab("Percent of plasma TP found in autopsy cores") +
#   theme_bw()
# 
# 
# ##################################################
# ##################### PEA310 #####################
# ##################################################
# # PEA310_10044875
# length(intersect(TP_PEA310_10044875,TP_PEA310_2018_10_17))*100/length(TP_PEA310_10044875)
# length(intersect(TP_PEA310_10044875,TP_PEA310_2021_01_26))*100/length(TP_PEA310_10044875)
# length(intersect(TP_PEA310_10044875,TP_PEA310_2022_03_30))*100/length(TP_PEA310_10044875)
# 
# # PEA310_10044916
# length(intersect(TP_PEA310_10044916,TP_PEA310_2018_10_17))*100/length(TP_PEA310_10044916)
# length(intersect(TP_PEA310_10044916,TP_PEA310_2021_01_26))*100/length(TP_PEA310_10044916)
# length(intersect(TP_PEA310_10044916,TP_PEA310_2022_03_30))*100/length(TP_PEA310_10044916)
# 
# # PEA310_10044940
# length(intersect(TP_PEA310_10044940,TP_PEA310_2018_10_17))*100/length(TP_PEA310_10044940)
# length(intersect(TP_PEA310_10044940,TP_PEA310_2021_01_26))*100/length(TP_PEA310_10044940)
# length(intersect(TP_PEA310_10044940,TP_PEA310_2022_03_30))*100/length(TP_PEA310_10044940)
# 
# # PEA310_10044953
# length(intersect(TP_PEA310_10044953,TP_PEA310_2018_10_17))*100/length(TP_PEA310_10044953)
# length(intersect(TP_PEA310_10044953,TP_PEA310_2021_01_26))*100/length(TP_PEA310_10044953)
# length(intersect(TP_PEA310_10044953,TP_PEA310_2022_03_30))*100/length(TP_PEA310_10044953)
# 
# # PEA310_10044958
# length(intersect(TP_PEA310_10044958,TP_PEA310_2018_10_17))*100/length(TP_PEA310_10044958)
# length(intersect(TP_PEA310_10044958,TP_PEA310_2021_01_26))*100/length(TP_PEA310_10044958)
# length(intersect(TP_PEA310_10044958,TP_PEA310_2022_03_30))*100/length(TP_PEA310_10044958)
# 
# 
# ##################################################
# ##################### PEA236 #####################
# ##################################################
# 
# # PEA236_10024218
# length(intersect(TP_PEA236_10024218,TP_PEA236_2019_11_11))*100/length(TP_PEA236_10024218)
# length(intersect(TP_PEA236_10024218,TP_PEA236_2019_12_11))*100/length(TP_PEA236_10024218)
# length(intersect(TP_PEA236_10024218,TP_PEA236_2020_01_08))*100/length(TP_PEA236_10024218)
# 
# # PEA236_10024219
# length(intersect(TP_PEA236_10024219,TP_PEA236_2019_11_11))*100/length(TP_PEA236_10024219)
# length(intersect(TP_PEA236_10024219,TP_PEA236_2019_12_11))*100/length(TP_PEA236_10024219)
# length(intersect(TP_PEA236_10024219,TP_PEA236_2020_01_08))*100/length(TP_PEA236_10024219)
# 
# # PEA236_10024225
# length(intersect(TP_PEA236_10024225,TP_PEA236_2019_11_11))*100/length(TP_PEA236_10024225)
# length(intersect(TP_PEA236_10024225,TP_PEA236_2019_12_11))*100/length(TP_PEA236_10024225)
# length(intersect(TP_PEA236_10024225,TP_PEA236_2020_01_08))*100/length(TP_PEA236_10024225)
# 
# # PEA236_10024236
# length(intersect(TP_PEA236_10024236,TP_PEA236_2019_11_11))*100/length(TP_PEA236_10024236)
# length(intersect(TP_PEA236_10024236,TP_PEA236_2019_12_11))*100/length(TP_PEA236_10024236)
# length(intersect(TP_PEA236_10024236,TP_PEA236_2020_01_08))*100/length(TP_PEA236_10024236)
# 
# ##################################################
# ##################### PEA307 #####################
# ##################################################
# 
# # PEA307_10052064   
# length(intersect(TP_PEA307_10052064,TP_PEA307_2020_11_04))*100/length(TP_PEA307_10052064)
# length(intersect(TP_PEA307_10052064,TP_PEA307_2021_04_16))*100/length(TP_PEA307_10052064)
# 
# # PEA307_10052074   
# length(intersect(TP_PEA307_10052074,TP_PEA307_2020_11_04))*100/length(TP_PEA307_10052074)
# length(intersect(TP_PEA307_10052074,TP_PEA307_2021_04_16))*100/length(TP_PEA307_10052074)
# 
# # PEA307_10052088   
# length(intersect(TP_PEA307_10052088,TP_PEA307_2020_11_04))*100/length(TP_PEA307_10052088)
# length(intersect(TP_PEA307_10052088,TP_PEA307_2021_04_16))*100/length(TP_PEA307_10052088)
# 
# # PEA307_100521111  
# length(intersect(TP_PEA307_10052111,TP_PEA307_2020_11_04))*100/length(TP_PEA307_10052111)
# length(intersect(TP_PEA307_10052111,TP_PEA307_2021_04_16))*100/length(TP_PEA307_10052111)
# 
# # PEA307_10052131   
# length(intersect(TP_PEA307_10052131,TP_PEA307_2020_11_04))*100/length(TP_PEA307_10052131)
# length(intersect(TP_PEA307_10052131,TP_PEA307_2021_04_16))*100/length(TP_PEA307_10052131)
# 
# # PEA307_10052145   
# length(intersect(TP_PEA307_10052145,TP_PEA307_2020_11_04))*100/length(TP_PEA307_10052145)
# length(intersect(TP_PEA307_10052145,TP_PEA307_2021_04_16))*100/length(TP_PEA307_10052145)
# 
# ##################################################
# ##################### PEA349 #####################
# ##################################################
# 
# # PEA349_10045150   
# length(intersect(TP_PEA349_10045150,TP_PEA349_2021_06_09))*100/length(TP_PEA349_10045150)
# length(intersect(TP_PEA349_10045150,TP_PEA349_2022_04_06))*100/length(TP_PEA349_10045150)
# 
# # PEA349_10045161   
# length(intersect(TP_PEA349_10045161,TP_PEA349_2021_06_09))*100/length(TP_PEA349_10045161)
# length(intersect(TP_PEA349_10045161,TP_PEA349_2022_04_06))*100/length(TP_PEA349_10045161)
# 
# # PEA349_10045168   
# length(intersect(TP_PEA349_10045168,TP_PEA349_2021_06_09))*100/length(TP_PEA349_10045168)
# length(intersect(TP_PEA349_10045168,TP_PEA349_2022_04_06))*100/length(TP_PEA349_10045168)
# 
# # PEA349_10045199   
# length(intersect(TP_PEA349_10045199,TP_PEA349_2021_06_09))*100/length(TP_PEA349_10045199)
# length(intersect(TP_PEA349_10045199,TP_PEA349_2022_04_06))*100/length(TP_PEA349_10045199)
# 
# 
# 
# 
# # 
# # for(i in seq(1,length(sampleNames))){
# #   model <- squaremodel(
# #     copyNumbersCalled,
# #     QDNAseqobjectsample = i, # In if (QDNAseqobjectsample) the condition has length > 1 and only the first element will be used
# #     prows=100,
# #     ptop=5,
# #     pbottom=1,
# #     method = 'RMSE',
# #     exclude = c("Y","MT"),
# #     penalty = 0.5,
# #     penploidy = 0.5,
# #     highlightminima = TRUE
# #   )
# #   
# #   PCE_list[[sampleNames[i]]]  <- model$minimadf[1,1:3]
# # }
# # 
# # 
# # 
# # ########################
# 
# 
# library(ComplexHeatmap)
# library(dplyr)
# 
# TP_plot <- rio::import("~/Documents/Research/prostate_cancer/PEACE/TP/TP_common_to_mets_and_plasma.csv")
# 
# # mat1 = matrix(TP_plot %>% 
# #                 filter(PEACE_ID == "PEA310") %>% 
# #                 pull(sample_ID) %>%
# #                 gsub("Plasma1_","",.) %>% 
# #                 gsub("Plasma2_","",.) %>% 
# #                 gsub("Plasma3_","",.) %>% 
# #                 unique(), ncol = 1)
# # 
# # rownames(mat1) <- TP_plot %>% 
# #   filter(PEACE_ID == "PEA310") %>% 
# #   pull(sample_ID) %>%
# #   gsub("Plasma1_","",.) %>% 
# #   gsub("Plasma2_","",.) %>% 
# #   gsub("Plasma3_","",.) %>% 
# #   unique()
# 
# 
# 
# #PEA236
# patientID="PEA236"
# 
# mat <- matrix(TP_plot %>% 
#                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>%
#                 pull(organ), ncol = 1)
# 
# 
# row_ha = rowAnnotation(plasma1 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm")),
#                        plasma2 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma2",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm")),
#                        plasma3 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma3",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm"))
# )
# 
# 
# Heatmap(mat, name = "Organs",  width = ncol(mat)*unit(5, "mm"), 
#         column_title = "organs", column_title_rot = 45,
#         col = setNames(TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ_colour), TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ)),
#         right_annotation = row_ha, row_labels = TP_plot %>% 
#           filter(PEACE_ID == patientID) %>% 
#           pull(sample_ID) %>%
#           gsub("Plasma1_","",.) %>% 
#           gsub("Plasma2_","",.) %>% 
#           gsub("Plasma3_","",.) %>% 
#           unique())
# 
# #PEA307
# 
# patientID="PEA307"
# mat <- matrix(TP_plot %>% 
#                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>%
#                 pull(organ), ncol = 1)
# 
# 
# row_ha = rowAnnotation(plasma1 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm")),
#                        plasma2 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma2",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm"))
# )
# 
# 
# Heatmap(mat, name = "Organs",  width = ncol(mat)*unit(5, "mm"), 
#         column_title = "organs", column_title_rot = 45,
#         col = setNames(TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ_colour), TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ)),
#         right_annotation = row_ha, row_labels = TP_plot %>% 
#           filter(PEACE_ID == patientID) %>% 
#           pull(sample_ID) %>%
#           gsub("Plasma1_","",.) %>% 
#           gsub("Plasma2_","",.) %>% 
#           gsub("Plasma3_","",.) %>% 
#           unique())
# 
# 
# 
# 
# #PEA310
# 
# patientID="PEA310"
# mat <- matrix(TP_plot %>% 
#                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>%
#                 pull(organ), ncol = 1)
# 
# 
# row_ha = rowAnnotation(plasma1 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm")),
#                        plasma2 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma2",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm")),
#                        plasma3 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma3",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm"))
# )
# 
# 
# Heatmap(mat, name = "Organs",  width = ncol(mat)*unit(5, "mm"), 
#         column_title = "organs", column_title_rot = 45,
#         col = setNames(TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ_colour), TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ)),
#         right_annotation = row_ha, row_labels = TP_plot %>% 
#           filter(PEACE_ID == patientID) %>% 
#           pull(sample_ID) %>%
#           gsub("Plasma1_","",.) %>% 
#           gsub("Plasma2_","",.) %>% 
#           gsub("Plasma3_","",.) %>% 
#           unique())
# 
# 
# 
# #PEA349
# 
# patientID="PEA349"
# mat <- matrix(TP_plot %>% 
#                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>%
#                 pull(organ), ncol = 1)
# 
# 
# row_ha = rowAnnotation(plasma1 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma1",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm")),
#                        plasma2 = anno_barplot(TP_plot %>% 
#                                                 filter(PEACE_ID == patientID & grepl("Plasma2",sample_ID)) %>% 
#                                                 pull(clonal_mutations_also_found_in_plasma), 
#                                               ylim = c(0,100), width = unit(20, "mm"))
# )
# 
# 
# Heatmap(mat, name = "Organs",  width = ncol(mat)*unit(5, "mm"), 
#         column_title = "organs", column_title_rot = 45,
#         col = setNames(TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ_colour), TP_plot %>% filter(PEACE_ID == patientID & Plasma == "Plasma1") %>% pull(organ)),
#         right_annotation = row_ha, row_labels = TP_plot %>% 
#           filter(PEACE_ID == patientID) %>% 
#           pull(sample_ID) %>%
#           gsub("Plasma1_","",.) %>% 
#           gsub("Plasma2_","",.) %>% 
#           gsub("Plasma3_","",.) %>% 
#           unique())
# 
# 
# 
# ######################
# 
# 
# # PEA236_10024218
# 
# splitstackshape:::charMat(list(intersect(TP_PEA236_10024218,TP_PEA236_2019_11_11),
#                                intersect(TP_PEA236_10024218,TP_PEA236_2019_12_11),
#                                intersect(TP_PEA236_10024218,TP_PEA236_2020_01_08)
# ), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# 
# # PEA236_10024219
# 
# splitstackshape:::charMat(list(intersect(TP_PEA236_10024219,TP_PEA236_2019_11_11),
#                                intersect(TP_PEA236_10024219,TP_PEA236_2019_12_11),
#                                intersect(TP_PEA236_10024219,TP_PEA236_2020_01_08)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# 
# # PEA236_10024225
# splitstackshape:::charMat(list(intersect(TP_PEA236_10024225,TP_PEA236_2019_11_11),
#                                intersect(TP_PEA236_10024225,TP_PEA236_2019_12_11),
#                                intersect(TP_PEA236_10024225,TP_PEA236_2020_01_08)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# 
# # PEA236_10024236
# splitstackshape:::charMat(list(intersect(TP_PEA236_10024236,TP_PEA236_2019_11_11),
#                                intersect(TP_PEA236_10024236,TP_PEA236_2019_12_11),
#                                intersect(TP_PEA236_10024236,TP_PEA236_2020_01_08)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# 
# ##################################################
# ##################### PEA307 #####################
# ##################################################
# 
# # PEA307_10052064 
# splitstackshape:::charMat(list(intersect(TP_PEA307_10052064,TP_PEA307_2020_11_04),
#                                intersect(TP_PEA307_10052064,TP_PEA307_2021_04_16)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA307_10052074   
# splitstackshape:::charMat(list(intersect(TP_PEA307_10052074,TP_PEA307_2020_11_04),
#                                intersect(TP_PEA307_10052074,TP_PEA307_2021_04_16)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# 
# # PEA307_10052088   
# splitstackshape:::charMat(list(intersect(TP_PEA307_10052088,TP_PEA307_2020_11_04),
#                                intersect(TP_PEA307_10052088,TP_PEA307_2021_04_16)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA307_100521111  
# splitstackshape:::charMat(list(intersect(TP_PEA307_10052111,TP_PEA307_2020_11_04),
#                                intersect(TP_PEA307_10052111,TP_PEA307_2021_04_16)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA307_10052131   
# splitstackshape:::charMat(list(intersect(TP_PEA307_10052131,TP_PEA307_2020_11_04),
#                                intersect(TP_PEA307_10052131,TP_PEA307_2021_04_16)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA307_10052145   
# splitstackshape:::charMat(list(intersect(TP_PEA307_10052145,TP_PEA307_2020_11_04),
#                                intersect(TP_PEA307_10052145,TP_PEA307_2021_04_16)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# ##################################################
# ##################### PEA349 #####################
# ##################################################
# 
# # PEA349_10045150   
# splitstackshape:::charMat(list(intersect(TP_PEA349_10045150,TP_PEA349_2021_06_09),
#                                intersect(TP_PEA349_10045150,TP_PEA349_2022_04_06)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA349_10045161   
# splitstackshape:::charMat(list(intersect(TP_PEA349_10045161,TP_PEA349_2021_06_09),
#                                intersect(TP_PEA349_10045161,TP_PEA349_2022_04_06)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA349_10045168   
# splitstackshape:::charMat(list(intersect(TP_PEA349_10045168,TP_PEA349_2021_06_09),
#                                intersect(TP_PEA349_10045168,TP_PEA349_2022_04_06)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 
# 
# # PEA349_10045199   
# splitstackshape:::charMat(list(intersect(TP_PEA349_10045199,TP_PEA349_2021_06_09),
#                                intersect(TP_PEA349_10045199,TP_PEA349_2022_04_06)), fill = 0L) %>% t() %>% ggparcoord(scale = "globalminmax")
# 



