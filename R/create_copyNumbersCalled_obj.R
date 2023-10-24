#' generate copy number call object
#'
#' @param BAM_path path of .bam files
#' @param pattern pattern for the project samples
#' @param bin_size bin size to calculate copy number by QDNAseq
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @import ACE
#' @import Biobase
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import CGHbase
#' @import CGHcall
#' @import tidyverse
#' @import dplyr
#' @import magrittr
#' @import QDNAseq
#' @import QDNAseq.hg19
#' @return a QDNAseq object
#' @export
#'
#' @examples
#' BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM_subset/"
#' bin_size=500
#' pattern = "PEA310"
#' generate_copyNumberCalled_obj(BAM_path,pattern, bin_size)
#' 
#' 
#'
generate_copyNumberCalled_obj <- function(BAM_path,pattern, bin_size){
  
  sampleFileNames <- setdiff(
    list.files(path = BAM_path, 
               pattern = pattern), 
    list.files(path = BAM_path, 
               pattern = ".bai")
  )
  
  
  bamFiles <- file.path(BAM_path, sampleFileNames)
  
  
  bins <- QDNAseq::getBinAnnotations(binSize=bin_size) 
  
  readCounts <- QDNAseq::binReadCounts(
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
  
  
  readCountsFiltered <- QDNAseq::applyFilters(
    readCounts, 
    residual=TRUE, 
    blacklist=TRUE,
    mappability=NA, 
    bases=NA,
    chromosomes=c("Y", "MT")
  )
  
  readCountsFiltered <- QDNAseq::estimateCorrection(
    readCountsFiltered, 
    span=0.65, 
    family="symmetric", 
    adjustIncompletes=TRUE,
    maxIter=1, 
    cutoff=4, 
    variables=c("gc", "mappability")
  )
  
  copyNumbers <- QDNAseq::correctBins(readCountsFiltered)
  
  copyNumbersNormalized <- QDNAseq::normalizeBins(
    copyNumbers,
    method="median", 
    force=FALSE)
  
  
  copyNumbersSmooth <- QDNAseq::smoothOutlierBins(
    copyNumbersNormalized,
    logTransform=TRUE, 
    force=FALSE)
  
  copyNumbersSegmented <- QDNAseq::segmentBins(
    copyNumbersSmooth, 
    transformFun="log2",
    smoothBy=FALSE, 
    alpha=1e-10, 
    undo.splits="sdundo", 
    undo.SD=1,
    force=FALSE
    # seeds=NULL
  )
  
  
  copyNumbersSegmented <- QDNAseq::normalizeSegmentedBins(
    copyNumbersSegmented,
    force = F
  )
  
  
  copyNumbersCalled <- QDNAseq::callBins(
    copyNumbersSegmented, 
    organism="human", 
    method="CGHcall",
    prior = "all"
    #cutoffs=log2(c(deletion = 0.5, loss = 1.5, gain = 2.5, amplification = 10)/2)
  )
  
  
  return(copyNumbersCalled)
}


