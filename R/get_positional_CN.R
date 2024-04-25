#' get copy number by provided position
#'
#' @param BAM_path path of the .bam files
#' @param bin_size bin size to calculate copy number by QDNAseq
#' @param pattern pattern for the project samples
#' @param copyNumbersCalled_obj a QDNAseq object
#' @param chr chromosome: single or multiple in a vector
#' @param pos position (in bp): single or multiple in a vector (according to the chromosome)
#' @importFrom magrittr %>%
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
#' @import splitstackshape
#' @return a data frame
#' @export
#' @examples
#' BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM_subset/"
#' bin_size <- 500
#' pattern <- "PEA310"
#' copyNumbersCalled_obj <- generate_copyNumberCalled_obj(BAM_path, pattern, bin_size)
#' chr <- c(1, 2, 13)
#' pos <- c(1000301, 14000000, 60000000)
#' get_positional_CN(BAM_path, pattern, bin_size, copyNumbersCalled_obj, chr, pos)
#'
get_positional_CN <- function(BAM_path, pattern, bin_size, copyNumbersCalled_obj, chr, pos) {
  name_tumour_samples <- tumour_samples(BAM_path, pattern)

  pos_CN_list <- list()

  for (i in name_tumour_samples) {
    test_df <- ACE::getadjustedsegments(ACE::objectsampletotemplate(copyNumbersCalled_obj, index = i),
      cellularity = ACE::squaremodel((ACE::objectsampletotemplate(copyNumbersCalled_obj, index = i)))[["minimadf"]]$cellularity[1],
      ploidy = ACE::squaremodel((ACE::objectsampletotemplate(copyNumbersCalled_obj, index = i)))[["minimadf"]]$ploidy[1]
    )

    pos_CN_df <- ACE::analyzegenomiclocations(test_df, chr, pos) %>%
      dplyr::mutate(sampleID = i)

    pos_CN_list[[i]] <- pos_CN_df
  }



  return(Reduce(rbind, pos_CN_list))
}
