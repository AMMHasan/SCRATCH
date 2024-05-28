#' create Transition Point binary matrics
#'
#' @param BAM_path path of the .bam files
##' @param bin_size bin size to calculate copy number by QDNAseq
#' @param pattern pattern for the project samples
#' @param copyNumbersCalled_obj a QDNAseq object
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
#' @return a binary matrix
#' @export
#' @examples
#' BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM_subset/"
#' bin_size <- 500
#' pattern <- "PEA310"
#' copyNumbersCalled_obj <- generate_copyNumberCalled_obj(BAM_path, pattern, bin_size)
#' generate_TP_CN_matrix(BAM_path, pattern, bin_size, copyNumbersCalled_obj)
#'
generate_TP_CN_matrix <- function(BAM_path, pattern, copyNumbersCalled_obj) {
  name_tumour_samples <- tumour_samples(BAM_path, pattern)

  TP_list <- list()

  for (i in name_tumour_samples) {
    test_df <- ACE::getadjustedsegments(ACE::objectsampletotemplate(copyNumbersCalled_obj, index = i),
      cellularity = ACE::squaremodel((ACE::objectsampletotemplate(copyNumbersCalled_obj, index = i)))[["minimadf"]]$cellularity[1],
      ploidy = ACE::squaremodel((ACE::objectsampletotemplate(copyNumbersCalled_obj, index = i)))[["minimadf"]]$ploidy[1]
    ) %>%
      dplyr::select(Chromosome, Start, End, Segment_Mean, Copies) %>%
      dplyr::mutate(TP_S = paste0(Chromosome, ":", Start), TP_E = paste0(Chromosome, ":", End))
    TP_list[[i]] <- rbind(
      test_df %>% dplyr::select(TP = TP_S, CN = Copies),
      test_df %>% dplyr::select(TP = TP_E, CN = Copies)
    )
  }

  TP_CN_matrix <- build_TP_CN_mat(TP_list)

  return(TP_CN_matrix)
}
