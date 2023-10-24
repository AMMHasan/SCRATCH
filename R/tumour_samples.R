#' tumour samples
#'
#' @param BAM_path character string for BAM directory
#' @param pattern unique pattern tumour samples
#'
#' @return character vector
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' pattern = c(paste0("PEA310","_1"))
#' BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM_subset/"
#' tumour_samples(BAM_path, pattern)

tumour_samples <- function(BAM_path, pattern) {
  name_tumour_samples <- list.files(BAM_path, pattern) %>% 
    gsub(".bam.bai", "", x=.) %>% 
    gsub(".bam", "", x=.)  %>% 
    unique()
  return(name_tumour_samples)
}
