#' create copy number matrix for transition points
#'
#' @param TP_CN_mat a matrix of copy-number at the transition points
#'
#' @importFrom magrittr %>%
#' @return a tree object
#' @export
#'
#' @examples
#' df1 <- data.frame(TP = c("TP1", "TP2", "TP3", "TP4"), CN = c(1, 2, 3, 4))
#' df2 <- data.frame(TP = c("TP2", "TP3", "TP4"), CN = c(1, 5, 6))
#' df3 <- data.frame(TP = c("TP1", "TP5", "TP3", "TP7"), CN = c(11, 22, 35, 14))
#' df_list <- list(df1, df2, df3)
#' names(df_list) <- c("df1", "df2", "df3")
#' TP_CN_mat <- build_TP_CN_mat(df_list)
#' build_tree(TP_CN_mat)
#'
build_tree <- function(TP_CN_mat) {
  met_tree <- TP_CN_mat %>%
    stats::cor() %>%
    tibble::as_tibble() %>%
    mutate_all(function(x){1-x}) %>% 
    stats::as.dist() %>%
    phangorn::upgma()

  return(met_tree)
}
