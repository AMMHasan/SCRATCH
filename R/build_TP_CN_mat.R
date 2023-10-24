#' create copy number matrix for transition points
#'
#' @param list_TP_CN list of data frame with TP position in the first column and the copy number of the TP in the second column
#'
#' @return a matrix
#' @export
#'
#' @examples
#' df1 <- data.frame(TP=c("TP1","TP2","TP3","TP4"), CN=c(1,2,3,4))
#' df2 <- data.frame(TP=c("TP2","TP3","TP4"), CN=c(1,5,6))
#' df3 <- data.frame(TP=c("TP1","TP5","TP3","TP7"), CN=c(11,22,35,14))
#' df_list <- list(df1,df2,df3)
#' names(df_list) <- c("df1","df2","df3")
#' build_TP_CN_mat(df_list)
#' 
#' 
build_TP_CN_mat <- function(list_TP_CN) {
  unique_elements <- sort(unique(unlist(lapply(list_TP_CN, function(df) df[,1]))))
  my_matrix <- matrix(NA, 
                   nrow = length(unique_elements), 
                   ncol = length(list_TP_CN), 
                   dimnames = list(unique_elements, names(list_TP_CN)))
  for (i in names(list_TP_CN)) {
    for(j in list_TP_CN[[i]][,1]){
      my_matrix[j,i] <- list_TP_CN[[i]][,2][list_TP_CN[[i]][,1] == j]
    }
  }
  
  my_matrix[is.na(my_matrix)] <- 0
  
  return(my_matrix)
  
}


