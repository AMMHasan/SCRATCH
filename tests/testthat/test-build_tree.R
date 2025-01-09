test_that("returns a phylo data class output", {
  # expect_equal(2 * 2, 4)
  df1 <- data.frame(TP = c("TP1", "TP2", "TP3", "TP4"), CN = c(1, 2, 3, 4))
  df2 <- data.frame(TP = c("TP2", "TP3", "TP4"), CN = c(1, 5, 6))
  df3 <- data.frame(TP = c("TP1", "TP5", "TP3", "TP7"), CN = c(11, 22, 35, 14))
  df_list <- list(df1, df2, df3)
  names(df_list) <- c("df1", "df2", "df3")
  TP_CN_mat <- build_TP_CN_mat(df_list)
  
  expect_identical(build_tree(TP_CN_mat) %>% class(), "phylo")
})
