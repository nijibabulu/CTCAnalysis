computeTPM <- function(mat, meta) {
  lengths_m <- meta %>% dplyr::select(Length) %>% purrr::flatten_dbl()
  tpl <- mat/lengths_m
  tpm <- 1e6*t(t(tpl)/colSums(tpl))
  tpm
}

computeFPKM <- function(mat, meta) {
  lengths_m <- meta %>% dplyr::select(Length) %>% purrr::flatten_dbl()
  10**9*t(t(mat/lengths_m)/colSums(mat))
}

computeDESeqNormalizedValues <- function(mat, meta) {
  sf <- DESeq2::estimateSizeFactorsForMatrix(mat)
  t(t(mat)/sf)
}

computeCov <- function(mat, meta) {
  lengths_m <- meta %>% dplyr::select(Length) %>% purrr::flatten_dbl()
  mat/lengths_m
}
