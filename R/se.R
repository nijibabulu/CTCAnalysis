seLongCounts <- function(se, id="Geneid", names_to="Sample", values_to="Value") {
  if(!"long_counts" %in% S4Vectors::metadata(se)) {
    long <- SummarizedExperiment::assay(se) %>%
      tibble::as_tibble(rownames = id) %>%
      tidyr::pivot_longer(-!!id, names_to = names_to, values_to = values_to)

    coldata_tbl <- SummarizedExperiment::colData(se) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      dplyr::mutate(dplyr::across(where(is.factor), as.character))
    long <- dplyr::left_join(long, coldata_tbl, by="Sample")
    S4Vectors::metadata(se)$long_counts <- long
  }
  S4Vectors::metadata(se)$long_counts
}

colTbl <- function(se, rowname_col = "Sample") {
  SummarizedExperiment::colData(se) %>% as.data.frame() %>% tibble::rownames_to_column(rowname_col) %>% tibble::as_tibble()
}

rowTbl <- function(se, include_rownames = F, rowname_col = "Geneid") {
  tbl <- SummarizedExperiment::rowData(se) %>% as.data.frame()
  if(include_rownames)
    tbl <- tbl %>%  tibble::rownames_to_column(rowname_col)
  tibble::as_tibble(tbl)
}
