read_featureCounts <- function(filename) {
  readr::read_tsv(filename,
                  comment =  "#",
                  col_types = readr::cols(
                    .default = readr::col_double(),
                    Geneid = readr::col_character(),
                    Chr = readr::col_character(),
                    Start = readr::col_character(),
                    End = readr::col_character(),
                    Strand = readr::col_character()
                  ))
}

read_featureCountsSummary <- function(filename, pivot = T) {
   summary_df <- readr::read_tsv(filename,
                                 col_types = readr::cols(
                                   .default = readr::col_double(),
                                   Status = readr::col_character())
                                 )
   if(pivot) {
     summary_df <- summary_df %>% tidyr::pivot_longer(-Status) %>% tidyr::pivot_wider(name, Status)
   }
   summary_df
}

merge_featureCounts <- function(..., verbose = F) {
  files <- purrr::flatten(list(...))
  if(!is.null(purrr::detect(files, purrr::negate(is.character)))) {
    stop(stringr::str_glue("One of merge_featureCounts files is not a filename: {stringr::str_c(files, collapse = ', ')}"))
  }

  fc_dfs <- purrr::map(files, read_featureCounts)
  combined_fc <- purrr::reduce(
    fc_dfs,  ~{
      if(verbose && any(.x[,1:6] != .y[,1:6])) {
        warning("Merging featureCounts with different gene columns")
      }
      dplyr::full_join(.x, .y, by=c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  })

  combined_fc
}

merge_featureCountsSummaries <- function(..., verbose = F) {
  files <- purrr::flatten(list(...))
  if(!is.null(purrr::detect(files, purrr::negate(is.character)))) {
    stop(stringr::str_glue("One of merge_featureCountsSummaries files is not a filename: {stringr::str_c(files, collapse = ', ')}"))
  }

  fcs_dfs <- purrr::map(files, read_featureCountsSummary)
  combined_fcs <- purrr::reduce(fcs_dfs, ~{
    if(verbose && (ncol(.x) != ncol(.y) || any(sort(colnames(.x)) != any(sort(colnames(.y)))))) {
      warning("Merging featureCountsSummaries with different gene columns")
    }
    dplyr::bind_rows(.x,.y)
  })

  combined_fcs
}

merge_featureCountsFiles <- function(..., target_base = NULL, merge_summaries = T, verbose = F) {
  if(is.null(target_base)) stop("Required target_base is missing")
  merged_fc <- merge_featureCounts(...)
  readr::write_tsv(merged_fc, file = stringr::str_glue("{target_base}"))
  if(merge_summaries) {
    fcs_files <- purrr::flatten(list(...)) %>% stringr::str_c("summary", sep = ".")
    merged_fcs <- merge_featureCountsSummaries(fcs_files)

    # merged summaries are pivoted, need to re-pivot to the file format
    merged_fcs %>%
      tidyr::pivot_longer(-name, names_to = "Status") %>%
      tidyr::pivot_wider(Status, name) %>%
      readr::write_tsv(stringr::str_glue("{target_base}.summary"))
  }
}

add_col_data <- function(se, tbl, join_col="Patient", all.y=F, all.x = T, ...) {
  coldf <- SummarizedExperiment::colData(se)
  coldf$old.row.names <- rownames(coldf)
  merged <- dplyr::left_join(as.data.frame(coldf), tbl, by=join_col)
  rownames(merged) <- merged$old.row.names
  merged$old.row.names <- NULL
  SummarizedExperiment::colData(se) <- as(merged, "DataFrame")
  se
}

alter_col_tbl <- function(se, f, ...)   {
  SummarizedExperiment::colData(se) <-
    SummarizedExperiment::colData(se) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    f(...) %>%
    tibble::column_to_rownames("Sample") %>%
    as("DataFrame")
  se
}
