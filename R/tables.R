
#' Create a workbook sheet out of a data frame
#'
#' This function is similar in spirit to \code{\link{xlsx::addDataFrame()}} but
#' has customized behavior to accommodate our standard formats.
#'
#' @param wb an excel workbook created with \code{\link{xlsx::createWorkbook()}}
#' @param df a data frame to place in the workbook
#' @param name the name of the sheet
#' @param auto_size_columns a vector of column names or indices which should be
#'   made auto-sized (i.e. the width of the column will equal the widest entry)
#' @param emph_style a style created with \code{\link{xlsx::CellStyle()}} which
#'   will be used for the emphasized data points (see \code{...} below for
#'   details).
#' @param head_style a style created with \code{\link{xlsx::CellStyle()}} which
#'   will be used for the column headers. Default is to add a bottom border.
#' @param ... predicates (functions which return \code{TRUE} or \code{FALSE})
#'   which will be used to determine if the cell in the column will be given
#'   the \code{emph_style} cell style. The names are used to determine which
#'   column to apply this predicate to.
#'
#' @return the \code{wb} with a new sheet added
#' @export
#'
#' @examples
#' \dontrun{
#'  wb <- xlsx::createWorkbook()
#'  table_add_df_sheet(wb, mtcars, "mtcars", am = ~.x == 1)
#'  mtcars$car <- rownames(mtcars)
#'  table_add_df_sheet(wb, mtcars, "mtcars_widened", auto_size_columns = c("car", "am", "vs"))
#' }
table_add_df_sheet <- function(wb, df, name, auto_size_columns = c(),
                               emph_style = xlsx::CellStyle(wb) + xlsx::Fill(foregroundColor = "#92D050"),
                               head_style = xlsx::CellStyle(wb) + xlsx::Border(position = "BOTTOM"),
                               ...) {
  emph_funcs <- list(...)
  df <- as.data.frame(df)
  row.names(df) <- NULL
  s <- xlsx::createSheet(wb, name)
  h <- xlsx::createRow(s, 1)
  hs <- xlsx::createCell(h, 1:ncol(df))
  purrr::walk2(hs, colnames(df), xlsx::setCellValue)
  purrr::walk(hs, xlsx::setCellStyle, head_style)
  rs <- xlsx::createRow(s, 2:(nrow(df)+1))
  cs <- xlsx::createCell(rs, 1:ncol(df))
  purrr::iwalk(colnames(df),
               ~{
                 purrr::walk2(cs[,.y], df[,.x, drop = T], xlsx::setCellValue)
                 if(.x %in% names(emph_funcs)) {
                   is <- purrr::map_lgl(df[,.x,drop = T], emph_funcs[[.x]])
                   emph_cs <- cs[which(is),.y]
                   purrr::walk(emph_cs, xlsx::setCellStyle, emph_style)
                 }
               })
  if(!is.numeric(auto_size_columns))
    auto_size_columns <- which(colnames(df) %in% auto_size_columns)
  purrr::walk(auto_size_columns, xlsx::autoSizeColumn, sheet = s)
  wb
}

#' Add a notes sheet to an excel workbook
#'
#' @param wb an excel workbook created with \code{\link{xlsx::createWorkbook()}}
#' @param notes a character vector of notes. Empty strings become spacer rows
#'   in the cell
#' @param name the name of the sheet.
#' @param width the column width of the notes column
#'
#' @return
#'   The workbook with a notes column added.
#' @export
#'
#' @examples
#' \dontrun{
#'  wb <- xlsx::createWorkbook()
#'  table_add_notes(wb, c("this is one note",
#'                        "this is another note",
#'                        "", # spacer row
#'                        "this is a summarizing note."))
#' }
table_add_notes <- function(wb, notes, name = "Notes", width = 100) {
  s <- xlsx::createSheet(wb, name)
  rs <- xlsx::createRow(s, seq_along(notes))
  cs <- xlsx::createCell(rs, 1)
  purrr::walk2(cs, notes, xlsx::setCellValue)
  xlsx::setColumnWidth(s, 1, width)
  wb
}



table_output_p116a_excel <- function(pc_bor, pc_mts_res, pc_surv_res, pc_mkr_pf, out) {
  pc_bor_out <-
    pc_bor %>%
    dplyr::filter(VISITNUM == 1.01) %>%
    dplyr::select(marker = PCTESTCD,
                  visit = VISITNUM,
                  n.bor.y = BORY_Num,
                  med.bor.y = BORY_median,
                  n.bor.n = BORN_Num,
                  med.bor.n = BORN_median,
                  p.wilcox = wilcox_p.value,
                  p.ttest = ttest_p.value,
                  p.reg = uv_binom_p.value,
                  p.reg.mv = mv_binom_p.value,
                  odds.ratio = mv_binom_estimate,
                  fdr.wilcox = fdr_wilcox,
                  fdr.ttest = fdr_ttest,
                  fdr.reg = uv_binom_fdr,
                  fdr.reg.mv = mv_binom_fdr) %>%
    dplyr::mutate(visit = "Baseline")
  # TODO get fdr for reg, reg.mv
  pc_mts_out <-
    pc_mts_res %>%
    dplyr::filter(VISITNUM == 1.01) %>%
    dplyr::select(marker = PCTESTCD,
                  visit = VISITNUM,
                  n = subjects_N,
                  cor.pearson = pearson_corr,
                  p.pearson = pearson_p,
                  cor.spearman = spearman_corr,
                  p.spearman = spearman_p,
                  p.reg = uv_p.value,
                  p.reg.mv = mv_p.value,
                  estimate = mv_estimate,
                  fdr.pearson = pearson_fdr,
                  fdr.spearman = spearman_fdr,
                  fdr.reg = uv_fdr,
                  fdr.reg.mv = mv_fdr) %>%
    dplyr::mutate(visit = "BASE")
  pc_surv_out <-
    pc_surv_res %>%
    dplyr::ungroup() %>%
    dplyr::filter(outcome == "PFS", cutoff_type == "median") %>%
    dplyr::select(Marker = PCTESTCD,
                  visit = VISIT,
                  n.wald = uv_n,
                  p.wald = uv_p.wald,
                  hr = uv_hr,
                  n.mv = mv_n,
                  p.reg.mv = mv_p.marker,
                  hr.adj = mv_hr,
                  fdr.wald = uv_fdr.wald,
                  fdr.reg.mv = mv_fdr.marker) %>%
    dplyr::mutate(visit = "BASE")
  # TODO get fdr for wald, reg.mv
  pc_mkr_pf_out <-
    pc_mkr_pf %>%
    tibble::rowid_to_column("SEQ") %>%
    dplyr::select(marker = PCTESTCD,
                  var = term,
                  Estimate = estimate,
                  `PR(>|t|)` = p.value,
                  SEQ) %>%
    dplyr::filter(var != "(Intercept)")
  notes <- c("In each sheet, p values of univariate regression (p.reg), and multivariable regression (p.reg.mv) are shown; odds ratio (BOR), estimate of coefficient (MTS), hazard ratio (PFS) are also shown",
             "Results of previous analysis (Wilcoxon, t-test for BOR; Spearman, Pearson for MTS; univariate cox regression for PFS) are shown for comparison",
             "p-values < 0.1 are colored",
             "",
             "The last sheet \"marker vs. var\" includes the results of evaluting what demographic and clinical variables impacted baseline serum marker levels (based on multivariable regression)"
  )
  pv_emph <- ~ .x < 0.05
  wb <- xlsx::createWorkbook()
  wb <- wb %>%
    table_add_notes(notes) %>%
    table_add_df_sheet(pc_bor_out, "BOR",
                       auto_size_columns = "marker",
                       p.wilcox = pv_emph, p.ttest = pv_emph,
                       p.reg = pv_emph, p.reg.mv = pv_emph) %>%
    table_add_df_sheet(pc_mts_out, "MTS",
                       auto_size_columns = "marker",
                       p.pearson = pv_emph, p.spearman = pv_emph,
                       p.reg = pv_emph, p.reg.mv = pv_emph) %>%
    table_add_df_sheet(pc_surv_out, "PFS",
                       auto_size_columns = "marker",
                       p.wald = pv_emph, p.reg.mv = pv_emph) %>%
    table_add_df_sheet(pc_mkr_pf_out, "marker vs. var",
                       auto_size_columns = c("marker", "var"),
                       emph_style = xlsx::CellStyle(wb) + xlsx::Fill(foregroundColor = "#D6E0F2"),
                       `PR(>|t|)` = pv_emph)

  xlsx::saveWorkbook(wb, out)

}
