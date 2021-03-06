enhanceDESeq2ResultTbl <- function(tbl, topn=20, fc=2, alpha=0.05, labels=NULL, filter_na = TRUE, rank=NULL) {
  if(is.null(rank)) {
    data <- tbl %>% dplyr::mutate(rank=sqrt(log2FoldChange**2 + (-log2(padj))**2))
  } else {
    data <- tbl %>% dplyr::mutate(rank=rank)
  }

  data <- data %>%
    dplyr::mutate(filtered=is.na(padj),
                  sig=dplyr::case_when(
                    abs(log2FoldChange) > log2(fc) & padj <= alpha ~ "s",
                    abs(log2FoldChange) > log2(fc) ~ "e" ,
                    TRUE ~ "n"
                  ),
                  ) %>%
    tidyr::replace_na(list(padj=1, sig="n")) %>%
    dplyr::filter(!is.na(log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(rank))


  if(is.null(labels) && topn > 0) {
    labels <- dplyr::bind_rows(
      data %>% dplyr::filter(log2FoldChange > 0 & sig == "s") %>% dplyr::slice(1:topn),
      data %>% dplyr::filter(log2FoldChange < 0 & sig == "s") %>% dplyr::slice(1:topn)
      )  %>%
      dplyr::select(Geneid) %>% purrr::flatten_chr()
  }

  if(filter_na) {
    data <- data %>% dplyr::filter(filtered == FALSE)
  }

  data %>% dplyr::mutate(label=dplyr::if_else(Geneid %in% labels, Geneid, ""),
                         Status=translated_factor(sig, c(s="p < 0.05, fc > 2", e="fc > 2", n="NS")))
}

DESeq2VolcanoPlot <- function(tbl, samples, name="Volcano Plot", alpha=0.05, fc=2, labels=NULL, topn=20, show_filtered = FALSE, rank=NULL, max.overlaps=100) {
  data <- enhanceDESeq2ResultTbl(tbl, topn = topn, labels = labels, filter_na = !show_filtered, rank=rank, fc = fc, alpha = alpha)

  aesthetic <- ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=Status, label=label)
  if(show_filtered) {
    aesthetic <- utils::modifyList(aesthetic, shape=filtered)
  }

  nplussig <- data %>% dplyr::filter(sig=="s" & log2FoldChange > 0) %>% nrow()
  nminussig <- data %>% dplyr::filter(sig=="s" & log2FoldChange < 0) %>% nrow()

  p <- ggplot2::ggplot(data, aesthetic) +
    ggplot2::annotate("rect", xmin = 0, ymin = 0, xmax = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_color_manual(values=purrr::set_names(c("red", "goldenrod", "grey"), levels(data$Status))) +
    ggrepel::geom_text_repel(show.legend = F, color = "black", max.overlaps = max.overlaps) +
    ggplot2::labs(title = stringr::str_glue("{samples[1]} versus {samples[2]} Volcano Plot"),
                  legend = "") +
    ggplot2::annotate("text", x = 1, y = Inf, vjust = "inward", hjust = 0, label = samples[2]) +
    ggplot2::annotate("text", x = -1, y = Inf, vjust = "inward", hjust = 1, label = samples[1]) +
    ggplot2::annotate("text", x = Inf, y = Inf, vjust = "inward", hjust = "inward", label = stringr::str_glue("n = {nplussig}")) +
    ggplot2::annotate("text", x = -Inf, y = Inf, vjust = "inward", hjust = "inward", label = stringr::str_glue("n = {nminussig}")) +
    ggplot2::theme_bw()

  if(show_filtered) {
    p <- p + ggplot2::scale_shape_manual(values=c(20,6))
  }

  p
}

DESeq2VolcanoPlots <- function(deres, alpha=0.05, fc=2, labels=NULL, topn=20, show_filtered = FALSE, rank = NULL) {
  purrr::map2(deres$tbl_results, stringr::str_split(names(deres$tbl_results), "_vs_"), DESeq2VolcanoPlot, topn=topn, alpha=alpha, fc=fc, labels=labels)
}

DESeq2MAPlot <- function(tbl, samples, alpha=0.05, fc=2, labels=NULL, topn=20, show_filtered = FALSE, rank = NULL) {
  data <- enhanceDESeq2ResultTbl(tbl, topn = topn, labels = labels, filter_na = !show_filtered, rank = rank)

  aesthetic <- aes(y=log2FoldChange, x=baseMean, color=Status, label=label)
  if(show_filtered) {
    aesthetic <- utils::modifyList(aesthetic, shape=filtered)
  }

  p <- ggplot(data, aesthetic) +
    geom_point(size = 0.5) +
    scale_color_manual(values=c("red", "goldenrod", "grey")) +
    geom_text_repel(show.legend = F, color = "black") +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(title = stringr::str_glue("{samples[1]} versus {samples[2]} MA Plot"),
         x = stringr::str_glue("Mean of {samples[1]}")) +
    theme_bw()

  if(show_filtered) {
    p <- p + ggplot2::scale_shape_manual(values=c(20,6))
  }

  p
}

DESeq2MAPlots <- function(deres, alpha=0.05, fc=2, labels=NULL, topn=20, show_filtered = FALSE, rank = NULL) {
  purrr::map2(deres$tbl_results, stringr::str_split(names(deres$tbl_results), "_vs_"), DESeq2MAPlot, alpha=alpha, fc=fc, labels=labels, rank = rank)
}

exprHeatmap <- function(se, genes = NULL, palette = ggthemes::stata_pal(), varInt="Type", ann_title = "Type",
                        result_tbl=NULL, topn = 20, rank = NULL, norm = computeTPM, transform = log1p,
                        scale=T, center=T, unplottable_action = c("remove", "zero"), ...) {
  unplottable_action = match.arg(unplottable_action)
  if(is.null(genes) & is.null(result_tbl)) {
    stop("Need either a subset of genes from the expression set or a results tbl")
  }

  #transform the matrix
  counts <- SummarizedExperiment::assay(se)
  if(!is.null(norm))
    counts <- norm(counts, rowTbl(se))
  if(!is.null(transform))
    counts <- transform(counts)

  if(scale) {
    counts <- t(scale(t(counts), center = center))
  }

  # subset the data
  if(!is.null(result_tbl)) {
    data <-  enhanceDESeq2ResultTbl(result_tbl, topn = topn, rank = rank)
    genes <- purrr::discard(data$label, ~stringr::str_length(.) == 0)
  }
  counts <- counts[genes, ]

  # handle values we cannot plot
  bad_rows <- apply(counts, 1, function(x) any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x)))
  if(sum(bad_rows))
    warning(stringr::str_c(rownames(counts)[bad_rows], collapse = ", "), " contain unplottable values")
  if(unplottable_action == "remove") { counts <- counts[!bad_rows, ] }
  else if(unplottable_action == "zero") { counts[bad_rows, ]  <- 0 }

  ann_data <- colTbl(se) %>% dplyr::select(dplyr::all_of(varInt))
  values <- purrr::map(ann_data, purrr::compose(unique, as.character)) %>% purrr::flatten_chr()
  colors <- palette(length(values))
  col_map <- list(purrr::set_names(colors, values))
  col_ann <- ComplexHeatmap::HeatmapAnnotation(df = ann_data %>% dplyr::mutate_if(is.factor, as.character) %>% as.data.frame(),
                                               which = "column",
                                               col = purrr::set_names(rep(col_map, length(varInt)), varInt))


  lim <- max(abs(quantile(counts, 0.02, na.rm = T)), abs(quantile(counts, 0.98, na.rm = T)))


  if(scale) {
    heat_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("Darkblue", "white", "red"))
  } else {
    cols <- c("white", "#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C", "purple", "black")
    heat_fun <- circlize::colorRamp2(c(0,.Machine$double.xmin, lim/rev(seq_along(cols[-c(1,2)]))), cols)
  }
  lgd <- list( title = "Expr.")


  grid::grid.grabExpr(ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(counts, col = heat_fun, top_annotation = col_ann,  heatmap_legend_param = lgd, ...)
  ))
}
