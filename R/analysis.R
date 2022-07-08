gseaAnalysis <- function(deres, pathways, nperm = 1e5) {
  purrr::map(deres$tbl_results, ~{
    dplyr::select(.x, Geneid, stat) %>%
      dplyr::filter(!is.na(stat)) %>%
      tibble::deframe() %>%
      fgsea::fgseaSimple(pathways = pathways, stats = ., nperm = nperm) %>%
      dplyr::mutate(leadingEdge = purrr::map(leadingEdge, stringr::str_c, collapse = ","))
    })
}

findSampleTypeClassifier.tbl <- function(counts, meta, target, x_genes, y_genes, x_group, y_group,
                                     varInt="group", x_filter="xy", y_filter="xy", sample_summary_f=max, classifier=median, norm=computeFPKM) {
  norm_counts <- norm(counts, meta)

  norm_counts <- ensureLongAttribute(norm_counts, meta)
  long_counts <- longCounts(norm_counts)

  long_counts <- long_counts %>% dplyr::mutate(signal=log10(Value+1)) %>% group_by(Sample)

  type_map <- target %>%  dplyr::mutate_if(is.factor, as.character) %>% select(label, Type=sym(varInt))
  signals <- full_join(
    long_counts %>% filter(Geneid %in% x_genes) %>% filter(signal == max(signal)) %>% select(Sample, x_gene=Geneid, x=signal),
    long_counts %>% filter(Geneid %in% y_genes) %>% filter(signal == max(signal)) %>% select(Sample, y_gene=Geneid, y=signal),
    by="Sample"
  ) %>% full_join(type_map, by=c(Sample="label"))

  classifiers <- list(x=classifier(signals$x), y=classifier(signals$y))

  target <- target %>% dplyr::mutate_if(is.factor, as.character)
  x_signals <- target %>% filter(!!sym(varInt) %in% x_group) %>% select(label) %>% left_join(signals, by=c(label="Sample"))
  #x_samples <- target[target[,varInt] == x_group,] %>% select(label) %>% flatten_chr()
  if(stringr::str_detect(x_filter, "x")) { x_signals <-  x_signals %>% filter(x > classifiers$x) }
  if(stringr::str_detect(x_filter, "y")) { x_signals <-  x_signals %>% filter(y < classifiers$y) }

  y_signals <- target %>% filter(!!sym(varInt) %in% y_group) %>% select(label) %>% left_join(signals, by=c(label="Sample"))
  if(stringr::str_detect(y_filter, "x")) { y_signals <-  y_signals %>% filter(x < classifiers$x) }
  if(stringr::str_detect(y_filter, "y")) { y_signals <-  y_signals %>% filter(y > classifiers$y) }

  pass <- bind_rows(x_signals, y_signals) %>% select(label) %>% flatten_chr()

  list(signals=signals, classifiers=classifiers, pass=pass)
}


findSampleTypeClassifier <- function(se, x_genes, y_genes, x_group, y_group, sample_column = "Sample",
                                     varInt="Type", x_filter="xy", y_filter="xy", sample_summary_f=max,
                                     classifier=median, norm=computeFPKM, transform=log1p) {
  counts <- SummarizedExperiment::assay(se)
  norm_counts <- transform(norm(counts, rowTbl(se)))

  signals <-
    list(tibble::enframe(apply(norm_counts[x_genes,], 2, max), name="Sample", value="x"),
         tibble::enframe(apply(norm_counts[y_genes,], 2, max), name="Sample", value="y"),
         tibble::enframe(apply(norm_counts[x_genes,], 2, function(x) x_genes[which.max(x)]), name="Sample", value="x_gene"),
         tibble::enframe(apply(norm_counts[y_genes,], 2, function(x) y_genes[which.max(x)]), name="Sample", value="y_gene"),
         colTbl(se) %>% dplyr::select(dplyr::all_of(c(sample_column, varInt)))) %>%
    purrr::reduce(dplyr::inner_join, by=sample_column)

  classifiers <- list(x=classifier(signals$x), y=classifier(signals$y))

  x_signals <- colTbl(se) %>%
    dplyr::filter(!!rlang::ensym(varInt) %in% x_group) %>%
    dplyr::select(dplyr::all_of(sample_column)) %>%
    dplyr::left_join(signals, by="Sample")
  #x_samples <- target[target[,varInt] == x_group,] %>% select(label) %>% flatten_chr()
  if(stringr::str_detect(x_filter, "x")) { x_signals <-  x_signals %>% dplyr::filter(x > classifiers$x) }
  if(stringr::str_detect(x_filter, "y")) { x_signals <-  x_signals %>% dplyr::filter(y < classifiers$y) }

  y_signals <- colTbl(se) %>%
    dplyr::filter(!!rlang::ensym(varInt) %in% y_group) %>%
    dplyr::select(dplyr::all_of(sample_column)) %>%
    dplyr::left_join(signals, by="Sample")
  if(stringr::str_detect(y_filter, "x")) { y_signals <-  y_signals %>% dplyr::filter(x < classifiers$x) }
  if(stringr::str_detect(y_filter, "y")) { y_signals <-  y_signals %>% dplyr::filter(y > classifiers$y) }

  pass <- dplyr::bind_rows(x_signals, y_signals) %>% dplyr::select(dplyr::all_of(sample_column)) %>% purrr::flatten_chr()

  list(signals=signals, classifiers=classifiers, pass=pass)
}


plotSampleTypeClassfier <- function(sample_type_classifier, palette=ggthemes::ptol_pal(), shape_col=NULL, linetype=3, shapes=c(21,24), xlab="x", ylab="y", label_pass=FALSE, label_fail=FALSE) {
  df <- sample_type_classifier$signals %>% dplyr::mutate(label="")
  if(label_pass) {
    df <- df %>% dplyr::mutate(label=dplyr::if_else(Sample %in% sample_type_classifier$pass, Sample, label))
  }
  if(label_fail) {
    df <- df %>% dplyr::mutate(label=dplyr::if_else(Sample %in% sample_type_classifier$pass, label, Sample))
  }

  aesthetic <- ggplot2::aes(x=x, y=y, label=label, color=Type)
  if(!is.null(shape_col)) {
    aesthetic <- utils::modifyList(aesthetic, ggplot2::aes(shape=!!ensym(shape_col)))
  }
  ggplot2::ggplot(df, aesthetic) +
    ggplot2::geom_point() +
    ggplot2::discrete_scale("color", "group", palette = palette) +
    ggplot2::labs(x=xlab,y=ylab) +
    ggrepel::geom_text_repel(show.legend = F) +
    ggplot2::geom_hline(yintercept = sample_type_classifier$classifiers$y, linetype = linetype) +
    ggplot2::geom_vline(xintercept = sample_type_classifier$classifiers$x, linetype = linetype) +
    ggplot2::theme_bw()
}

if(0) {
findSampleTypeClassifier(prostateCounts, prostateMeta, prostateTarget, x_genes=c(epithelial_markers, prostate_gene_aliases), y_genes=leukocyte_marker_aliases,
                         x_group="ctc", y_group="fi")
}
