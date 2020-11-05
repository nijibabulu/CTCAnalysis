
findSampleTypeClassifier <- function(counts, meta, target, x_genes, y_genes, x_group, y_group,
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

plotSampleTypeClassfier <- function(sample_type_classifier, palette=ptol_pal(), shape_col=NULL, linetype=3, shapes=c(21,24), xlab="x", ylab="y", label_pass=FALSE, label_fail=FALSE) {
  df <- sample_type_classifier$signals %>% mutate(label="")
  if(label_pass) {
    df <- df %>% mutate(label=dplyr::if_else(Sample %in% sample_type_classifier$pass, Sample, label))
  }
  if(label_fail) {
    df <- df %>% mutate(label=dplyr::if_else(Sample %in% sample_type_classifier$pass, label, Sample))
  }

  aesthetic <- aes(x=x, y=y, label=label, color=Type)
  if(!is.null(shape_col)) {
    aesthetic <- utils::modifyList(aesthetic, aes(shape=!!ensym(shape_col)))
  }
  ggplot2::ggplot(df, aesthetic) +
    ggplot2::geom_point() +
    ggplot2::discrete_scale("color", "group", palette = palette) +
    labs(x=xlab,y=ylab) +
    ggrepel::geom_text_repel(show.legend = F) +
    geom_hline(yintercept = sample_type_classifier$classifiers$y, linetype = linetype) +
    geom_vline(xintercept = sample_type_classifier$classifiers$x, linetype = linetype) +
    theme_bw()
}

if(0) {
findSampleTypeClassifier(prostateCounts, prostateMeta, prostateTarget, x_genes=c(epithelial_markers, prostate_gene_aliases), y_genes=leukocyte_marker_aliases,
                         x_group="ctc", y_group="fi")
}
