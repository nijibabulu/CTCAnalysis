summaryPlot <- function(counts, x, y, geom=ggplot2::geom_col(),
                        meta=NULL, group="group", facet=NULL,
                        facet_scales = "fixed", facet_space = "fixed",
                        xlab = "Sample",
                        palette=ggthemes::ptol_pal()) {

  summaries <- geneSummaries(counts, meta)
  aesthetic <- aes(x=!!ensym(x), y=!!ensym(y))
  if(!is.null(group))
    aesthetic <- utils::modifyList(aesthetic, aes(fill=!!ensym(group)))

  p <- ggplot2::ggplot(summaries, aesthetic) +
    geom +
    ggplot2::discrete_scale("fill", "group", palette = palette) +
    ggplot2::labs(x = xlab) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90),
                   legend.title = ggplot2::element_blank())

  if(!is.null(facet)) {
    p <- p + ggplot2::facet_grid(facet, scales = facet_scales, space=facet_space)
  }
  p

}


barplotNull <- function(counts, meta, ...) {
  summaryPlot(counts, meta = meta, "Sample", "NumberOfNulls", ...)
}
barplotGeneCounts <- function(counts, meta, ...) {
  summaryPlot(counts, meta = meta, "Sample", "NumberOfGenes", ...)
}

boxplotNull <- function(counts, meta, geom=ggplot2::geom_boxplot(), ...) {
  summaryPlot(counts, meta= meta, "group", "NumberOfNulls", geom = geom, xlab = "Type", ...) +
    ggplot2::theme(legend.position = "none")
}
boxplotGeneCounts <- function(counts, meta, geom=ggplot2::geom_boxplot(), ...) {
  summaryPlot(counts, meta = meta, "group", "NumberOfGenes", geom=geom, xlab = "Type",  ...) +
    ggplot2::theme(legend.position = "none")
}

maxGenePlot <- function(counts, meta=NULL, norm=NULL, group="group",
                        facet=NULL, facet_scales="fixed", facet_space="fixed",
                        palette=ggthemes::ptol_pal()) {
  maxGeneTbl <- maxGene(counts, meta, norm)
  aesthetic <- aes(x=Sample, y = Value, label=Geneid)
  if(!is.null(group))
    aesthetic <- utils::modifyList(aesthetic, aes(fill=!!ensym(group)))
  p <- ggplot2::ggplot(maxGeneTbl, aesthetic) +
    geom_col() +
    ggplot2::discrete_scale("fill", "group", palette = palette) +
    ggplot2::geom_text(angle=90) +
    ggplot2::labs(x = "Sample", y = "Frequency", fill = "Type") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))

  p
}

genePCAPlot <- function(counts, meta=NULL, norm=NULL, group="group", palette=ggthemes::ptol_pal(),
                        pcs=3, nrow=1, n_genes = 500, label_points = F) {
  pca <- genePCA(counts, meta)
  aesthetic <- aes(x=Sample, y = Value, label=Geneid)
  if(is.numeric(pcs)) {
    pc_names <- stringr::str_glue("PC{seq(pcs)}")
    pcs <- purrr::cross2(pc_names, pc_names, .filter = `>=`)
  }
  ps <- pcs %>%
    purrr::map(function(pcs) {
      pc1 <- pcs[[1]]; pc2 <- pcs[[2]]
      aesthetic <- aes(x=!!ensym(pc1), y=!!ensym(pc2), label=Sample)
      if(!is.null(group))
        aesthetic <- utils::modifyList(aesthetic, aes(color=!!ensym(group)))
      p <- ggplot2::ggplot(pca$tbl, aesthetic) +
        ggplot2::geom_point() +
        discrete_scale("color", "group", palette = palette) +
        labs(x=str_glue("{pc1} ({pca$prp[pc1]} % Var.)"),
             y=str_glue("{pc2} ({pca$prp[pc2]} % Var.)"),
             color="Type")

      if(label_points) {
        p <- p + ggrepel::geom_text_repel(show.legend = F)
      }
      pc_vals <- pca$tbl %>% dplyr::select(!!!pc1, !!!pc2) %>% purrr::flatten_dbl()
      limits <- c(min(pc_vals), max(pc_vals))
      p + coord_fixed(xlim=limits, ylim=limits) + theme_bw()
    })
  if(length(ps) > 1) {
    patchwork::wrap_plots(ps, patchwork::guide_area(), nrow = nrow, widths = 1, heights = 1, guides = "collect")
  } else {
    ps[[1]]
  }
}

sampleClustPlot <- function(counts, meta=NULL, title = "Cluster Dendrogram", method="ward.D",
                            vst = DESeq2::varianceStabilizingTransformation) {
  hc <- sampleHclust(counts, meta, vst = vst, clust_method = method)

  data <- ggdendro::dendro_data(hc, type = "rectangle")
  ggplot2::ggplot(ggdendro::segment(data), aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_segment() +
    ggplot2::scale_x_continuous(breaks = seq_along(data$labels$label),
                                labels = data$labels$label) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, y = "Distance", x = "Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())

}

