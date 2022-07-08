summaryPlot <- function(se, x, y, geom=ggplot2::geom_col(),
                        group="Type", facet=NULL,
                        facet_scales = "fixed", facet_space = "fixed",
                        xlab = "Sample", color="black",
                        palette=ggthemes::ptol_pal()) {

  summaries <- geneSummaries(se)
  aesthetic <- ggplot2::aes(x=!!rlang::ensym(x), y=!!rlang::ensym(y))
  if(!is.null(group))
    aesthetic <- utils::modifyList(aesthetic, ggplot2::aes(fill=!!rlang::ensym(group)))

  p <- ggplot2::ggplot(summaries, aesthetic, color=color) +
    geom +
    ggplot2::discrete_scale("fill", "group", palette = palette) +
    ggplot2::labs(x = xlab) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90),
                   legend.title = ggplot2::element_blank())

  if(!is.null(facet)) {
    p <- p + ggplot2::facet_grid(facet, scales = facet_scales, space=facet_space)
  }
  p

}

barplotNull <- function(se,  ...) {
  summaryPlot(se,  "Sample", "NumberOfNulls", ...)
}
barplotGeneCounts <- function(se,  ...) {
  summaryPlot(se,  "Sample", "NumberOfGenes", ...)
}

boxplotNull <- function(counts, meta, geom=ggplot2::geom_boxplot(), ...) {
  summaryPlot(counts, "Type", "NumberOfNulls", geom = geom, xlab = "Type", ...) +
    ggplot2::theme(legend.position = "none")
}
boxplotGeneCounts <- function(counts, meta, geom=ggplot2::geom_boxplot(), ...) {
  summaryPlot(counts, "Type", "NumberOfGenes", geom=geom, xlab = "Type",  ...) +
    ggplot2::theme(legend.position = "none")
}


plot_summary_stats.tbl_df <- function(tbl, x="Type", facet="Stage~Patient",
                               levels=c("Unmapped", "NoFeatures", "Ambiguity", "Multimapped", "Assigned"),
                               colors=ggthemes::ptol_pal()(length(levels)),
                               facet_scales = "fixed", facet_space = "fixed") {
  p <- ggplot2::ggplot(tbl, ggplot2::aes(x=.data[[x]], y = value, fill=factor(name, levels=levels))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::theme_bw() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
          axis.title.y = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.title = ggplot2::element_blank())
  if(!rlang::is_null(facet)) {
    p <- p + ggplot2::facet_grid(stringr::str_glue("{facet}"),  scales = facet_scales, space=facet_space)
  }
  p
}

plot_summary_stats <- function(se, x="Type", facet="Stage~Patient",
                               levels=c("Unmapped", "NoFeatures", "Ambiguity", "MultiMapping", "Assigned"),
                               colors=ggthemes::ptol_pal()(length(levels)),
                               facet_scales = "fixed", facet_space = "fixed") {
  tbl <- SummarizedExperiment::colData(se) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::rename_with( ~stringr::str_remove(.x, "Unassigned_"), starts_with("Unassigned_")) %>%
    tidyr::pivot_longer(levels)
  p <- ggplot2::ggplot(tbl, ggplot2::aes(x=!!rlang::ensym(x), y = value, fill=factor(name, levels=levels))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::theme_bw() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
  if(!rlang::is_null(facet)) {
    p <- p + ggplot2::facet_grid(stringr::str_glue("{facet}"),  scales = facet_scales, space=facet_space)
  }
  p
}

gseaPathwayPlot <- function(pathway_name, gsea_results, pathways, stats, gseaParam = 1, ticksSize = 0.2, linecolor="red",
                            label_stats=c("NES", `p-value`="pval", `q-value`="padj"), label_digits = 2,
                            title=pathway_name, titleSize=NULL, titleHjust=0.5) {

  pathway <- pathways[[pathway_name]]
  rnk <- rank(-stats)
  ord <- order(rnk)

  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0), label = "", stringsAsFactors = F)
  es.x = which.max(toPlot$y)
  es = toPlot[es.x,]$y
  gsea_result <- gsea_results[gsea_results$pathway == pathway_name,]
  toPlot[es.x,'label'] = stringr::str_glue("ES={prettyNum(es, digits={label_digits})}")


  diff <- (max(tops) - min(bottoms)) / 8

  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot2::ggplot(toPlot, ggplot2::aes(x=x, y=y))


  g <- g + ggplot2::geom_point(color=linecolor, size=0.1)

  if(! is.null(label_stats)) {
    label_stats_names <- dplyr::if_else(stringr::str_length(names(label_stats)) == 0, label_stats, names(label_stats))
    stats_info <- purrr::map2_chr(label_stats, label_stats_names,
                                  ~stringr::str_glue("{.y}={prettyNum(gsea_result[.x], digits={label_digits})}")) %>%
      stringr::str_c(collapse = "\n")
    g <- g + ggplot2::annotate("text", x = length(stats), y = es, vjust = 1, hjust = 1, label=stats_info)
             #label = stringr::str_glue("NES={prettyNum(gsea_result$NES, digits=2)}\np-value = {prettyNum(gsea_result$pval, digits=2, format='E')}\nq-value = {prettyNum(gsea_result$padj, digits=2)}")) +
  }
    #geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    #geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
  g <- g +
    ggplot2::geom_hline(yintercept=0, colour="black") +
    ggplot2::geom_line(color=linecolor) + ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits=c(0, max(toPlot$x)),
                       expand=c(0, 0)) +
    ggplot2::geom_segment(data=data.frame(x=pathway),
                          mapping=ggplot2::aes(x=x, y=-diff/2,
                                      xend=x, yend=diff/2),
                          size=ticksSize) +

    ggplot2::theme_bw()  +
    ggplot2::theme(panel.border=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 0, l = 5.5),
                   plot.title = ggplot2::element_text(size = titleSize,
                                                      hjust = titleHjust)) +

    ggplot2::labs( y="Enrichment Score", title = title)


  rankData <- tibble::tibble(x=seq_along(statsAdj), xend=x,
                             y=0, yend=statsAdj,
                             PathwayGene=x %in% pathway)

  size = 0.01
  e <- ggplot2::ggplot(rankData) +
    ggplot2::geom_segment(ggplot2::aes(x=x, xend=xend, y=0, yend=yend),
                 color = "grey80", size=size) +
    ggplot2::geom_segment(ggplot2::aes(x=x, xend=xend, y=0, yend=yend, color=PathwayGene),
                 data=rankData %>% dplyr::filter(PathwayGene), color = "black", size=size) +
    ggplot2::scale_x_continuous(limits=c(0, max(toPlot$x)),
                       expand=c(0, 0)) +
    ggplot2::scale_y_continuous(limits=c(-1, 1),
                       expand=c(0, 0)) +
    ggplot2::labs(x = "Rank", y = "Ranking Metric") +
    ggplot2::scale_color_manual(values=c("grey80", "black")) +
    ggplot2::scale_alpha_manual(values=c(0.5,1)) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          legend.position = "none",
          #axis.line=element_blank(),
          panel.grid = ggplot2::element_blank(),
          plot.margin =  ggplot2::margin(t = 0.1, r = 5.5, b = 5.5, l = 5.5),
          panel.spacing = rep(ggplot2::unit(0,"null"),4)
    )
  patchwork::wrap_plots(g,e,ncol=1)
}

enhancedGseaResultsTbl <- function(gsea_results, threshold, ethresh) {
  gsea_results %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sig=dplyr::case_when(
                    abs(NES) > ethresh & padj <= threshold ~ "s",
                    abs(NES) > ethresh ~ "e" ,
                    TRUE ~ "n"
                  ))
}

gseaVolcanoPlot <- function(gsea_results, threshold=0.2, ethresh=1, bell=F) {
  enhanced_results <- eenhancedGseaResultsTbl(gsea_results, threshold, ethresh) %>%
    dplyr::mutate(padj=dplyr::case_when(bell ~ padj, TRUE ~ -log10(padj)))
  ylab <- if(bell) "Adjusted p-value"  else "-log10(Adjusted p-value)"

  ggplot(enhanced_results, aes(x=NES, y = padj, color=sig)) +
    geom_point() +
    scale_color_manual(values = c("goldenrod", "grey", "red")) +
    labs(x = "Normalized Enrichment Score", y = ylab) +
    theme_bw() +
    theme(legend.position="none")
}

gseaDotPlot <- function(gsea_results, threshold=0.2, ethresh=1, limits = NULL) {
  enhanced_results <- enhancedGseaResultsTbl(gsea_results, threshold, ethresh)
  ranked_results <-
    enhanced_results %>%
    dplyr::filter(sig == "s") %>%
    dplyr::mutate(gseaResultRank = sign(NES)*-1*log10(padj)) %>%
    dplyr::mutate(pathway = forcats::fct_reorder(.$pathway, .$gseaResultRank))
  if(nrow(ranked_results) < 1)
    rlang::abort("No significant results.")
  if(is.null(limits))
    limits <- range(ranked_results$NES)

  cr <- c(min(limits[1], 0), max(limits[2], 0))
  colors <- c("white")
  if(cr[1] < 0) colors <- c("red", "orange", colors)
  if(cr[2] > 0) colors <- c(colors, "purple", "black")

  ggplot2::ggplot(ranked_results, ggplot2::aes(x=1, y=pathway, size=-log10(padj), color=NES)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add=-.5)) +
    ggplot2::scale_color_gradientn(colors=colors, limits = cr) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())
}



maxGenePlot <- function(se, meta=NULL, norm=NULL, group="Type",
                        facet=NULL, facet_scales="fixed", facet_space="fixed",
                        palette=ggthemes::ptol_pal()) {
  maxGeneTbl <- maxGene(se, norm)
  aesthetic <- ggplot2::aes(x=Sample, y = Value, label=Geneid)
  if(!is.null(group))
    aesthetic <- utils::modifyList(aesthetic, ggplot2::aes(fill=!!rlang::ensym(group)))
  p <- ggplot2::ggplot(maxGeneTbl, aesthetic) +
    ggplot2::geom_col() +
    ggplot2::discrete_scale("fill", "group", palette = palette) +
    ggplot2::geom_text(angle=90) +
    ggplot2::labs(x = "Sample", y = "Frequency", fill = "Type") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))

  p
}

genePCAPlot <- function(se, norm=NULL, sample_column = "Sample", group="Type", palette=ggthemes::ptol_pal(),
                        pcs=2, nrow=1, n_genes = 500, label_points = T, fix_coords = T, label_col = "Patient", max.overlaps = 50,  ...) {
  pca <- genePCA(se, sample_column=sample_column,...)
  aesthetic <- ggplot2::aes(x=Sample, y = Value, label=Geneid)
  if(is.numeric(pcs)) {
    pc_names <- stringr::str_glue("PC{seq(pcs)}")
    pcs <- purrr::cross2(pc_names, pc_names, .filter = `>=`)
  }
  ps <- pcs %>%
    purrr::map(function(pcs) {
      pc1 <- pcs[[1]]; pc2 <- pcs[[2]]
      aesthetic <- ggplot2::aes(x=!!rlang::ensym(pc1), y=!!rlang::ensym(pc2), label=!!rlang::ensym(label_col))
      if(!is.null(group))
        aesthetic <- utils::modifyList(aesthetic, ggplot2::aes(color=!!rlang::ensym(group)))
      p <- ggplot2::ggplot(pca$tbl, aesthetic) +
        ggplot2::geom_point() +
        ggplot2::discrete_scale("color", "group", palette = palette) +
        ggplot2::labs(x=stringr::str_glue("{pc1} ({pca$prp[pc1]} % Var.)"),
                      y=stringr::str_glue("{pc2} ({pca$prp[pc2]} % Var.)"),
                      color="Type")

      if(label_points) {
        p <- p + ggrepel::geom_text_repel(show.legend = F, max.overlaps = max.overlaps)
      }
      pc_vals <- pca$tbl %>% dplyr::select(!!!pc1, !!!pc2) %>% purrr::flatten_dbl()
      limits <- c(min(pc_vals), max(pc_vals))
      if(fix_coords) {
        p <- p + ggplot2::coord_fixed(xlim=limits, ylim=limits)
      }
      p + ggplot2::theme_bw()
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

