geneSummaries <- function(mat, meta=NULL) {
  mat <- ensureLongAttribute(mat, meta)
  target <- countsTarget(mat)
  countsLong <- longCounts(mat)
  countsLong %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarize(NumberOfGenes = sum(Value > 0), NumberOfNulls = sum(Value == 0)) %>%
    left_join(target %>% dplyr::mutate_if(is.factor, as.character), by=c(Sample="label"))
}

maxGene <- function(mat, meta=NULL, norm=NULL, frequency=T) {
  mat <- ensureLongAttribute(mat, meta)
  target <- countsTarget(mat)
  if(!is.null(norm)) {
    mat <- norm(mat, meta)
    mat <- ensureLongAttribute(mat)
  }

  longCounts(mat) %>%
    dplyr::group_by(Sample) %>%
    dplyr::arrange(Value) %>%
    {if(frequency) dplyr::mutate(.,Value=Value/sum(Value)) else {.}} %>%
    dplyr::slice(n())
}

genePCA <- function(mat, meta=NULL, norm=NULL, vst=DESeq2::varianceStabilizingTransformation,
                    n_genes = min(500, nrow(mat))) {
  if(!is.null(norm)) {
    mat <- norm(mat, meta)
    mat <- ensureLongAttribute(mat)
  }

  trans_mat <- vst(mat)
  rowVar <- apply(trans_mat, 1, var, na.rm = T)
  trans_matFilt <- trans_mat[order(rowVar, decreasing = T), ][1:n_genes, ]
  pca <- stats::prcomp(t(trans_matFilt))
  prp <- pca$sdev^2 * 100/sum(pca$sdev^2)
  prp <- round(prp, 2)
  names(prp) <- colnames(pca$x)
  pca_tbl <- pca$x %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Sample = rownames(pca$x)) %>%
    dplyr::left_join(countsTarget(lungCounts)  %>% dplyr::mutate_if(is.factor, as.character),
                     by=c(Sample="label"))
  list(pca=pca, prp=prp, tbl=pca_tbl)
}

sampleHclust <- function(mat, meta=NULL, vst=DESeq2::varianceStabilizingTransformation,
                         clust_method = "ward.D") {
  mat <- vst(mat)
  hc <- hclust(dist(t(mat)), method = clust_method)
  hc
}
