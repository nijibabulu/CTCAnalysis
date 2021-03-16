geneSummaries.tbl <- function(mat, meta=NULL) {
  mat <- ensureLongAttribute(mat, meta)
  target <- countsTarget(mat)
  countsLong <- longCounts(mat)
  countsLong %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarize(NumberOfGenes = sum(Value > 0), NumberOfNulls = sum(Value == 0)) %>%
    dplyr::left_join(target %>% dplyr::mutate_if(is.factor, as.character), by=c(Sample="label"))
}

geneSummaries <- function(se) {
  .summary_tbl <- function(f, v) SummarizedExperiment::assay(se) %>% apply(2,f) %>% tibble::enframe(name = "Sample", value=v)
  dplyr::inner_join(.summary_tbl(function(x) sum(x > 0), "NumberOfGenes"),
                    .summary_tbl(function(x) sum(x == 0), "NumberOfNulls"),
                    by="Sample") %>%
    dplyr::left_join(colTbl(lung_data), by="Sample")
}

maxGene.tbl <- function(mat, meta=NULL, norm=NULL, frequency=T) {
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
    dplyr::slice(dplyr::n())
}

maxGene<- function(se, norm=NULL, frequency=T) {
  mat <- SummarizedExperiment::assay(se)
  if(!is.null(norm)) {
    mat <- norm(mat, colTbl(se))
  }

  tibble::enframe(apply(mat, 2, which.max), name="Sample", value="Index") %>%
    dplyr::inner_join(tibble::enframe(apply(mat, 2, max), name="Sample", value="Value"), by="Sample") %>%
    dplyr::inner_join(tibble::enframe(colSums(mat), name="Sample", value="LibSize"), by="Sample") %>%
    dplyr::mutate(Geneid=rownames(mat)[.$Index]) %>%
    {if(frequency) dplyr::mutate(.,Value=Value/LibSize) else {.}} %>%
    dplyr::left_join(colTbl(se), by="Sample")
}

genePCA.tbl <- function(mat, meta=NULL, norm=NULL, sample_column = "label",
                    vst=DESeq2::varianceStabilizingTransformation,
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
    dplyr::left_join(meta  %>% dplyr::mutate_if(is.factor, as.character),
                     by=c(Sample = sample_column))
  list(pca=pca, prp=prp, tbl=pca_tbl)
}

genePCA <- function(se, norm=NULL, sample_column = "Sample",
                    vst=DESeq2::varianceStabilizingTransformation,
                    n_genes = min(500, nrow(mat))) {
  mat <- SummarizedExperiment::assay(se)
  if(!is.null(norm))
    mat <- norm(mat, colTbl(se))

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
    dplyr::left_join(colTbl(se),
                     by=c(Sample = sample_column))
  list(pca=pca, prp=prp, tbl=pca_tbl)
}


sampleHclust <- function(mat, meta=NULL, vst=DESeq2::varianceStabilizingTransformation,
                         clust_method = "ward.D") {
  mat <- vst(mat)
  hc <- hclust(dist(t(mat)), method = clust_method)
  hc
}
