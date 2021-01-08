
run_DESeq2 <- function (counts, target, varInt = "group", batch = NULL, locfunc = "median",
          fitType = "parametric", pAdjustMethod = "BH", cooksCutoff = TRUE, lfcShrinkMethod = "apeglm",
          independentFiltering = TRUE, alpha = 0.05, verbose=F, ...)
{
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = target,
                                        design = formula(paste("~", ifelse(!is.null(batch),
                                                                   paste(batch, "+"), ""), varInt)))
  if(verbose)  {
    cat("Design of the statistical model:\n")
    cat(paste(as.character(design(dds)), collapse = " "), "\n")
  }
  dds <- DESeq2::estimateSizeFactors(dds, locfunc = eval(as.name(locfunc)))

  if(verbose) {
    cat("\nNormalization factors:\n")
    print(sizeFactors(dds))
  }

  dds <- DESeq2::estimateDispersions(dds, fitType = fitType)
  conditions <- levels(SummarizedExperiment::colData(dds)[, varInt])

  # Perform the test with different base conditions in order to allow
  # for shrinkage estimates for all pairs.
  results <- conditions %>% head(-1) %>% purrr::map(function(refCond) {
      SummarizedExperiment::colData(dds)[[varInt]] <-
        stats::relevel(SummarizedExperiment::colData(dds)[[varInt]], ref=refCond)
      dds <- DESeq2::nbinomWaldTest(dds, ...)

      comparisons <- tail(conditions, -purrr::detect_index(conditions, `==`, refCond))

      results <- comparisons %>% purrr::map(function(testCond) {
         res <- DESeq2::results(
          dds, contrast = c(varInt, testCond, refCond), pAdjustMethod = pAdjustMethod,
          cooksCutoff = cooksCutoff, independentFiltering = independentFiltering, alpha = alpha)
          if(!is.null(lfcShrinkMethod)) {
            res <- DESeq2::lfcShrink(dds, coef = stringr::str_glue("{varInt}_{testCond}_vs_{refCond}"), type = lfcShrinkMethod)
          }
         res
        }) %>% purrr::set_names(stringr::str_glue("{refCond}_vs_{comparisons}"))
    }) %>% purrr::flatten()
  tbl_results <- purrr::map(results, tibble::as_tibble, rownames = "Geneid")

  return(list(dds = dds, results = results, tbl_results = tbl_results, sf = DESeq2::sizeFactors(dds)))
}

#' Compute differential expression between two conditions using the wilcoxon rank sum test
#'
#' @param expr a matrix of genes (rows) by samples
#' @param condition identifies the condition of the columns. Only supports two conditions.
#'                  The first of the  two conditions in sort order will be the base condition.
#'                  An ordered factor can be used to force the baseline condition if its name
#'                  is not alphabetically first.
#' @param adjustMethod method for adjusting the p-values
#'
#' @return A tibble of differential expression.
wilcoxon_de <- function(expr, condition, adjustMethod="BH") {
  if(length(unique(condition)) != 2) stop("wilcoxon_de only works with 2 conditions.")

  base_cond <- as.character(sort(condition)[1])
  test_cond <- as.character(condition[condition != base_cond][1])
  base_ind <- which(condition == base_cond)
  test_ind <- which(condition == test_cond)

  c(M,N) %<-% dim(expr)
  N1 <- length(test_ind)
  N0 <- N-N1

  # minimum value if all the values in the test group are less than the base
  minW <- sum(1:N1)

  ranks <- apply(expr, 1, rank)
  ranksum <- apply(ranks[test_ind,], 2, sum)

  wilcox.p <- pwilcox(ranksum-minW, N0, N1)
  p <- ifelse(ranksum-minW < N1*N0/2,  2*wilcox.p, 2*(1-wilcox.p))
  padj <- p.adjust(p, method=adjustMethod)

  base_mean <- apply(expr[,base_ind], 1, mean)
  test_mean <- apply(expr[,test_ind], 1, mean)
  bottom <- min(c(base_mean, test_mean), na.rm = T)
  corr <- min(bottom, 0) - 1
  l2fc <- log2(test_mean-corr)-log2(base_mean-corr)

  tibble::tibble(Geneid=rownames(expr), base_mean=base_mean, test_mean=test_mean, log2FoldChange=l2fc, p=p, padj=padj)
}
