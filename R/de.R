run_DESeq2_results <- function(testCond, refCond, dds, varInt = "Type", pAdjustMethod = "BH", cooksCutoff = TRUE,
                               lfcShrinkMethod = "apeglm", independentFiltering = TRUE, alpha = 0.05) {
  res <- DESeq2::results(
    dds, contrast = c(varInt, testCond, refCond), pAdjustMethod = pAdjustMethod,
    cooksCutoff = cooksCutoff, independentFiltering = independentFiltering, alpha = alpha)
  if(is.null(lfcShrinkMethod)) {
    res
  } else {
    res_shrink <- DESeq2::lfcShrink(dds, coef = stringr::str_glue("{varInt}_{testCond}_vs_{refCond}"), type = lfcShrinkMethod)
    if(!all(rownames(res_shrink) == rownames(res))) {
      stop("shrunk DESeqRresults are not compatible. Cannot add stats in")
    }
    res_shrink$stat <- res$stat
    res_shrink
  }
}
run_DESeq2 <- function (se, varInt = "Type", batch = NULL, locfunc = "median",
          fitType = "parametric", pAdjustMethod = "BH", cooksCutoff = TRUE,
          lfcShrinkMethod = "apeglm", lrt = T,
          independentFiltering = TRUE, alpha = 0.05, verbose=F, ...)
{
  # suppress the DESeq2 warning about converting to a factor
  SummarizedExperiment::colData(se)[,varInt] <- factor(SummarizedExperiment::colData(se)[,varInt])

  # costruct a formula
  design_formula <- stringr::str_c("~ ", stringr::str_c(batch, varInt, sep = " + "))
  reduced_formula <- stringr::str_c("~ ", ifelse(is.null(batch), 1, batch))
  dds <- DESeq2::DESeqDataSet(se, design = formula(design_formula))
  if(verbose)  {
    cat("Design of the statistical model:\n")
    cat(paste(as.character(DESeq2::design(dds)), collapse = " "), "\n")
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
      if(lrt)
        dds_lrt <- DESeq2::nbinomLRT(dds, reduced = formula(reduced_formula), ...)

      comparisons <- tail(conditions, -purrr::detect_index(conditions, `==`, refCond))

      run_DESeq2_results_closure <- purrr::partial(run_DESeq2_results, refCond = refCond, varInt = varInt,
                                                   pAdjustMethod = pAdjustMethod, cooksCutoff = cooksCutoff,
                                                   lfcShrinkMethod = lfcShrinkMethod,  independentFiltering = independentFiltering,
                                                   alpha = alpha)
      res <- comparisons %>% purrr::map(run_DESeq2_results_closure, dds = dds) %>% purrr::set_names(stringr::str_glue("{refCond}_vs_{comparisons}"))
      if(lrt) {
        res_lrt <-
          comparisons %>% purrr::map(run_DESeq2_results_closure, dds = dds_lrt) %>% purrr::set_names(stringr::str_glue("{refCond}_vs_{comparisons}_lrt"))
        res <- c( res, res_lrt )
      }
      res
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
