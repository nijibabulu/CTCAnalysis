
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
            res <- DESeq2::lfcShrink(dds, coef = str_glue("{varInt}_{testCond}_vs_{refCond}"), type = lfcShrinkMethod)
          }
         res
        }) %>% purrr::set_names(str_glue("{refCond}_vs_{comparisons}"))
    }) %>% purrr::flatten()
  tbl_results <- purrr::map(results, tibble::as_tibble, rownames = "Geneid")

  return(list(dds = dds, results = results, tbl_results = tbl_results, sf = DESeq2::sizeFactors(dds)))
}
