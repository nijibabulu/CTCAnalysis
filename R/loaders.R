# TODO: move parsing summaries to a separate function
loadFeatureCountData <- function(countFile, pdata=NULL, pdataNameColumn="Sample", parseSummary=F, keep_positions=F, summaryFile=paste0(countFile, ".summary")) {
  counts <- readr::read_tsv(countFile,
                            comment =  "#",
                            col_types = readr::cols(
                              .default = readr::col_double(),
                              Geneid = readr::col_character(),
                              Chr = readr::col_character(),
                              Start = readr::col_character(),
                              End = readr::col_character(),
                              Strand = readr::col_character()
                            ))

  geneInfo <- counts[,1:6]
  chroms <- stringr::str_split(geneInfo$Chr, ";") %>% purrr::map_chr(~names(sort(table(.x), decreasing = T))[1])
  starts <- stringr::str_split(geneInfo$Start, ";") %>% purrr::map_int(purrr::compose(min, as.integer))
  ends <- stringr::str_split(geneInfo$End, ";") %>% purrr::map_int(purrr::compose(max, as.integer))
  strands <- stringr::str_split(geneInfo$Strand, ";") %>% purrr::map_chr(~names(sort(table(.x), decreasing = T))[1])
  geneRanges <- GenomicRanges::GRanges(seqnames = chroms,
                                       ranges = IRanges::IRanges(start= starts, end = ends),
                                       strand = strands,
                                       Geneid = geneInfo$Geneid,
                                       Length = geneInfo$Length)

  count_mat <- counts[,-c(2:6)] %>%
    dplyr::mutate(Geneid=make.names(Geneid, unique=T)) %>%
    tibble::column_to_rownames("Geneid") %>%
    as.matrix()
  if(!is.null(pdata)) {
    pdata <- pdata %>% dplyr::rename(name=pdataNameColumn)
  }
  if(parseSummary) {
    summary <- readr::read_tsv(summaryFile,
                               col_types = readr::cols( .default = readr::col_double(), Status = readr::col_character())) %>%
      tidyr::pivot_longer(-Status) %>% tidyr::pivot_wider(name, Status)
    if(is.null(pdata)) {
      pdata <- summary
    } else {
      common_colnames <- intersect(colnames(summary), colnames(pdata))
      if(length(common_colnames)) {
        stop(stringr::str_glue("pdata has columns in common with summary {stringr::str_c(common_colnames, collapse=', ')}"))
      }
      pdata <- dplyr::full_join(summary, pdata, by="name")
    }
  }

  if(is.null(pdata)){
    pdata <- S4Vectors::DataFrame(x = seq(ncol(count_mat)), row.names = colnames(count_mat))
  } else {
    pdata <- pdata %>% tibble::column_to_rownames("name")
  }

  SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=count_mat),
                                             #rowData=geneInfo,
                                             rowRanges=geneRanges,
                                             colData=pdata)
}

loadLibrarySizes <- function(se, libsz_file, compute_unmapped=F,
                             unmapped_col = "Unassigned_Unmapped",
                             mapped_cols = c("Assigned", "Unassigned_Read_Type", "Unassigned_Singleton", "Unassigned_MappingQuality",
                                             "Unassigned_Chimera", "Unassigned_FragmentLength", "Unassigned_Duplicate", "Unassigned_MultiMapping",
                                             "Unassigned_Secondary", "Unassigned_NonSplit",          "Unassigned_NoFeatures",        "Unassigned_Overlapping_Length",
                                             "Unassigned_Ambiguity" )) {
  libsz_tbl <- readr::read_tsv(libsz_file, col_names = c("name", "LibrarySize"), col_types = "ci")
  coldata_tbl <- SummarizedExperiment::colData(se) %>% as.data.frame() %>% tibble::rownames_to_column("name") %>% dplyr::left_join(libsz_tbl, by="name")
  if(compute_unmapped) {
    coldata_tbl <- dplyr::rowwise(coldata_tbl) %>%
      dplyr::mutate(!!unmapped_col :=  LibrarySize - sum(!!as.name(mapped_cols)))
  }
  SummarizedExperiment::colData(se) <-  coldata_tbl %>% tibble::column_to_rownames("name") %>% as("DataFrame")
  se
}


loadTarget <-
  function (targetFile, varInt, condRef, batch=NULL, verbose=F, group_translation=NULL, keep_group_code=T)
{
  target <- read.table(targetFile, header = TRUE, sep = "\t",
                       na.strings = "")
  if (!I(varInt %in% names(target)))
    stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target)))
    stop(paste("The batch effect", batch, "is not in the target file"))
  target[, varInt] <- as.factor(target[, varInt])
  if (!I(condRef %in% as.character(target[, varInt])))
    stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
  lev <- c(condRef, unique(as.character(target[, varInt])))
  lev <- lev[!duplicated(lev)]
  target[, varInt] <- factor(target[, varInt], levels = lev)
  target <- target[order(target[, varInt]), ]
  rownames(target) <- as.character(target[, 1])
  if (min(table(target[, varInt])) < 2)
    stop(paste("The factor of interest", varInt, "has a level without replicates"))
  if (any(is.na(cbind(target[, c(varInt, batch)], target[,
                                                         1:2]))))
    stop("NA are present in the target file")
  if (!is.null(batch) && is.numeric(target[, batch]))
    warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  if (any(grepl("[[:punct:]]", as.character(target[, varInt]))))
    stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
  if(!is.null(group_translation)) {
    if(keep_group_code) {
      target <- dplyr::mutate(target, group_code = target[, varInt])
    }
    target <- dplyr::mutate(target, !!varInt := translated_factor(group, group_translation))
  }
  if(verbose) {
    cat("Target file:\n")
    print(target)
  }
  return(target)
  }

joinCountFile <-
  function (tbl, count_file, label, by, delim="\t")
{
    counts <- readr::read_delim(count_file, delim = delim, col_names = c(by, label), col_types = "ci")
    dplyr::left_join(tbl, counts, by = by)
}

loadSummaries <-
  function (target, rawDir = "raw", files_col = "file", labels_col = "label",
            library_sizes_file = NULL, multimapped_reads_file = NULL,
            mapped_reads_file = NULL, duplicate_reads_file = NULL, skip = 0)
{
  labels <- as.character(target[, labels_col])
  files <- as.character(target[, files_col])
  if(fs::file_exists(fs::path(rawDir, stringr::str_glue("{files[1]}.summary")))) {
    table <-
      purrr::map2(labels, files,
                 ~readr::read_delim(file = fs::path(rawDir, stringr::str_glue("{..2}.summary")),
                             delim="\t",
                             col_names = c("stat", "val"),
                             skip = 1, col_types="ci") %>%
                   dplyr::mutate(Library = as.character(..1))) %>%
      dplyr::bind_rows() %>%
      dplyr::filter(val > 0) %>%
      tidyr::pivot_wider(names_from = "stat", values_from = "val")

    if(!purrr::is_null(library_sizes_file)) {
      table <- joinCountFile(table, library_sizes_file, "LibrarySize", "Library")
    }
    if(!purrr::is_null(multimapped_reads_file)) {
      table <- joinCountFile(table, multimapped_reads_file, "Multimapped", "Library")
    }
    if(!purrr::is_null(mapped_reads_file)) {
      table <- joinCountFile(table, mapped_reads_file, "Mapped", "Library")
    }
    if(!purrr::is_null(duplicate_reads_file)) {
      table <- readr::read_delim(duplicate_reads_file, delim="\t", col_types=readr::cols()) %>% dplyr::full_join(table, ., by=c("Library"="name"))
    }
    table
  } else {
    stop("Only featureCounts count summaries are currently supported.")
  }
}

ensureTargetAttribute <- function(counts, target) {
  attr <- attributes(counts)
  attr$target <- target
  attributes(counts) <- attr
  counts
}

loadCounts <-
  function (target, rawDir = "raw", skip = 0, files_col = "file", labels_col = "label",
            featuresToRemove = c("alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual"),
            verbose = F)
{
  labels <- as.character(target[, labels_col])
  files <- as.character(target[, files_col])
  f1 <- read.table(file.path(rawDir, files[1]), sep = "\t",
                   quote = "\"", header = FALSE, skip = 5, nrows = 5, stringsAsFactors = FALSE)
  if (ncol(f1) >= 7 && is.numeric(f1[, 7])) {
    idCol <- 1
    countsCol <- 7
    header <- TRUE
  }
  else {
    if (ncol(f1) >= 2 && is.numeric(f1[, 2])) {
      idCol <- 1
      countsCol <- 2
      header <- FALSE
    }
    else {
      stop("Can't determine if count files come from HTSeq-count or featureCounts")
    }
  }
  rawCounts <- read.table(file.path(rawDir, files[1]), sep = "\t",
                          quote = "\"", header = header, skip = skip, stringsAsFactors = FALSE)
  rawCounts <- rawCounts[, c(idCol, countsCol)]
  colnames(rawCounts) <- c("Id", labels[1])
  if (any(duplicated(rawCounts$Id))) {
    stop("Duplicated feature names in ", files[1], ": ",
         paste(unique(rawCounts$Id[duplicated(rawCounts$Id)]),
               collapse = ", "))
  }
  if(verbose) {
    cat("Loading files:\n")
    cat(files[1], ": ", length(rawCounts[, labels[1]]), " rows and ",
        sum(rawCounts[, labels[1]] == 0), " null count(s)\n",
        sep = "")
  }
  for (i in 2:length(files)) {
    tmp <- read.table(file.path(rawDir, files[i]), sep = "\t",
                      quote = "\"", header = header, skip = skip, stringsAsFactors = FALSE)
    tmp <- tmp[, c(idCol, countsCol)]
    colnames(tmp) <- c("Id", labels[i])
    if (any(duplicated(tmp$Id))) {
      stop("Duplicated feature names in ", files[i], ": ",
           paste(unique(tmp$Id[duplicated(tmp$Id)]), collapse = ", "))
    }
    rawCounts <- merge(rawCounts, tmp, by = "Id", all = TRUE)
    if(verbose) {
      cat(files[i], ": ", length(tmp[, labels[i]]), " rows and ",
          sum(tmp[, labels[i]] == 0), " null count(s)\n",
          sep = "")
    }
  }
  rawCounts[is.na(rawCounts)] <- 0
  counts <- as.matrix(rawCounts[, -1])
  rownames(counts) <- rawCounts[, 1]
  counts <- counts[order(rownames(counts)), ]
  if (any(counts%%1 != 0))
    stop("Input counts are not integer values as required by DESeq2 and edgeR.")
  if(verbose)
    cat("\nFeatures removed:\n")
  for (f in setdiff(featuresToRemove, "")) {
    match <- grep(f, rownames(counts))
    if (length(match) > 0) {
      if(verbose)
        cat(rownames(counts)[match], sep = "\n")
      counts <- counts[-match, ]
    }
  }
  if(verbose) {
    cat("\nTop of the counts matrix:\n")
    print(head(counts))
    cat("\nBottom of the counts matrix:\n")
    print(tail(counts))
  }

  counts <- ensureTargetAttribute(counts, target)


  return(counts)
}

loadMeta <-
  function (target, rawDir = "raw", files_col = "file")
{
  files <- as.character(target[, files_col])
  f1 <- readr::read_tsv(file.path(rawDir, files[1]),
                        quote = "\"", comment = "#")
  if(ncol(f1) >= 7) {
    count_col <- colnames(f1) %>% tail(n=1)
    return(f1 %>% dplyr::arrange(Geneid) %>% dplyr::select(-!!count_col))
  } else {
    stop("Can only load metadata from a featureCounts file")
  }
}


matrixToLong <-
  function(mat, meta=NULL, id="Geneid", names_to="Sample", values_to="Value") {
  long <- mat %>%
      tibble::as_tibble(rownames=id) %>%
      tidyr::pivot_longer(-!!id, names_to=names_to, values_to=values_to)

  target <- countsTarget(mat)
  long <- dplyr::left_join(long, target %>% dplyr::mutate_if(is.factor, as.character), by=c(Sample="label"))

  # XXX should we check if there is an update?
  if(!is.null(meta)) {
    long %>% dplyr::left_join(meta, by=id)
  } else {
    long
  }
}


ensureLongAttribute <- function(mat, meta=NULL, force_update=F, id="Geneid") {
    attr <- attributes(mat)
    if(!"long" %in% names(attr) || force_update) {
      attr$long <- matrixToLong(mat, meta)
    }
    if(!is.null(meta)) {
      if(!"long" %in% names(attr) || force_update) {
        attr$meta <- meta
      }
    }
    attributes(mat) <- attr
    mat
}

longCounts <- function(mat) {
  attr <- attributes(mat)
  if(!"long" %in% names(attr)) {
    stop("Cannot get long counts from the attributes of the counts matrix")
  }
  attr$long
}

countsTarget <- function(mat) {
  attr <- attributes(mat)
  if(!"target" %in% names(attr)) {
    stop("Cannot get target from the counts matrix. Were these counts loaded with loadCounts()?")
  }
  attr$target
}
