translated_factor <- function(vec, dict) factor(vec, levels=names(dict), labels=dict)

# add the sample variables (patient, celltype, stage) into the data
expand_sample_names <- function(tbl, col_name="Sample", stage_names=NULL, treatment_names=NULL, patient_type_names=NULL) {
  tbl %>%
    tidyr::separate(!!col_name, c("Study", "Type", "Stage", "Patient"), "_", remove=F) %>%
    tidyr::separate(Patient, c("PatientType", "PatientId"), "-", remove=F) %>%
    {if(!purrr::is_null(stage_names))
      dplyr::mutate(., Stage=translated_factor(.$Stage, stage_names)) else .} %>%
    {if(!purrr::is_null(treatment_names))
      dplyr::mutate(., Treatment=translated_factor(.$Treatment, treatment_names)) else .} %>%
    {if(!purrr::is_null(patient_type_names))
      dplyr::mutate(., PatientType=translated_factor(.$PatientType, patient_type_names)) else .}
}

# add the column data to the summarized experiment structure
expand_sample_names_to_col_data <- function(se, ...)   {
  SummarizedExperiment::colData(se) <-
    SummarizedExperiment::colData(se) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    expand_sample_names(...) %>%
    tibble::column_to_rownames("Sample") %>%
    as("DataFrame")
  se
}



drop_samples <- function(counts, target, sample_names, sample_var="label") {
  is_factor <- target %>% select(!!sample_var) %>% summarize_all(class) == "factor"
  target <- target %>%
    mutate(!!sample_var := as.character(!!sym(sample_var)))

  target <- target[!target[,sample_var] %in% sample_names,]

  if(is_factor) {
    target[,sample_var] <- factor(target[,sample_var])
  }

  counts <- counts[, !colnames(counts) %in% sample_names]
  counts <- ensureTargetAttribute(counts, target)

  list(counts = counts, target = target)
}
