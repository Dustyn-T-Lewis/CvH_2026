suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

normalize_cvh_subject_id <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_replace("^CR0+([0-9]+)$", "CR\\1") |>
    stringr::str_replace("^PPS0+([0-9]+)$", "PPS\\1")
}

normalize_cvh_col_id <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_replace("^CR0+([0-9]+)(_.+)$", "CR\\1\\2") |>
    stringr::str_replace("^PPS0+([0-9]+)(_.+)?$", "PPS\\1\\2")
}

derive_cvh_analysis_metadata <- function(pheno_meta) {
  required_cols <- c("sample_id", "pid", "timepoint", "supp", "cancer")
  missing_cols <- setdiff(required_cols, names(pheno_meta))
  if (length(missing_cols) > 0) {
    stop("CRm_meta missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  derived <- pheno_meta |>
    tibble::as_tibble() |>
    dplyr::transmute(
      Col_ID = normalize_cvh_col_id(sample_id),
      Subject_ID = normalize_cvh_subject_id(pid),
      Timepoint = stringr::str_trim(timepoint),
      Supplement = dplyr::na_if(stringr::str_trim(supp), ""),
      Cancer = stringr::str_trim(cancer)
    ) |>
    dplyr::mutate(
      Group = dplyr::case_when(
        Cancer == "SURV" & Supplement == "CRE" ~ "CR_CRE",
        Cancer == "SURV" & Supplement == "PLA" ~ "CR_PLA",
        Cancer == "CTL"                        ~ "PPS",
        TRUE                                   ~ NA_character_
      ),
      Group_Time = dplyr::case_when(
        Group == "CR_CRE" ~ paste0("CRE_", Timepoint),
        Group == "CR_PLA" ~ paste0("PLA_", Timepoint),
        Group == "PPS" & Timepoint == "T1" ~ "H_T1",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::select(Col_ID, Subject_ID, Group, Timepoint, Group_Time, Supplement) |>
    dplyr::distinct()

  assert_cvh_design_rules(derived, context = "derived CRm_meta")
  derived
}

assert_cvh_design_rules <- function(meta,
                                    context = "CvH metadata",
                                    allow_t2_only_singletons = FALSE) {
  required_cols <- c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time", "Supplement")
  missing_cols <- setdiff(required_cols, names(meta))
  if (length(missing_cols) > 0) {
    stop(context, " missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  meta_tbl <- meta |>
    tibble::as_tibble() |>
    dplyr::mutate(
      Col_ID = normalize_cvh_col_id(Col_ID),
      Subject_ID = normalize_cvh_subject_id(Subject_ID),
      Timepoint = stringr::str_trim(Timepoint),
      Group = stringr::str_trim(Group),
      Group_Time = stringr::str_trim(Group_Time),
      Supplement = dplyr::na_if(stringr::str_trim(as.character(Supplement)), "")
    )

  if (anyDuplicated(meta_tbl$Col_ID)) {
    dup_ids <- unique(meta_tbl$Col_ID[duplicated(meta_tbl$Col_ID)])
    stop(context, " has duplicated Col_ID values: ", paste(dup_ids, collapse = ", "))
  }
  if (anyDuplicated(meta_tbl$Subject_ID)) {
    dup_subjects <- meta_tbl |>
      dplyr::count(Subject_ID) |>
      dplyr::filter(n > 2) |>
      dplyr::pull(Subject_ID)
    if (length(dup_subjects) > 0) {
      stop(context, " has subjects with >2 samples: ", paste(dup_subjects, collapse = ", "))
    }
  }

  allowed_timepoints <- c("T1", "T2")
  bad_timepoints <- setdiff(unique(meta_tbl$Timepoint), allowed_timepoints)
  if (length(bad_timepoints) > 0) {
    stop(context, " has unsupported Timepoint values: ", paste(bad_timepoints, collapse = ", "))
  }

  allowed_groups <- c("CR_CRE", "CR_PLA", "PPS")
  bad_groups <- setdiff(unique(meta_tbl$Group), allowed_groups)
  if (length(bad_groups) > 0) {
    stop(context, " has unsupported Group values: ", paste(bad_groups, collapse = ", "))
  }

  allowed_group_time <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
  bad_group_time <- setdiff(unique(meta_tbl$Group_Time), allowed_group_time)
  if (length(bad_group_time) > 0) {
    stop(context, " has unsupported Group_Time values: ", paste(bad_group_time, collapse = ", "))
  }

  expected_group <- dplyr::case_when(
    meta_tbl$Group_Time %in% c("CRE_T1", "CRE_T2") ~ "CR_CRE",
    meta_tbl$Group_Time %in% c("PLA_T1", "PLA_T2") ~ "CR_PLA",
    meta_tbl$Group_Time == "H_T1"                  ~ "PPS",
    TRUE                                           ~ NA_character_
  )
  if (!all(meta_tbl$Group == expected_group)) {
    bad_rows <- meta_tbl$Col_ID[meta_tbl$Group != expected_group]
    stop(context, " has inconsistent Group vs Group_Time for: ", paste(bad_rows, collapse = ", "))
  }

  expected_supp <- dplyr::case_when(
    meta_tbl$Group == "CR_CRE" ~ "CRE",
    meta_tbl$Group == "CR_PLA" ~ "PLA",
    meta_tbl$Group == "PPS"    ~ NA_character_,
    TRUE                       ~ NA_character_
  )
  same_supp <- (is.na(meta_tbl$Supplement) & is.na(expected_supp)) |
    dplyr::coalesce(meta_tbl$Supplement == expected_supp, FALSE)
  if (!all(same_supp)) {
    bad_rows <- meta_tbl$Col_ID[!same_supp]
    stop(context, " has inconsistent Supplement values for: ", paste(bad_rows, collapse = ", "))
  }

  expected_timepoint <- dplyr::case_when(
    meta_tbl$Group_Time %in% c("CRE_T1", "PLA_T1", "H_T1") ~ "T1",
    meta_tbl$Group_Time %in% c("CRE_T2", "PLA_T2")         ~ "T2",
    TRUE                                                   ~ NA_character_
  )
  if (!all(meta_tbl$Timepoint == expected_timepoint)) {
    bad_rows <- meta_tbl$Col_ID[meta_tbl$Timepoint != expected_timepoint]
    stop(context, " has inconsistent Timepoint vs Group_Time for: ", paste(bad_rows, collapse = ", "))
  }

  healthy_bad <- meta_tbl |>
    dplyr::filter(Group == "PPS", Timepoint != "T1")
  if (nrow(healthy_bad) > 0) {
    stop(context, " includes healthy samples outside T1: ", paste(healthy_bad$Col_ID, collapse = ", "))
  }

  cr_by_subject <- meta_tbl |>
    dplyr::filter(Group != "PPS") |>
    dplyr::group_by(Subject_ID) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_supp = dplyr::n_distinct(Supplement),
      has_t1 = any(Timepoint == "T1"),
      has_t2 = any(Timepoint == "T2"),
      .groups = "drop"
    )

  inconsistent_cr <- cr_by_subject |>
    dplyr::filter(
      n_supp != 1 |
        n_samples > 2 |
        (!allow_t2_only_singletons & !has_t1 & has_t2)
    )
  if (nrow(inconsistent_cr) > 0) {
    stop(
      context, " has inconsistent CR repeated-measures structure for subjects: ",
      paste(inconsistent_cr$Subject_ID, collapse = ", ")
    )
  }

  invisible(meta_tbl)
}

validate_cvh_sample_ids <- function(sample_ids,
                                    meta,
                                    context = "CvH sample match",
                                    allow_t2_only_singletons = FALSE) {
  if (anyDuplicated(sample_ids)) {
    dup_ids <- unique(sample_ids[duplicated(sample_ids)])
    stop(context, " has duplicated matrix sample IDs: ", paste(dup_ids, collapse = ", "))
  }

  meta_tbl <- assert_cvh_design_rules(
    meta,
    context = context,
    allow_t2_only_singletons = allow_t2_only_singletons
  )

  missing_in_meta <- setdiff(sample_ids, meta_tbl$Col_ID)
  missing_in_mat <- setdiff(meta_tbl$Col_ID, sample_ids)
  if (length(missing_in_meta) > 0 || length(missing_in_mat) > 0) {
    msg <- c()
    if (length(missing_in_meta) > 0) {
      msg <- c(msg, paste0("matrix-only: ", paste(missing_in_meta, collapse = ", ")))
    }
    if (length(missing_in_mat) > 0) {
      msg <- c(msg, paste0("metadata-only: ", paste(missing_in_mat, collapse = ", ")))
    }
    stop(context, " mismatch between matrix and metadata sample IDs (", paste(msg, collapse = " | "), ")")
  }

  meta_tbl |>
    dplyr::mutate(.order = match(Col_ID, sample_ids)) |>
    dplyr::arrange(.order) |>
    dplyr::select(-.order)
}

load_cvh_analysis_metadata <- function(meta_file,
                                       pheno_file = NULL,
                                       raw_sample_ids = NULL,
                                       allow_t2_only_singletons = FALSE) {
  meta_tbl <- readr::read_csv(meta_file, show_col_types = FALSE) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      Col_ID = normalize_cvh_col_id(Col_ID),
      Subject_ID = normalize_cvh_subject_id(Subject_ID)
    )

  meta_tbl <- assert_cvh_design_rules(
    meta_tbl,
    context = basename(meta_file),
    allow_t2_only_singletons = allow_t2_only_singletons
  )

  if (!is.null(pheno_file) && file.exists(pheno_file)) {
    derived_tbl <- derive_cvh_analysis_metadata(readr::read_csv(pheno_file, show_col_types = FALSE))
    compare_cols <- c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time", "Supplement")
    analysis_cmp <- meta_tbl |>
      dplyr::select(dplyr::all_of(compare_cols)) |>
      dplyr::arrange(Col_ID)
    derived_cmp <- derived_tbl |>
      dplyr::select(dplyr::all_of(compare_cols)) |>
      dplyr::arrange(Col_ID)

    if (!identical(analysis_cmp, derived_cmp)) {
      joined <- dplyr::full_join(
        analysis_cmp |> dplyr::rename_with(~ paste0(.x, "_analysis"), -Col_ID),
        derived_cmp |> dplyr::rename_with(~ paste0(.x, "_derived"), -Col_ID),
        by = "Col_ID"
      )
      stop(
        basename(meta_file), " does not match metadata derived from ",
        basename(pheno_file), ". First mismatched rows:\n",
        paste(capture.output(print(head(joined, 10))), collapse = "\n")
      )
    }
  }

  if (!is.null(raw_sample_ids)) {
      meta_tbl <- validate_cvh_sample_ids(
        raw_sample_ids,
        meta_tbl,
        context = basename(meta_file),
        allow_t2_only_singletons = allow_t2_only_singletons
      )
  }

  meta_tbl
}
