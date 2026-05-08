# Unwanted pregnancy analysis using CSA analysis pipeline
# ------------------------------------------------------

source("R_files/CSA_analysis.R")

read_unwanted_covariate_table <- function(path, ids) {
  df <- utils::read.csv(path, check.names = FALSE)
  cov_idx <- 2:6
  cov_idx <- cov_idx[cov_idx <= ncol(df)]
  cov_cols <- names(df)[cov_idx]
  if (length(ids) != nrow(df)) {
    stop(sprintf("ID length mismatch for %s: %d ids vs %d rows.", path, length(ids), nrow(df)), call. = FALSE)
  }
  out <- data.frame(id = as.numeric(ids), stringsAsFactors = FALSE)
  for (j in seq_along(cov_cols)) {
    out[[make.names(cov_cols[j])]] <- as.numeric(df[[cov_idx[j]]])
  }
  out <- out[is.finite(out$id), , drop = FALSE]
  out <- out[!duplicated(out$id), , drop = FALSE]
  out
}

run_unwanted_preg_analysis <- function(
  data_dir = "Unwanted_pregnancy_data",
  alpha = 0.05,
  gamma_search_grid = seq(1, 5, by = 0.1),
  output_csv = "PDF_images/unwanted_preg_results.csv",
  n_cores = detect_available_cores()
) {
  assert_required_packages()
  assert_real_data_new_packages()

  pair_map <- utils::read.csv(file.path(data_dir, "matched_pair_whole_data.csv"), check.names = FALSE)
  trt_cov <- read_unwanted_covariate_table(file.path(data_dir, "trt_var.csv"), pair_map[[1]])
  cont_cov <- read_unwanted_covariate_table(file.path(data_dir, "cont_var.csv"), pair_map[[2]])

  specs <- data.frame(
    outcome = c(
      "cesd", "depressed", "inter_personal", "low_positive", "somatic",
      "pwb", "autonomy", "environmental_mastery", "positive_relation", "purpose_in_life", "self_acceptance",
      "alcohol_dep", "drinking_days", "regular_at_risk",
      "SF_12", "smoking", "income_to_poverty_level"
    ),
    file = c(
      "cesd.csv", "depressed.csv", "inter_personal.csv", "low_positive.csv", "somatic.csv",
      "pwb.csv", "autonomy.csv", "environmental_mastery.csv", "positive_relation.csv", "purpose_in_life.csv", "self_acceptance.csv",
      "alcohol_dep.csv", "drinking_days.csv", "regular_at_risk.csv",
      "SF_12.csv", "smoking.csv", "poverty_level.csv"
    ),
    direction = c(
      rep("trt>cont", 5),
      rep("trt<cont", 6),
      rep("trt>cont", 3),
      "trt<cont",
      "trt>cont",
      "trt<cont"
    ),
    stringsAsFactors = FALSE
  )

  task_fun <- function(i) {
    analyze_one_outcome_new(
      outcome_name = specs$outcome[i],
      outcome_file = file.path(data_dir, specs$file[i]),
      direction = specs$direction[i],
      trt_cov = trt_cov,
      cont_cov = cont_cov,
      alpha = alpha,
      gamma_search_grid = gamma_search_grid,
      K = 5,
      seed = 9000 + i
    )
  }

  n_cores <- as.integer(n_cores)
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  n_cores <- min(n_cores, nrow(specs))

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    cat(sprintf("Running with %d cores.\n", n_cores))
    flush.console()
    rows <- parallel::mclapply(seq_len(nrow(specs)), task_fun, mc.cores = n_cores, mc.preschedule = TRUE)
  } else {
    rows <- vector("list", nrow(specs))
    for (i in seq_len(nrow(specs))) {
      cat(sprintf("Analyzing %d/%d: %s\n", i, nrow(specs), specs$outcome[i]))
      flush.console()
      rows[[i]] <- task_fun(i)
    }
  }

  tab <- do.call(rbind, rows)
  utils::write.csv(tab, output_csv, row.names = FALSE)
  tab
}

# Example usage:
# tab <- run_unwanted_preg_analysis(
#   data_dir = "Unwanted_pregnancy_data",
#   alpha = 0.05,
#   gamma_search_grid = seq(1, 5, by = 0.1),
#   output_csv = "PDF_images/unwanted_preg_results.csv",
#   n_cores = detect_available_cores()
# )
# print(tab)
