# CSA analysis (both-sided tests)
# -------------------------------

source("R_files/CSA_analysis.R")

analyze_one_outcome_both <- function(
  outcome_name,
  outcome_file,
  trt_cov,
  cont_cov,
  alpha = 0.05,
  gamma_search_grid = seq(1, 5, by = 0.1),
  K = 5,
  seed = 123
) {
  out_df <- read_outcome_file_new(outcome_file)
  data_obj <- build_analysis_data_new(out_df, trt_cov, cont_cov)
  complete_pairs <- data_obj$complete_pairs
  aipw_data <- data_obj$aipw_data
  counts <- data_obj$counts
  x_cols <- names(trt_cov)[-1]

  if (nrow(complete_pairs) == 0) {
    stop(sprintf("Outcome '%s' has no complete pairs.", outcome_name), call. = FALSE)
  }
  if (nrow(aipw_data) < K || length(unique(aipw_data$T)) < 2) {
    stop(sprintf("Outcome '%s' has insufficient AIPW sample.", outcome_name), call. = FALSE)
  }

  match_side <- "two.sided"
  inc_side <- "two.sided"
  frac_complete <- counts[["n"]] / counts[["total_pairs"]]

  match_res <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = 1,
    side = match_side,
    score = "huber",
    huber_c = 1.345,
    mc_reps = 0,
    seed = seed
  )
  prop_res <- combined_weighted_test_new(
    complete_pairs = complete_pairs,
    aipw_data = aipw_data,
    counts = counts,
    alpha = alpha,
    Gamma_M = 1,
    Gamma_A = 1,
    match_side = match_side,
    incomplete_side = inc_side,
    K = K,
    x_cols = x_cols,
    seed = seed
  )

  match_gamma <- if (match_res$p_value <= alpha) {
    find_min_gamma_accept_new(function(g) {
      matched_sensitivity_pvalue(
        complete_pairs = complete_pairs,
        Gamma = g,
        side = match_side,
        score = "huber",
        huber_c = 1.345,
        mc_reps = 0
      )$p_value
    }, alpha = alpha, gamma_grid = gamma_search_grid)
  } else NA_real_

  prop_gamma <- if (prop_res$p_value <= alpha) {
    find_min_gamma_accept_new(function(g) {
      combined_weighted_test_new(
        complete_pairs = complete_pairs,
        aipw_data = aipw_data,
        counts = counts,
        alpha = alpha,
        Gamma_M = g,
        Gamma_A = g,
        match_side = match_side,
        incomplete_side = inc_side,
        K = K,
        x_cols = x_cols
      )$p_value
    }, alpha = alpha, gamma_grid = gamma_search_grid)
  } else NA_real_

  data.frame(
    Outcome = outcome_name,
    Complete_pair_fraction = formatC(frac_complete, format = "f", digits = 3),
    Matching = fmt_pval_decision_new(match_res$p_value, alpha = alpha),
    Matching_min_Gamma_accept = ifelse(is.na(match_gamma), "NA", formatC(match_gamma, format = "f", digits = 2)),
    Proposed = fmt_pval_decision_new(prop_res$p_value, alpha = alpha),
    Proposed_min_Gamma_accept = ifelse(is.na(prop_gamma), "NA", formatC(prop_gamma, format = "f", digits = 2)),
    stringsAsFactors = FALSE
  )
}

run_CSA_both_analysis <- function(
  data_dir = "CSA_data",
  alpha = 0.05,
  gamma_search_grid = seq(1, 5, by = 0.1),
  output_csv = "PDF_images/CSA_both_results.csv",
  n_cores = detect_available_cores()
) {
  assert_required_packages()
  assert_real_data_new_packages()

  trt_cov <- read_covariate_table_new(file.path(data_dir, "trt_var.csv"))
  cont_cov <- read_covariate_table_new(file.path(data_dir, "cont_var.csv"))

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
    stringsAsFactors = FALSE
  )

  task_fun <- function(i) {
    analyze_one_outcome_both(
      outcome_name = specs$outcome[i],
      outcome_file = file.path(data_dir, specs$file[i]),
      trt_cov = trt_cov,
      cont_cov = cont_cov,
      alpha = alpha,
      gamma_search_grid = gamma_search_grid,
      K = 5,
      seed = 8000 + i
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
# tab <- run_CSA_both_analysis(
#   data_dir = "CSA_data",
#   alpha = 0.05,
#   gamma_search_grid = seq(1, 5, by = 0.1),
#   output_csv = "PDF_images/CSA_both_results.csv",
#   n_cores = detect_available_cores()
# )
# print(tab)
