# Sensitivity analysis under unobserved covariates Z9, Z10, Z11
# --------------------------------------------------------------
# This script reuses core functions from combined_test_simulation_new.R and runs:
#   - eta = 0 (null)
#   - Gamma_M = Gamma_A = Gamma
#   - Gamma grid in {1, 3, 5, 9, 11}
# while creating matched pairs and fitting nuisance models using only observed
# covariates Z2..Z8 (Z9..Z11 treated as unobserved).

source("R_files/combined_test_simulation_new.R")

build_optmatch_pairs_partial_obs <- function(
  pop_df,
  match_z_cols = paste0("Z", 2:8),
  all_z_cols = paste0("Z", 2:11)
) {
  assert_required_packages()
  dat <- pop_df

  # Estimate propensity score with observed covariates only.
  ps_formula <- stats::as.formula(
    sprintf("T ~ %s", paste(match_z_cols, collapse = " + "))
  )
  ps_fit <- stats::glm(ps_formula, data = dat, family = stats::binomial())
  dat$logit_ps <- as.numeric(stats::predict(ps_fit, type = "link"))

  distmat <- optmatch::match_on(T ~ logit_ps, data = dat)
  pm <- optmatch::pairmatch(distmat, data = dat)
  dat$pair_id <- as.character(pm)
  dat <- dat[!is.na(dat$pair_id), , drop = FALSE]

  pair_ids <- unique(dat$pair_id)
  out_list <- vector("list", length(pair_ids))
  keep <- logical(length(pair_ids))

  for (idx in seq_along(pair_ids)) {
    pid <- pair_ids[idx]
    pair_df <- dat[dat$pair_id == pid, , drop = FALSE]
    if (nrow(pair_df) != 2 || !all(sort(pair_df$T) == c(0, 1))) next
    tr <- pair_df[pair_df$T == 1, , drop = FALSE]
    ct <- pair_df[pair_df$T == 0, , drop = FALSE]
    rec <- data.frame(
      pair_id = pid,
      treated_id = tr$id,
      control_id = ct$id,
      Y_t = tr$Y,
      Y_c = ct$Y,
      stringsAsFactors = FALSE
    )
    # Keep all Z columns in pair records for missingness generation.
    for (nm in c("Z1", all_z_cols)) {
      rec[[paste0("t_", nm)]] <- tr[[nm]]
      rec[[paste0("c_", nm)]] <- ct[[nm]]
      rec[[paste0("avg_", nm)]] <- mean(c(tr[[nm]], ct[[nm]]))
    }
    out_list[[idx]] <- rec
    keep[idx] <- TRUE
  }

  pair_df <- do.call(rbind, out_list[keep])
  rownames(pair_df) <- NULL
  pair_df
}

run_one_sensitivity_experiment <- function(
  Gamma,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  score = "huber",
  huber_c = 1.345,
  seed = NULL,
  match_z_cols = paste0("Z", 2:8),
  adjust_z_cols = paste0("Z", 2:8)
) {
  pop <- simulate_population(N = N, eta = 0, seed = seed)
  pairs <- build_optmatch_pairs_partial_obs(pop_df = pop, match_z_cols = match_z_cols)
  missing_obj <- induce_pair_level_missingness(
    pairs,
    gamma0_missing = gamma0_missing,
    seed = if (is.null(seed)) NULL else seed + 1L
  )

  complete_pairs <- missing_obj$complete_pairs
  incomplete_data <- missing_obj$incomplete_data

  proposed <- combined_partial_missing_test(
    complete_pairs = complete_pairs,
    incomplete_data = incomplete_data,
    alpha = alpha,
    Gamma_M = Gamma,
    Gamma_A = Gamma,
    match_side = "two.sided",
    incomplete_side = "two.sided",
    K = K,
    x_cols = adjust_z_cols,
    score = score,
    huber_c = huber_c,
    mc_reps_match = 0,
    seed = if (is.null(seed)) NULL else seed + 2L
  )

  match_only <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = Gamma,
    side = "two.sided",
    score = score,
    huber_c = huber_c,
    mc_reps = 0,
    seed = if (is.null(seed)) NULL else seed + 3L
  )

  list(
    Gamma = Gamma,
    proposed_p = proposed$p_value,
    proposed_reject = as.numeric(proposed$p_value <= alpha),
    match_p = match_only$p_value,
    match_reject = as.numeric(match_only$p_value <= alpha),
    counts = missing_obj$counts
  )
}

summarize_gamma_rejection_rates <- function(trial_df, conf_level = 0.90) {
  agg <- stats::aggregate(
    reject ~ Gamma + method,
    data = trial_df,
    FUN = function(x) c(rate = mean(x), n_rep = length(x))
  )
  summary_df <- data.frame(
    Gamma = agg$Gamma,
    method = agg$method,
    rejection_rate = agg$reject[, "rate"],
    n_rep = agg$reject[, "n_rep"]
  )
  z <- stats::qnorm(0.5 + conf_level / 2)
  se <- sqrt(summary_df$rejection_rate * (1 - summary_df$rejection_rate) / summary_df$n_rep)
  summary_df$conf_level <- conf_level
  summary_df$lower_ci <- pmax(0, summary_df$rejection_rate - z * se)
  summary_df$upper_ci <- pmin(1, summary_df$rejection_rate + z * se)
  summary_df
}

run_gamma_sensitivity_simulation <- function(
  gamma_grid = c(1, 3, 5, 9, 11),
  n_rep = 200,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  score = "huber",
  huber_c = 1.345,
  base_seed = 12345,
  n_cores = detect_available_cores(),
  show_progress = TRUE,
  save_path = "R_files/sensitivity_results.rds"
) {
  assert_required_packages()

  grid_df <- expand.grid(
    gamma_idx = seq_along(gamma_grid),
    replication = seq_len(n_rep),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n_tasks <- nrow(grid_df)
  n_cores <- as.integer(n_cores)
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  n_cores <- min(n_cores, n_tasks)
  batch_size <- max(1L, 2L * n_cores)
  records <- vector("list", n_tasks)

  report_progress <- function(done_tasks) {
    if (!isTRUE(show_progress)) return(invisible(NULL))
    pct <- 100 * done_tasks / n_tasks
    cat(sprintf("Progress: %d/%d (%.1f%%)\n", done_tasks, n_tasks, pct))
    flush.console()
  }

  run_task <- function(task_id) {
    j <- grid_df$gamma_idx[task_id]
    r <- grid_df$replication[task_id]
    Gamma <- gamma_grid[j]
    seed_here <- base_seed + 100000L * j + r
    one <- run_one_sensitivity_experiment(
      Gamma = Gamma,
      N = N,
      alpha = alpha,
      gamma0_missing = gamma0_missing,
      K = K,
      score = score,
      huber_c = huber_c,
      seed = seed_here
    )
    data.frame(
      Gamma = Gamma,
      replication = r,
      method = c("Proposed test", "Matching test"),
      p_value = c(one$proposed_p, one$match_p),
      reject = c(one$proposed_reject, one$match_reject),
      n_complete = one$counts[["n"]],
      n_treated_only = one$counts[["n1"]],
      n_control_only = one$counts[["n2"]]
    )
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    for (start_idx in seq.int(1L, n_tasks, by = batch_size)) {
      end_idx <- min(start_idx + batch_size - 1L, n_tasks)
      task_chunk <- seq.int(start_idx, end_idx)
      records[task_chunk] <- parallel::mclapply(
        X = task_chunk,
        FUN = run_task,
        mc.cores = n_cores,
        mc.preschedule = TRUE
      )
      report_progress(end_idx)
    }
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, iseed = base_seed)
    parallel::clusterExport(
      cl = cl,
      varlist = c(
        "grid_df", "gamma_grid", "base_seed", "N", "alpha", "gamma0_missing",
        "K", "score", "huber_c", "run_one_sensitivity_experiment"
      ),
      envir = environment()
    )
    for (start_idx in seq.int(1L, n_tasks, by = batch_size)) {
      end_idx <- min(start_idx + batch_size - 1L, n_tasks)
      task_chunk <- seq.int(start_idx, end_idx)
      records[task_chunk] <- parallel::parLapply(cl, task_chunk, run_task)
      report_progress(end_idx)
    }
  } else {
    for (task_id in seq_len(n_tasks)) {
      records[[task_id]] <- run_task(task_id)
      report_progress(task_id)
    }
  }

  trial_df <- do.call(rbind, records)
  summary_df <- summarize_gamma_rejection_rates(trial_df, conf_level = 0.90)

  out <- list(
    trials = trial_df,
    summary = summary_df,
    settings = list(
      gamma_grid = gamma_grid,
      n_rep = n_rep,
      N = N,
      alpha = alpha,
      gamma0_missing = gamma0_missing,
      K = K,
      score = score,
      huber_c = huber_c,
      base_seed = base_seed,
      n_cores = n_cores,
      show_progress = show_progress
    )
  )

  if (!is.null(save_path)) saveRDS(out, file = save_path)
  out
}

plot_sensitivity_rejection_rates <- function(
  results_obj,
  x_lab = expression(Gamma),
  y_lab = "Type-I error",
  color_lab = "Test",
  palette = c("Proposed test" = "#1b9e77", "Matching test" = "#d95f02"),
  line_width = 0.9,
  point_size = 2.0,
  legend_position = "bottom",
  axis_title_size = 18,
  axis_text_size = 15,
  legend_title_size = 16,
  legend_text_size = 14,
  show_ci_band = TRUE,
  ci_alpha = 0.2,
  alpha_line = 0.05,
  alpha_line_width = 0.6
) {
  assert_required_packages()
  if (is.list(results_obj) && !is.null(results_obj$summary)) {
    summary_df <- results_obj$summary
    if (show_ci_band && (!all(c("lower_ci", "upper_ci") %in% names(summary_df))) && !is.null(results_obj$trials)) {
      summary_df <- summarize_gamma_rejection_rates(results_obj$trials, conf_level = 0.90)
    }
  } else {
    summary_df <- results_obj
  }

  y_vals <- summary_df$rejection_rate
  y_breaks <- sort(unique(c(pretty(stats::na.omit(y_vals), n = 6), alpha_line)))

  p <- ggplot2::ggplot(summary_df, ggplot2::aes_string(x = "Gamma", y = "rejection_rate", color = "method"))

  if (show_ci_band && all(c("lower_ci", "upper_ci") %in% names(summary_df))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes_string(ymin = "lower_ci", ymax = "upper_ci", fill = "method"),
      alpha = ci_alpha,
      linewidth = 0,
      color = NA
    ) + ggplot2::scale_fill_manual(values = palette, guide = "none")
  }

  p <- p +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::labs(x = x_lab, y = y_lab, color = color_lab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position,
      axis.title = ggplot2::element_text(face = "bold", size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.title = ggplot2::element_text(face = "bold", size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      plot.title = ggplot2::element_blank()
    )

  if (!is.null(alpha_line)) {
    p <- p + ggplot2::geom_hline(
      yintercept = alpha_line,
      linetype = "dotted",
      linewidth = alpha_line_width,
      color = "black"
    )
  }

  p
}

run_sensitivity_study <- function(
  gamma_grid = c(1, 3, 5, 9, 11),
  n_rep = 200,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  base_seed = 12345,
  n_cores = detect_available_cores(),
  save_results_path = "R_files/sensitivity_results.rds",
  save_plot_path = "PDF_images/sensitivity_analysis.pdf"
) {
  results <- run_gamma_sensitivity_simulation(
    gamma_grid = gamma_grid,
    n_rep = n_rep,
    N = N,
    alpha = alpha,
    gamma0_missing = gamma0_missing,
    K = K,
    score = "huber",
    huber_c = 1.345,
    base_seed = base_seed,
    n_cores = n_cores,
    show_progress = TRUE,
    save_path = save_results_path
  )

  p <- plot_sensitivity_rejection_rates(
    results,
    x_lab = expression(Gamma),
    y_lab = "Type-I error",
    alpha_line = alpha
  )

  ggplot2::ggsave(
    filename = save_plot_path,
    plot = p,
    width = 10,
    height = 6,
    units = "in"
  )

  results
}

# Example usage:
# res <- run_sensitivity_study(
#   gamma_grid = c(1, 3, 5, 9, 11),
#   n_rep = 200,
#   N = 2000,
#   alpha = 0.05,
#   save_results_path = "R_files/sensitivity_results.rds",
#   save_plot_path = "PDF_images/sensitivity_analysis.pdf"
# )
