# Sensitivity analysis under favorable alternative:
# plot Gamma_crit(alpha) vs eta

source("R_files/combined_test_simulation_new.R")

find_min_gamma_accept <- function(p_values, gamma_grid, alpha = 0.05) {
  ord <- order(gamma_grid)
  g <- gamma_grid[ord]
  p <- p_values[ord]
  idx <- which(is.finite(p) & p > alpha)
  if (length(idx) == 0) return(NA_real_)
  g[min(idx)]
}

run_one_gamma_crit_experiment <- function(
  eta,
  gamma_grid,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  score = "huber",
  huber_c = 1.345,
  seed = NULL
) {
  pop <- simulate_population(N = N, eta = eta, seed = seed)
  pairs <- build_optmatch_pairs(pop)
  missing_obj <- induce_pair_level_missingness(
    pairs,
    gamma0_missing = gamma0_missing,
    seed = if (is.null(seed)) NULL else seed + 1L
  )

  complete_pairs <- missing_obj$complete_pairs
  aipw_data <- build_aipw_data_with_complete_pairs(missing_obj)
  x_cols <- paste0("Z", 2:11)
  n <- missing_obj$counts[["n"]]
  n1 <- missing_obj$counts[["n1"]]
  n2 <- missing_obj$counts[["n2"]]

  p_match <- numeric(length(gamma_grid))
  p_prop <- numeric(length(gamma_grid))
  for (j in seq_along(gamma_grid)) {
    g <- gamma_grid[j]
    match_only <- matched_sensitivity_pvalue(
      complete_pairs = complete_pairs,
      Gamma = g,
      side = "two.sided",
      score = score,
      huber_c = huber_c,
      mc_reps = 0,
      seed = if (is.null(seed)) NULL else seed + 100L + j
    )
    prop <- combined_partial_missing_test_new(
      complete_pairs = complete_pairs,
      aipw_data = aipw_data,
      n_complete = n,
      n1 = n1,
      n2 = n2,
      alpha = alpha,
      Gamma_M = g,
      Gamma_A = g,
      match_side = "two.sided",
      incomplete_side = "two.sided",
      K = K,
      x_cols = x_cols,
      score = score,
      huber_c = huber_c,
      seed = if (is.null(seed)) NULL else seed + 200L + j
    )
    p_match[j] <- match_only$p_value
    p_prop[j] <- prop$p_value
  }

  list(
    eta = eta,
    proposed_gamma_crit = find_min_gamma_accept(p_prop, gamma_grid, alpha = alpha),
    match_gamma_crit = find_min_gamma_accept(p_match, gamma_grid, alpha = alpha),
    counts = missing_obj$counts
  )
}

summarize_gamma_crit <- function(trial_df, conf_level = 0.90) {
  mean_agg <- stats::aggregate(
    gamma_crit ~ eta + method,
    data = trial_df,
    FUN = function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      mean(x)
    }
  )
  sd_agg <- stats::aggregate(
    gamma_crit ~ eta + method,
    data = trial_df,
    FUN = function(x) {
      x <- x[is.finite(x)]
      if (length(x) <= 1) return(NA_real_)
      stats::sd(x)
    }
  )
  n_rep_agg <- stats::aggregate(
    gamma_crit ~ eta + method,
    data = trial_df,
    FUN = length
  )
  n_finite_agg <- stats::aggregate(
    gamma_crit ~ eta + method,
    data = trial_df,
    FUN = function(x) sum(is.finite(x))
  )

  summary_df <- data.frame(
    eta = mean_agg$eta,
    method = mean_agg$method,
    mean_gamma_crit = mean_agg$gamma_crit,
    sd_gamma_crit = sd_agg$gamma_crit,
    n_rep = n_rep_agg$gamma_crit,
    n_finite = n_finite_agg$gamma_crit,
    stringsAsFactors = FALSE
  )
  z <- stats::qnorm(0.5 + conf_level / 2)
  se <- summary_df$sd_gamma_crit / sqrt(pmax(summary_df$n_finite, 1))
  summary_df$conf_level <- conf_level
  summary_df$lower_ci <- summary_df$mean_gamma_crit - z * se
  summary_df$upper_ci <- summary_df$mean_gamma_crit + z * se
  summary_df
}

run_eta_gamma_crit_simulation <- function(
  eta_grid = seq(-1, 1, by = 0.25),
  gamma_grid = seq(1, 5, by = 0.1),
  n_rep = 100,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  score = "huber",
  huber_c = 1.345,
  base_seed = 4242,
  n_cores = detect_available_cores(),
  show_progress = TRUE,
  save_path = "R_files/sensitivity_fav_alt_results.rds"
) {
  assert_required_packages()

  grid_df <- expand.grid(
    eta_idx = seq_along(eta_grid),
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
    cat(sprintf("Progress: %d/%d (%.1f%%)\n", done_tasks, n_tasks, 100 * done_tasks / n_tasks))
    flush.console()
  }

  run_task <- function(task_id) {
    j <- grid_df$eta_idx[task_id]
    r <- grid_df$replication[task_id]
    eta <- eta_grid[j]
    seed_here <- base_seed + 100000L * j + r
    one <- run_one_gamma_crit_experiment(
      eta = eta,
      gamma_grid = gamma_grid,
      N = N,
      alpha = alpha,
      gamma0_missing = gamma0_missing,
      K = K,
      score = score,
      huber_c = huber_c,
      seed = seed_here
    )
    data.frame(
      eta = eta,
      replication = r,
      method = c("Proposed test", "Matching test"),
      gamma_crit = c(one$proposed_gamma_crit, one$match_gamma_crit),
      n_complete = one$counts[["n"]],
      n_treated_only = one$counts[["n1"]],
      n_control_only = one$counts[["n2"]]
    )
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    for (start_idx in seq.int(1L, n_tasks, by = batch_size)) {
      end_idx <- min(start_idx + batch_size - 1L, n_tasks)
      task_chunk <- seq.int(start_idx, end_idx)
      records[task_chunk] <- parallel::mclapply(task_chunk, run_task, mc.cores = n_cores, mc.preschedule = TRUE)
      report_progress(end_idx)
    }
  } else {
    for (task_id in seq_len(n_tasks)) {
      records[[task_id]] <- run_task(task_id)
      report_progress(task_id)
    }
  }

  trial_df <- as.data.frame(do.call(rbind, records), stringsAsFactors = FALSE)
  summary_df <- summarize_gamma_crit(trial_df, conf_level = 0.90)

  out <- list(
    trials = trial_df,
    summary = summary_df,
    settings = list(
      eta_grid = eta_grid,
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

plot_gamma_crit_vs_eta <- function(
  results_obj,
  x_lab = expression(eta),
  y_lab = expression(Gamma[plain("crit")](alpha)),
  legend_title = "Test",
  line_width = 1.4,
  point_size = 2.2,
  ribbon_alpha = 0.15,
  axis_title_size = 15,
  axis_text_size = 13,
  legend_title_size = 14,
  legend_text_size = 12
) {
  df <- results_obj$summary
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = eta, y = mean_gamma_crit, color = method, fill = method, group = method)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
      alpha = ribbon_alpha,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      color = legend_title,
      fill = legend_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title = ggplot2::element_text(face = "bold", size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.title = ggplot2::element_text(face = "bold", size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size)
    )
  p
}

run_sensitivity_fav_alt_study <- function(
  eta_grid = seq(-1, 1, by = 0.25),
  gamma_grid = seq(1, 5, by = 0.1),
  n_rep = 100,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  base_seed = 4242,
  n_cores = detect_available_cores(),
  show_progress = TRUE,
  save_path = "R_files/sensitivity_fav_alt_results.rds"
) {
  run_eta_gamma_crit_simulation(
    eta_grid = eta_grid,
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
    show_progress = show_progress,
    save_path = save_path
  )
}

regenerate_sensitivity_fav_alt_plot <- function(
  cache_path = "R_files/sensitivity_fav_alt_results.rds",
  output_pdf = "PDF_images/sensitivity_fav_alt.pdf",
  width = 10,
  height = 6
) {
  results <- readRDS(cache_path)
  p <- plot_gamma_crit_vs_eta(results)
  ggplot2::ggsave(filename = output_pdf, plot = p, width = width, height = height, units = "in")
  invisible(p)
}

# Example usage:
# results <- run_sensitivity_fav_alt_study(
#   eta_grid = seq(-1, 1, by = 0.25),
#   gamma_grid = seq(1, 5, by = 0.1),
#   n_rep = 100,
#   n_cores = detect_available_cores(),
#   save_path = "R_files/sensitivity_fav_alt_results.rds"
# )
# p <- plot_gamma_crit_vs_eta(results)
# ggplot2::ggsave("PDF_images/sensitivity_fav_alt.pdf", p, width = 10, height = 6, units = "in")
# regenerate_sensitivity_fav_alt_plot(
#   cache_path = "R_files/sensitivity_fav_alt_results.rds",
#   output_pdf = "PDF_images/sensitivity_fav_alt.pdf"
# )
