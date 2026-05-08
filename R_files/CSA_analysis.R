# CSA data analysis (RF only + expanded AIPW sample + weighted combination)
# --------------------------------------------------------------------------
# Changes vs real_data_analysis.R:
# 1) No SuperLearner; use random forests (ranger) for nuisance fits.
# 2) Cross-fit AIPW sample uses unmatched observations PLUS treated/control from complete pairs.
# 3) Combined p-value uses:
#    min{1, 2(w1*P1 + w2*P2)},
#    w1 = 2n / (4n + n1 + n2), w2 = (2n + n1 + n2) / (4n + n1 + n2)
#    where n = #complete pairs, n1 = #treated-only pairs, n2 = #control-only pairs.

source("R_files/combined_test_simulation_new.R")

assert_real_data_new_packages <- function() {
  needed <- c("ranger")
  missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required packages: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
}

read_covariate_table_new <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  id_col <- names(df)[2]
  cov_cols <- names(df)[3:7]
  out <- df[, c(id_col, cov_cols), drop = FALSE]
  names(out) <- c("id", make.names(cov_cols))
  out$id <- as.numeric(out$id)
  out <- out[!duplicated(out$id), , drop = FALSE]
  out
}

read_outcome_file_new <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  data.frame(
    y_t = as.numeric(df[[2]]),
    id_t = as.numeric(df[[3]]),
    id_c = as.numeric(df[[4]]),
    y_c = as.numeric(df[[5]])
  )
}

is_observed_outcome_new <- function(y) {
  !is.na(y) & y >= 0
}

build_analysis_data_new <- function(outcome_df, trt_cov, cont_cov) {
  m <- merge(outcome_df, trt_cov, by.x = "id_t", by.y = "id", all.x = TRUE, sort = FALSE)
  names(m)[5:9] <- paste0("t_", names(trt_cov)[-1])
  m <- merge(m, cont_cov, by.x = "id_c", by.y = "id", all.x = TRUE, sort = FALSE)
  names(m)[10:14] <- paste0("c_", names(cont_cov)[-1])

  has_t <- is_observed_outcome_new(m$y_t)
  has_c <- is_observed_outcome_new(m$y_c)
  complete_idx <- which(has_t & has_c)
  only_t_idx <- which(has_t & !has_c)
  only_c_idx <- which(!has_t & has_c)

  complete_pairs <- data.frame(
    Y_t = m$y_t[complete_idx],
    Y_c = m$y_c[complete_idx]
  )

  cov_base <- names(trt_cov)[-1]

  make_arm_df <- function(idx, arm = c("treated", "control")) {
    arm <- match.arg(arm)
    if (length(idx) == 0) {
      out <- data.frame(T = integer(0), Y = numeric(0))
      for (nm in cov_base) out[[nm]] <- numeric(0)
      return(out)
    }
    out <- data.frame(
      T = if (arm == "treated") rep(1, length(idx)) else rep(0, length(idx)),
      Y = if (arm == "treated") m$y_t[idx] else m$y_c[idx]
    )
    for (nm in cov_base) {
      out[[nm]] <- if (arm == "treated") m[[paste0("t_", nm)]][idx] else m[[paste0("c_", nm)]][idx]
    }
    out
  }

  # Unmatched-only sample (for counts n1, n2)
  unmatched_treated <- make_arm_df(only_t_idx, "treated")
  unmatched_control <- make_arm_df(only_c_idx, "control")
  unmatched_incomplete_data <- rbind(unmatched_treated, unmatched_control)
  rownames(unmatched_incomplete_data) <- NULL

  # Expanded AIPW sample: complete-pair arms + unmatched arms
  complete_treated <- make_arm_df(complete_idx, "treated")
  complete_control <- make_arm_df(complete_idx, "control")
  aipw_data <- rbind(complete_treated, complete_control, unmatched_treated, unmatched_control)
  rownames(aipw_data) <- NULL

  list(
    complete_pairs = complete_pairs,
    aipw_data = aipw_data,
    counts = c(
      n = nrow(complete_pairs),
      n1 = nrow(unmatched_treated),
      n2 = nrow(unmatched_control),
      total_pairs = nrow(outcome_df)
    )
  )
}

fit_rf_regression_new <- function(df, y_col, x_cols, num_trees = 500, min_node_size = 5) {
  x_cols <- nonconstant_cols(df, x_cols)
  if (length(x_cols) == 0 || nrow(df) < 10) {
    return(list(type = "constant", mean = mean(df[[y_col]], na.rm = TRUE)))
  }
  fit_df <- df[, c(y_col, x_cols), drop = FALSE]
  fit <- ranger::ranger(
    formula = stats::as.formula(sprintf("%s ~ .", y_col)),
    data = fit_df,
    num.trees = num_trees,
    min.node.size = min_node_size,
    seed = 1
  )
  list(type = "ranger", fit = fit, x_cols = x_cols)
}

predict_rf_regression_new <- function(obj, newdata) {
  if (obj$type == "constant") return(rep(obj$mean, nrow(newdata)))
  as.numeric(stats::predict(obj$fit, data = newdata[, obj$x_cols, drop = FALSE])$predictions)
}

fit_propensity_rf_new <- function(df, x_cols, num_trees = 500, min_node_size = 5) {
  x_cols <- nonconstant_cols(df, x_cols)
  if (length(x_cols) == 0 || nrow(df) < 10) {
    return(list(type = "constant", p = mean(df$T == 1)))
  }
  fit_df <- df[, c("T", x_cols), drop = FALSE]
  fit_df$T <- factor(fit_df$T, levels = c(0, 1))
  fit <- ranger::ranger(
    formula = T ~ .,
    data = fit_df,
    probability = TRUE,
    num.trees = num_trees,
    min.node.size = min_node_size,
    seed = 2
  )
  list(type = "ranger", fit = fit, x_cols = x_cols)
}

predict_propensity_rf_new <- function(obj, newdata, clip_eps = 1e-6) {
  if (obj$type == "constant") return(rep(clip_prob(obj$p, eps = clip_eps), nrow(newdata)))
  pred <- stats::predict(obj$fit, data = newdata[, obj$x_cols, drop = FALSE])$predictions
  p1 <- as.numeric(pred[, "1"])
  clip_prob(p1, eps = clip_eps)
}

fit_theta_nu_models_rf_new <- function(train_arm, x_cols, kappa_neg, seed = NULL) {
  if (nrow(train_arm) < 6) {
    theta_fit <- fit_rf_regression_new(train_arm, "Y", x_cols)
    theta_pred <- predict_rf_regression_new(theta_fit, train_arm)
    train_arm$nu_target <- 1 + (kappa_neg - 1) * as.numeric(train_arm$Y < theta_pred)
    nu_fit <- fit_rf_regression_new(train_arm, "nu_target", x_cols)
    return(list(theta_fit = theta_fit, nu_fit = nu_fit))
  }
  split <- split_in_half(nrow(train_arm), seed = seed)
  df_theta <- train_arm[split$first, , drop = FALSE]
  df_nu <- train_arm[split$second, , drop = FALSE]
  if (nrow(df_nu) < 3) df_nu <- df_theta

  theta_fit <- fit_rf_regression_new(df_theta, "Y", x_cols)
  theta_on_nu <- predict_rf_regression_new(theta_fit, df_nu)
  df_nu$nu_target <- 1 + (kappa_neg - 1) * as.numeric(df_nu$Y < theta_on_nu)
  nu_fit <- fit_rf_regression_new(df_nu, "nu_target", x_cols)
  list(theta_fit = theta_fit, nu_fit = nu_fit)
}

crossfit_aipw_sensitivity_rf_new <- function(
  aipw_data,
  Gamma = 1,
  side = c("two.sided", "right", "left"),
  K = 5,
  x_cols,
  seed = NULL,
  clip_eps = 1e-6
) {
  side <- match.arg(side)
  if (nrow(aipw_data) < K) stop("Need at least K AIPW observations.", call. = FALSE)
  if (length(unique(aipw_data$T)) < 2) stop("AIPW data must include both arms.", call. = FALSE)

  x_cols <- nonconstant_cols(aipw_data, x_cols)
  N <- nrow(aipw_data)
  folds <- make_stratified_folds(aipw_data$T, K = K, seed = seed)

  ehat <- rep(NA_real_, N)
  theta1_L <- nu1_L <- theta0_L <- nu0_L <- rep(NA_real_, N)
  theta1_U <- nu1_U <- theta0_U <- nu0_U <- rep(NA_real_, N)

  for (k in seq_len(K)) {
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)
    train_df <- aipw_data[train_idx, , drop = FALSE]
    test_df <- aipw_data[test_idx, , drop = FALSE]
    train_t1 <- train_df[train_df$T == 1, , drop = FALSE]
    train_t0 <- train_df[train_df$T == 0, , drop = FALSE]

    prop_fit <- fit_propensity_rf_new(train_df, x_cols)
    ehat[test_idx] <- predict_propensity_rf_new(prop_fit, test_df, clip_eps = clip_eps)

    seed_base <- if (is.null(seed)) NULL else seed + 1000L * k
    fit_1L <- fit_theta_nu_models_rf_new(train_t1, x_cols, kappa_neg = Gamma, seed = seed_base)
    fit_0L <- fit_theta_nu_models_rf_new(train_t0, x_cols, kappa_neg = 1 / Gamma, seed = if (is.null(seed_base)) NULL else seed_base + 1L)
    fit_1U <- fit_theta_nu_models_rf_new(train_t1, x_cols, kappa_neg = 1 / Gamma, seed = if (is.null(seed_base)) NULL else seed_base + 2L)
    fit_0U <- fit_theta_nu_models_rf_new(train_t0, x_cols, kappa_neg = Gamma, seed = if (is.null(seed_base)) NULL else seed_base + 3L)

    theta1_L[test_idx] <- predict_rf_regression_new(fit_1L$theta_fit, test_df)
    theta0_L[test_idx] <- predict_rf_regression_new(fit_0L$theta_fit, test_df)
    theta1_U[test_idx] <- predict_rf_regression_new(fit_1U$theta_fit, test_df)
    theta0_U[test_idx] <- predict_rf_regression_new(fit_0U$theta_fit, test_df)
    nu1_L[test_idx] <- pmax(predict_rf_regression_new(fit_1L$nu_fit, test_df), clip_eps)
    nu0_L[test_idx] <- pmax(predict_rf_regression_new(fit_0L$nu_fit, test_df), clip_eps)
    nu1_U[test_idx] <- pmax(predict_rf_regression_new(fit_1U$nu_fit, test_df), clip_eps)
    nu0_U[test_idx] <- pmax(predict_rf_regression_new(fit_0U$nu_fit, test_df), clip_eps)
  }

  Tt <- aipw_data$T
  Yy <- aipw_data$Y
  ehat <- clip_prob(ehat, eps = clip_eps)

  r1L_pos <- pmax(Yy - theta1_L, 0)
  r1L_neg <- pmax(theta1_L - Yy, 0)
  r0L_pos <- pmax(Yy - theta0_L, 0)
  r0L_neg <- pmax(theta0_L - Yy, 0)
  r1U_pos <- pmax(Yy - theta1_U, 0)
  r1U_neg <- pmax(theta1_U - Yy, 0)
  r0U_pos <- pmax(Yy - theta0_U, 0)
  r0U_neg <- pmax(theta0_U - Yy, 0)

  phiL <-
    Tt * Yy +
    (1 - Tt) * theta1_L +
    Tt * ((r1L_pos - Gamma * r1L_neg) * (1 - ehat) / (nu1_L * ehat)) -
    ((1 - Tt) * Yy + Tt * theta0_L +
       (1 - Tt) * ((r0L_pos - (1 / Gamma) * r0L_neg) * ehat / (nu0_L * (1 - ehat))))

  phiU <-
    Tt * Yy +
    (1 - Tt) * theta1_U +
    Tt * ((r1U_pos - (1 / Gamma) * r1U_neg) * (1 - ehat) / (nu1_U * ehat)) -
    ((1 - Tt) * Yy + Tt * theta0_U +
       (1 - Tt) * ((r0U_pos - Gamma * r0U_neg) * ehat / (nu0_U * (1 - ehat))))

  tauLhat <- mean(phiL)
  tauUhat <- mean(phiU)
  sigmaLhat <- sqrt(mean((phiL - tauLhat)^2))
  sigmaUhat <- sqrt(mean((phiU - tauUhat)^2))
  zL <- if (sigmaLhat > 0) sqrt(N) * tauLhat / sigmaLhat else 0
  zU <- if (sigmaUhat > 0) sqrt(N) * tauUhat / sigmaUhat else 0
  p_right <- 1 - stats::pnorm(zL)
  p_left <- stats::pnorm(zU)
  p_two <- min(1, 2 * min(p_right, p_left))
  p_value <- switch(side, "right" = p_right, "left" = p_left, "two.sided" = p_two)
  list(p_value = p_value, p_right = p_right, p_left = p_left, p_two = p_two)
}

combined_weighted_test_new <- function(
  complete_pairs,
  aipw_data,
  counts,
  alpha = 0.05,
  Gamma_M = 1,
  Gamma_A = 1,
  match_side = c("two.sided", "right", "left"),
  incomplete_side = c("two.sided", "right", "left"),
  K = 5,
  x_cols,
  seed = NULL
) {
  match_side <- match.arg(match_side)
  incomplete_side <- match.arg(incomplete_side)

  pM <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = Gamma_M,
    side = match_side,
    score = "huber",
    huber_c = 1.345,
    mc_reps = 0,
    seed = seed
  )
  pA <- crossfit_aipw_sensitivity_rf_new(
    aipw_data = aipw_data,
    Gamma = Gamma_A,
    side = incomplete_side,
    K = K,
    x_cols = x_cols,
    seed = if (is.null(seed)) NULL else seed + 99L
  )

  n <- counts[["n"]]
  n1 <- counts[["n1"]]
  n2 <- counts[["n2"]]
  denom <- n + n1 + n2
  w1 <- if (denom > 0) n / denom else 0.5
  w2 <- if (denom > 0) (n1 + n2) / denom else 0.5
  p_comb <- min(1, min(pM$p_value / w1, pA$p_value / w2))

  list(
    p_value = p_comb,
    p_match = pM,
    p_incomplete = pA
  )
}

find_min_gamma_accept_new <- function(eval_fun, alpha = 0.05, gamma_grid = seq(1, 5, by = 0.1)) {
  pvals <- vapply(gamma_grid, eval_fun, FUN.VALUE = numeric(1))
  idx <- which(pvals > alpha)
  if (length(idx) == 0) return(NA_real_)
  gamma_grid[min(idx)]
}

fmt_pval_decision_new <- function(p, alpha = 0.05, digits = 4) {
  d <- ifelse(p <= alpha, "R", "A")
  paste0(formatC(p, format = "f", digits = digits), " (", d, ")")
}

analyze_one_outcome_new <- function(
  outcome_name,
  outcome_file,
  direction,
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

  if (nrow(complete_pairs) == 0) stop(sprintf("Outcome '%s' has no complete pairs.", outcome_name), call. = FALSE)
  if (nrow(aipw_data) < K || length(unique(aipw_data$T)) < 2) stop(sprintf("Outcome '%s' has insufficient AIPW sample.", outcome_name), call. = FALSE)

  match_side <- if (direction == "trt>cont") "right" else "left"
  inc_side <- match_side
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

run_real_data_analysis_new <- function(
  data_dir = "CSA_data",
  alpha = 0.05,
  gamma_search_grid = seq(1, 5, by = 0.1),
  output_csv = "PDF_images/CSA_results.csv",
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
      seed = 7000 + i
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
# tab <- run_real_data_analysis_new(
#   data_dir = "CSA_data",
#   alpha = 0.05,
#   gamma_search_grid = seq(1, 5, by = 0.1),
#   output_csv = "PDF_images/CSA_results.csv",
#   n_cores = detect_available_cores()
# )
# print(tab)
