# Real data analysis on CSA_data
# ------------------------------
# Outputs a 6-column table for 17 outcomes with:
# 1) outcome
# 2) fraction of pairs with both outcomes observed
# 3) Matching p-value + (A/R)
# 4) Min Gamma accepted for Matching (or NA)
# 5) Proposed p-value + (A/R)
# 6) Min Gamma accepted for Proposed (or NA)
#
# The proposed test uses SuperLearner-based nuisance fitting for the incomplete component.

source("R_files/combined_test_simulation.R")

assert_real_data_packages <- function() {
  needed <- c("SuperLearner", "ranger", "xgboost", "nnet")
  missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      sprintf("Missing required packages: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

get_sl_library <- function() {
  list(
    c("SL.ranger", "sl_screen_all"),
    c("SL.nnet", "sl_screen_all"),
    c("SL.xgboost", "sl_screen_all")
  )
}

sl_screen_all <- function(Y, X, family, obsWeights, id, ...) {
  rep(TRUE, ncol(X))
}

initialize_superlearner_wrappers <- function() {
  for (nm in c("SL.ranger", "SL.nnet", "SL.xgboost")) {
    if (!exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
      assign(nm, get(nm, envir = asNamespace("SuperLearner")), envir = .GlobalEnv)
    }
  }
  if (!exists("sl_screen_all", envir = .GlobalEnv, inherits = FALSE)) {
    assign("sl_screen_all", sl_screen_all, envir = .GlobalEnv)
  }
}

read_covariate_table <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  # File has an index-like first column, then id, then 5 covariates.
  id_col <- names(df)[2]
  cov_cols <- names(df)[3:7]
  out <- df[, c(id_col, cov_cols), drop = FALSE]
  names(out) <- c("id", make.names(cov_cols))
  out$id <- as.numeric(out$id)
  out <- out[!duplicated(out$id), , drop = FALSE]
  out
}

read_outcome_file <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  # Format: col2 treated outcome, col3 treated id, col4 control id, col5 control outcome
  out <- data.frame(
    y_t = as.numeric(df[[2]]),
    id_t = as.numeric(df[[3]]),
    id_c = as.numeric(df[[4]]),
    y_c = as.numeric(df[[5]])
  )
  out
}

is_observed_outcome <- function(y) {
  !is.na(y) & y >= 0
}

build_analysis_data <- function(outcome_df, trt_cov, cont_cov) {
  m <- merge(outcome_df, trt_cov, by.x = "id_t", by.y = "id", all.x = TRUE, sort = FALSE)
  names(m)[5:9] <- paste0("t_", names(trt_cov)[-1])
  m <- merge(m, cont_cov, by.x = "id_c", by.y = "id", all.x = TRUE, sort = FALSE)
  names(m)[10:14] <- paste0("c_", names(cont_cov)[-1])

  has_t <- is_observed_outcome(m$y_t)
  has_c <- is_observed_outcome(m$y_c)

  complete_idx <- which(has_t & has_c)
  only_t_idx <- which(has_t & !has_c)
  only_c_idx <- which(!has_t & has_c)

  complete_pairs <- data.frame(
    Y_t = m$y_t[complete_idx],
    Y_c = m$y_c[complete_idx]
  )

  cov_base <- names(trt_cov)[-1]
  treated_only <- data.frame(
    T = rep(1, length(only_t_idx)),
    Y = m$y_t[only_t_idx]
  )
  control_only <- data.frame(
    T = rep(0, length(only_c_idx)),
    Y = m$y_c[only_c_idx]
  )
  for (nm in cov_base) {
    treated_only[[nm]] <- m[[paste0("t_", nm)]][only_t_idx]
    control_only[[nm]] <- m[[paste0("c_", nm)]][only_c_idx]
  }
  incomplete_data <- rbind(treated_only, control_only)
  rownames(incomplete_data) <- NULL

  list(
    complete_pairs = complete_pairs,
    incomplete_data = incomplete_data
  )
}

fit_rf_regression <- function(df, y_col, x_cols) {
  x_cols <- nonconstant_cols(df, x_cols)
  if (length(x_cols) == 0) {
    return(list(type = "constant", mean = mean(df[[y_col]], na.rm = TRUE)))
  }
  if (nrow(df) < 10) {
    return(list(type = "constant", mean = mean(df[[y_col]], na.rm = TRUE)))
  }
  fit <- SuperLearner::SuperLearner(
    Y = df[[y_col]],
    X = df[, x_cols, drop = FALSE],
    family = stats::gaussian(),
    SL.library = get_sl_library(),
    method = "method.NNLS",
    cvControl = list(V = 5)
  )
  list(type = "superlearner", fit = fit, x_cols = x_cols)
}

predict_rf_regression <- function(obj, newdata) {
  if (obj$type == "constant") {
    return(rep(obj$mean, nrow(newdata)))
  }
  pred <- SuperLearner::predict.SuperLearner(obj$fit, newdata = newdata[, obj$x_cols, drop = FALSE])$pred
  as.numeric(pred)
}

fit_propensity_rf <- function(df, x_cols) {
  x_cols <- nonconstant_cols(df, x_cols)
  if (length(x_cols) == 0) {
    p <- mean(df$T == 1)
    return(list(type = "constant", p = p))
  }
  if (nrow(df) < 10) {
    p <- mean(df$T == 1)
    return(list(type = "constant", p = p))
  }
  fit <- SuperLearner::SuperLearner(
    Y = as.numeric(df$T),
    X = df[, x_cols, drop = FALSE],
    family = stats::binomial(),
    SL.library = get_sl_library(),
    method = "method.NNLS",
    cvControl = list(V = 5)
  )
  list(type = "superlearner", fit = fit, x_cols = x_cols)
}

predict_propensity_rf <- function(obj, newdata, clip_eps = 1e-6) {
  if (obj$type == "constant") {
    return(rep(clip_prob(obj$p, clip_eps), nrow(newdata)))
  }
  p1 <- SuperLearner::predict.SuperLearner(obj$fit, newdata = newdata[, obj$x_cols, drop = FALSE])$pred
  p1 <- as.numeric(p1)
  clip_prob(p1, eps = clip_eps)
}

fit_theta_nu_models_rf <- function(train_arm, x_cols, kappa_neg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(train_arm)
  if (n < 6) {
    theta_fit <- fit_rf_regression(train_arm, "Y", x_cols)
    theta_pred <- predict_rf_regression(theta_fit, train_arm)
    train_arm$nu_target <- 1 + (kappa_neg - 1) * as.numeric(train_arm$Y < theta_pred)
    nu_fit <- fit_rf_regression(train_arm, "nu_target", x_cols)
    return(list(theta_fit = theta_fit, nu_fit = nu_fit))
  }

  split <- split_in_half(n, seed = seed)
  df_theta <- train_arm[split$first, , drop = FALSE]
  df_nu <- train_arm[split$second, , drop = FALSE]
  if (nrow(df_nu) < 3) df_nu <- df_theta

  theta_fit <- fit_rf_regression(df_theta, "Y", x_cols)
  theta_on_nu <- predict_rf_regression(theta_fit, df_nu)
  df_nu$nu_target <- 1 + (kappa_neg - 1) * as.numeric(df_nu$Y < theta_on_nu)
  nu_fit <- fit_rf_regression(df_nu, "nu_target", x_cols)
  list(theta_fit = theta_fit, nu_fit = nu_fit)
}

crossfit_aipw_sensitivity_rf <- function(
  incomplete_data,
  Gamma = 1,
  side = c("two.sided", "right", "left"),
  K = 5,
  x_cols,
  seed = NULL,
  clip_eps = 1e-6
) {
  side <- match.arg(side)
  if (nrow(incomplete_data) < K) stop("Need at least K incomplete observations.", call. = FALSE)
  if (length(unique(incomplete_data$T)) < 2) stop("Incomplete data must contain both arms.", call. = FALSE)

  x_cols <- nonconstant_cols(incomplete_data, x_cols)
  N <- nrow(incomplete_data)
  folds <- make_stratified_folds(incomplete_data$T, K = K, seed = seed)

  ehat <- rep(NA_real_, N)
  theta1_L <- nu1_L <- theta0_L <- nu0_L <- rep(NA_real_, N)
  theta1_U <- nu1_U <- theta0_U <- nu0_U <- rep(NA_real_, N)

  for (k in seq_len(K)) {
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)
    train_df <- incomplete_data[train_idx, , drop = FALSE]
    test_df <- incomplete_data[test_idx, , drop = FALSE]
    train_t1 <- train_df[train_df$T == 1, , drop = FALSE]
    train_t0 <- train_df[train_df$T == 0, , drop = FALSE]

    prop_fit <- fit_propensity_rf(train_df, x_cols)
    ehat[test_idx] <- predict_propensity_rf(prop_fit, test_df, clip_eps = clip_eps)

    seed_base <- if (is.null(seed)) NULL else seed + 1000L * k
    fit_1L <- fit_theta_nu_models_rf(train_t1, x_cols, kappa_neg = Gamma, seed = seed_base)
    fit_0L <- fit_theta_nu_models_rf(train_t0, x_cols, kappa_neg = 1 / Gamma, seed = if (is.null(seed_base)) NULL else seed_base + 1L)
    fit_1U <- fit_theta_nu_models_rf(train_t1, x_cols, kappa_neg = 1 / Gamma, seed = if (is.null(seed_base)) NULL else seed_base + 2L)
    fit_0U <- fit_theta_nu_models_rf(train_t0, x_cols, kappa_neg = Gamma, seed = if (is.null(seed_base)) NULL else seed_base + 3L)

    theta1_L[test_idx] <- predict_rf_regression(fit_1L$theta_fit, test_df)
    theta0_L[test_idx] <- predict_rf_regression(fit_0L$theta_fit, test_df)
    theta1_U[test_idx] <- predict_rf_regression(fit_1U$theta_fit, test_df)
    theta0_U[test_idx] <- predict_rf_regression(fit_0U$theta_fit, test_df)

    nu1_L[test_idx] <- pmax(predict_rf_regression(fit_1L$nu_fit, test_df), clip_eps)
    nu0_L[test_idx] <- pmax(predict_rf_regression(fit_0L$nu_fit, test_df), clip_eps)
    nu1_U[test_idx] <- pmax(predict_rf_regression(fit_1U$nu_fit, test_df), clip_eps)
    nu0_U[test_idx] <- pmax(predict_rf_regression(fit_0U$nu_fit, test_df), clip_eps)
  }

  Tt <- incomplete_data$T
  Yy <- incomplete_data$Y
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
    ((1 - Tt) * Yy +
       Tt * theta0_L +
       (1 - Tt) * ((r0L_pos - (1 / Gamma) * r0L_neg) * ehat / (nu0_L * (1 - ehat))))

  phiU <-
    Tt * Yy +
    (1 - Tt) * theta1_U +
    Tt * ((r1U_pos - (1 / Gamma) * r1U_neg) * (1 - ehat) / (nu1_U * ehat)) -
    ((1 - Tt) * Yy +
       Tt * theta0_U +
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

combined_partial_missing_test_rf <- function(
  complete_pairs,
  incomplete_data,
  alpha = 0.05,
  Gamma_M = 1,
  Gamma_A = 1,
  match_side = c("two.sided", "right", "left"),
  incomplete_side = c("two.sided", "right", "left"),
  K = 5,
  x_cols,
  score = c("huber", "sign", "identity"),
  huber_c = 1.345,
  seed = NULL
) {
  match_side <- match.arg(match_side)
  incomplete_side <- match.arg(incomplete_side)
  score <- match.arg(score)

  pM <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = Gamma_M,
    side = match_side,
    score = score,
    huber_c = huber_c,
    mc_reps = 0,
    seed = seed
  )
  pA <- crossfit_aipw_sensitivity_rf(
    incomplete_data = incomplete_data,
    Gamma = Gamma_A,
    side = incomplete_side,
    K = K,
    x_cols = x_cols,
    seed = if (is.null(seed)) NULL else seed + 99L
  )
  p_comb <- 1 - (1 - min(pM$p_value, pA$p_value))^2
  list(
    p_value = p_comb,
    reject = as.logical(p_comb <= alpha),
    p_match = pM,
    p_incomplete = pA
  )
}

find_min_gamma_accept <- function(eval_fun, alpha = 0.05, gamma_grid = seq(1, 5, by = 0.1)) {
  pvals <- vapply(gamma_grid, eval_fun, FUN.VALUE = numeric(1))
  idx <- which(pvals > alpha)
  if (length(idx) == 0) return(NA_real_)
  gamma_grid[min(idx)]
}

fmt_pval_with_decision <- function(p, alpha = 0.05, digits = 4) {
  decision <- ifelse(p <= alpha, "R", "A")
  paste0(formatC(p, format = "f", digits = digits), " (", decision, ")")
}

analyze_one_outcome <- function(
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
  out_df <- read_outcome_file(outcome_file)
  data_obj <- build_analysis_data(out_df, trt_cov, cont_cov)
  complete_pairs <- data_obj$complete_pairs
  incomplete_data <- data_obj$incomplete_data
  x_cols <- names(trt_cov)[-1]

  if (nrow(complete_pairs) == 0) {
    stop(sprintf("Outcome '%s' has no complete pairs.", outcome_name), call. = FALSE)
  }
  use_incomplete <- nrow(incomplete_data) >= K && length(unique(incomplete_data$T)) == 2
  frac_complete <- nrow(complete_pairs) / nrow(out_df)

  match_side <- if (direction == "trt>cont") "right" else "left"
  inc_side <- match_side

  match_res <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = 1,
    side = match_side,
    score = "huber",
    huber_c = 1.345,
    mc_reps = 0,
    seed = seed
  )
  prop_res <- if (use_incomplete) {
    combined_partial_missing_test_rf(
      complete_pairs = complete_pairs,
      incomplete_data = incomplete_data,
      alpha = alpha,
      Gamma_M = 1,
      Gamma_A = 1,
      match_side = match_side,
      incomplete_side = inc_side,
      K = K,
      x_cols = x_cols,
      score = "huber",
      huber_c = 1.345,
      seed = seed
    )
  } else {
    list(p_value = match_res$p_value)
  }

  match_gamma <- NA_real_
  if (match_res$p_value <= alpha) {
    match_gamma <- find_min_gamma_accept(
      eval_fun = function(g) {
        matched_sensitivity_pvalue(
          complete_pairs = complete_pairs,
          Gamma = g,
          side = match_side,
          score = "huber",
          huber_c = 1.345,
          mc_reps = 0
        )$p_value
      },
      alpha = alpha,
      gamma_grid = gamma_search_grid
    )
  }

  prop_gamma <- NA_real_
  if (prop_res$p_value <= alpha) {
    prop_gamma <- find_min_gamma_accept(
      eval_fun = function(g) {
        if (use_incomplete) {
          combined_partial_missing_test_rf(
            complete_pairs = complete_pairs,
            incomplete_data = incomplete_data,
            alpha = alpha,
            Gamma_M = g,
            Gamma_A = g,
            match_side = match_side,
            incomplete_side = inc_side,
            K = K,
            x_cols = x_cols,
            score = "huber",
            huber_c = 1.345
          )$p_value
        } else {
          matched_sensitivity_pvalue(
            complete_pairs = complete_pairs,
            Gamma = g,
            side = match_side,
            score = "huber",
            huber_c = 1.345,
            mc_reps = 0
          )$p_value
        }
      },
      alpha = alpha,
      gamma_grid = gamma_search_grid
    )
  }

  data.frame(
    Outcome = outcome_name,
    Complete_pair_fraction = formatC(frac_complete, format = "f", digits = 3),
    Matching = fmt_pval_with_decision(match_res$p_value, alpha = alpha),
    Matching_min_Gamma_accept = ifelse(is.na(match_gamma), "NA", formatC(match_gamma, format = "f", digits = 2)),
    Proposed = fmt_pval_with_decision(prop_res$p_value, alpha = alpha),
    Proposed_min_Gamma_accept = ifelse(is.na(prop_gamma), "NA", formatC(prop_gamma, format = "f", digits = 2)),
    stringsAsFactors = FALSE
  )
}

run_real_data_analysis <- function(
  data_dir = "CSA_data",
  alpha = 0.05,
  gamma_search_grid = seq(1, 5, by = 0.1),
  output_csv = "PDF_images/real_data_analysis_results.csv",
  n_cores = detect_available_cores()
) {
  assert_required_packages()
  assert_real_data_packages()
  initialize_superlearner_wrappers()

  trt_cov <- read_covariate_table(file.path(data_dir, "trt_var.csv"))
  cont_cov <- read_covariate_table(file.path(data_dir, "cont_var.csv"))

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
    analyze_one_outcome(
      outcome_name = specs$outcome[i],
      outcome_file = file.path(data_dir, specs$file[i]),
      direction = specs$direction[i],
      trt_cov = trt_cov,
      cont_cov = cont_cov,
      alpha = alpha,
      gamma_search_grid = gamma_search_grid,
      K = 5,
      seed = 5000 + i
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
# res_tab <- run_real_data_analysis(
#   data_dir = "CSA_data",
#   alpha = 0.05,
#   gamma_search_grid = seq(1, 5, by = 0.1),
#   output_csv = "PDF_images/real_data_analysis_results.csv"
# )
# print(res_tab)
