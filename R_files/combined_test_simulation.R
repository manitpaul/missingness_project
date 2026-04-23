# Combined matched-pair / cross-fit AIPW sensitivity test
# ------------------------------------------------------
# This script implements Algorithm 1 from the draft for data consisting of
#   (i) complete matched pairs and
#   (ii) incomplete observations coming from partially observed matched pairs.
#
# Main exported functions:
#   1. combined_partial_missing_test()
#   2. run_eta_simulation()
#   3. plot_rejection_rates()
#
# Required packages: optmatch, ggplot2

assert_required_packages <- function() {
  needed <- c("optmatch", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      sprintf(
        "The following packages are required but not installed: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

expit <- function(x) 1 / (1 + exp(-x))

detect_available_cores <- function() {
  n <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.finite(n) && !is.na(n) && n >= 1) return(as.integer(n))

  getconf_n <- suppressWarnings(tryCatch(
    as.integer(system("getconf _NPROCESSORS_ONLN", intern = TRUE)),
    error = function(e) NA_integer_
  ))
  if (is.finite(getconf_n) && !is.na(getconf_n) && getconf_n >= 1) return(getconf_n)

  n_env <- suppressWarnings(as.integer(Sys.getenv("NCPUS", unset = "1")))
  if (is.finite(n_env) && !is.na(n_env) && n_env >= 1) return(n_env)

  1L
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

nonconstant_cols <- function(df, x_cols) {
  keep <- vapply(x_cols, function(nm) {
    x <- df[[nm]]
    if (is.factor(x)) {
      nlevels(droplevels(x)) > 1
    } else {
      stats::var(x, na.rm = TRUE) > 0
    }
  }, FUN.VALUE = logical(1))
  x_cols[keep]
}

build_rhs_formula <- function(x_cols) {
  if (length(x_cols) == 0) {
    return("1")
  }
  paste(x_cols, collapse = " + ")
}

make_model_formula <- function(response, x_cols) {
  stats::as.formula(sprintf("%s ~ %s", response, build_rhs_formula(x_cols)))
}

make_stratified_folds <- function(treat, K = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(treat)
  folds <- integer(n)
  for (tt in sort(unique(treat))) {
    idx <- which(treat == tt)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_ids <- rep(seq_len(K), length.out = length(idx))
    folds[idx] <- fold_ids
  }
  folds
}

make_psi_fun <- function(score = c("huber", "sign", "identity"), huber_c = 1.345) {
  score <- match.arg(score)
  if (score == "huber") {
    return(function(u) pmin(u, huber_c))
  }
  if (score == "sign") {
    return(function(u) as.numeric(u > 0))
  }
  function(u) u
}

matched_sensitivity_pvalue <- function(
  complete_pairs,
  Gamma = 1,
  side = c("two.sided", "right", "left"),
  score = c("huber", "sign", "identity"),
  huber_c = 1.345,
  mc_reps = 20000,
  seed = NULL
) {
  side <- match.arg(side)
  score <- match.arg(score)
  if (!all(c("Y_t", "Y_c") %in% names(complete_pairs))) {
    stop("complete_pairs must contain columns 'Y_t' and 'Y_c'.", call. = FALSE)
  }
  D <- complete_pairs$Y_t - complete_pairs$Y_c
  n <- length(D)
  if (n < 1) stop("At least one complete pair is required.", call. = FALSE)

  h0 <- stats::median(abs(D), na.rm = TRUE)
  if (!is.finite(h0) || h0 <= 0) {
    h0 <- max(stats::mad(D, center = 0, constant = 1, na.rm = TRUE), 1)
  }
  psi_fun <- make_psi_fun(score = score, huber_c = huber_c)
  psi_k <- psi_fun(abs(D) / h0)
  # Use analytic normal approximation to avoid Monte Carlo resampling.
  T_M <- sum(sign(D) * psi_k)
  mu_psi <- ((Gamma - 1) / (Gamma + 1)) * sum(psi_k)
  sigma2_psi <- (4 * Gamma / (Gamma + 1)^2) * sum(psi_k^2)
  sigma_psi <- sqrt(max(sigma2_psi, 0))

  if (!is.finite(sigma_psi) || sigma_psi <= 0) {
    p_right <- as.numeric(mu_psi >= T_M)
    p_left <- as.numeric(mu_psi >= -T_M)
  } else {
    p_right <- stats::pnorm(T_M, mean = mu_psi, sd = sigma_psi, lower.tail = FALSE)
    p_left <- stats::pnorm(-T_M, mean = mu_psi, sd = sigma_psi, lower.tail = FALSE)
  }
  p_two   <- min(1, 2 * min(p_right, p_left))

  p_value <- switch(
    side,
    "right" = p_right,
    "left" = p_left,
    "two.sided" = p_two
  )

  list(
    p_value = p_value,
    p_right = p_right,
    p_left = p_left,
    p_two = p_two,
    T_M = T_M,
    psi_k = psi_k,
    h0 = h0,
    n_complete = n,
    Gamma = Gamma,
    side = side,
    score = score,
    mc_reps = mc_reps,
    mu_psi = mu_psi,
    sigma_psi = sigma_psi
  )
}

fit_linear_reg <- function(df, y_col, x_cols, weights = NULL) {
  x_cols <- nonconstant_cols(df, x_cols)
  formula <- make_model_formula(y_col, x_cols)
  if (is.null(weights)) {
    fit <- stats::lm(formula, data = df)
  } else {
    fit <- stats::lm(formula, data = df, weights = weights)
  }
  list(fit = fit, x_cols = x_cols, formula = formula)
}

predict_linear_reg <- function(obj, newdata) {
  as.numeric(stats::predict(obj$fit, newdata = newdata))
}

fit_propensity_glm <- function(df, t_col = "T", x_cols) {
  x_cols <- nonconstant_cols(df, x_cols)
  formula <- make_model_formula(t_col, x_cols)
  fit <- stats::glm(formula, data = df, family = stats::binomial())
  list(fit = fit, x_cols = x_cols, formula = formula)
}

predict_propensity_glm <- function(obj, newdata, clip_eps = 1e-6) {
  p <- as.numeric(stats::predict(obj$fit, newdata = newdata, type = "response"))
  clip_prob(p, eps = clip_eps)
}

fit_asymmetric_ls <- function(
  df,
  y_col,
  x_cols,
  kappa_neg,
  maxit = 100,
  tol = 1e-7,
  ridge = 1e-8
) {
  x_cols <- nonconstant_cols(df, x_cols)
  formula <- make_model_formula(y_col, x_cols)
  X <- stats::model.matrix(formula, data = df)
  y <- df[[y_col]]
  p <- NCOL(X)

  beta <- rep(0, p)
  XtX <- crossprod(X) + diag(ridge, p)
  Xty <- crossprod(X, y)
  beta <- tryCatch(
    solve(XtX, Xty),
    error = function(e) qr.solve(XtX, Xty)
  )

  for (iter in seq_len(maxit)) {
    r <- as.numeric(y - X %*% beta)
    w <- ifelse(r >= 0, 1, kappa_neg)
    WX <- X * sqrt(w)
    Wy <- y * sqrt(w)
    XtX_w <- crossprod(WX) + diag(ridge, p)
    Xty_w <- crossprod(WX, Wy)
    beta_new <- tryCatch(
      solve(XtX_w, Xty_w),
      error = function(e) qr.solve(XtX_w, Xty_w)
    )
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }

  list(beta = beta, formula = formula, x_cols = x_cols)
}

predict_asymmetric_ls <- function(obj, newdata) {
  X <- stats::model.matrix(obj$formula, data = newdata)
  as.numeric(X %*% obj$beta)
}

split_in_half <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- sample.int(n, size = n, replace = FALSE)
  m1 <- ceiling(n / 2)
  list(first = idx[seq_len(m1)], second = idx[-seq_len(m1)])
}

fit_theta_nu_models <- function(
  train_arm,
  x_cols,
  kappa_neg,
  seed = NULL,
  clip_eps = 1e-6,
  maxit = 100,
  tol = 1e-7
) {
  if (nrow(train_arm) < 6) {
    stop("Need at least 6 observations per arm in the outer-training sample.", call. = FALSE)
  }
  split <- split_in_half(nrow(train_arm), seed = seed)
  idx_theta <- split$first
  idx_nu <- split$second
  if (length(idx_nu) < 3) idx_nu <- idx_theta

  df_theta <- train_arm[idx_theta, , drop = FALSE]
  df_nu <- train_arm[idx_nu, , drop = FALSE]

  theta_fit <- fit_asymmetric_ls(
    df = df_theta,
    y_col = "Y",
    x_cols = x_cols,
    kappa_neg = kappa_neg,
    maxit = maxit,
    tol = tol
  )

  theta_on_nu <- predict_asymmetric_ls(theta_fit, df_nu)
  pseudo_nu <- 1 + (kappa_neg - 1) * as.numeric(df_nu$Y < theta_on_nu)
  df_nu_aug <- df_nu
  df_nu_aug$nu_target <- pseudo_nu
  nu_fit <- fit_linear_reg(df_nu_aug, y_col = "nu_target", x_cols = x_cols)

  list(theta_fit = theta_fit, nu_fit = nu_fit)
}

crossfit_aipw_sensitivity <- function(
  incomplete_data,
  Gamma = 1,
  side = c("two.sided", "right", "left"),
  K = 5,
  x_cols = NULL,
  seed = NULL,
  clip_eps = 1e-6,
  maxit = 100,
  tol = 1e-7
) {
  side <- match.arg(side)
  if (!all(c("Y", "T") %in% names(incomplete_data))) {
    stop("incomplete_data must contain columns 'Y' and 'T'.", call. = FALSE)
  }
  if (is.null(x_cols)) {
    x_cols <- setdiff(names(incomplete_data), c("Y", "T"))
  }
  x_cols <- nonconstant_cols(incomplete_data, x_cols)
  N <- nrow(incomplete_data)
  if (N < K) stop("Need at least K incomplete observations.", call. = FALSE)
  if (length(unique(incomplete_data$T)) < 2) stop("Incomplete data must contain both treatment arms.", call. = FALSE)

  folds <- make_stratified_folds(incomplete_data$T, K = K, seed = seed)

  ehat <- rep(NA_real_, N)
  theta1_L <- nu1_L <- theta0_L <- nu0_L <- rep(NA_real_, N)
  theta1_U <- nu1_U <- theta0_U <- nu0_U <- rep(NA_real_, N)

  for (k in seq_len(K)) {
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)
    train_df <- incomplete_data[train_idx, , drop = FALSE]
    test_df <- incomplete_data[test_idx, , drop = FALSE]

    prop_fit <- fit_propensity_glm(train_df, t_col = "T", x_cols = x_cols)
    ehat[test_idx] <- predict_propensity_glm(prop_fit, test_df, clip_eps = clip_eps)

    train_t1 <- train_df[train_df$T == 1, , drop = FALSE]
    train_t0 <- train_df[train_df$T == 0, , drop = FALSE]

    seed_base <- if (is.null(seed)) NULL else seed + 1000L * k

    fit_1L <- fit_theta_nu_models(
      train_arm = train_t1,
      x_cols = x_cols,
      kappa_neg = Gamma,
      seed = seed_base,
      clip_eps = clip_eps,
      maxit = maxit,
      tol = tol
    )
    fit_0L <- fit_theta_nu_models(
      train_arm = train_t0,
      x_cols = x_cols,
      kappa_neg = 1 / Gamma,
      seed = if (is.null(seed_base)) NULL else seed_base + 1L,
      clip_eps = clip_eps,
      maxit = maxit,
      tol = tol
    )
    fit_1U <- fit_theta_nu_models(
      train_arm = train_t1,
      x_cols = x_cols,
      kappa_neg = 1 / Gamma,
      seed = if (is.null(seed_base)) NULL else seed_base + 2L,
      clip_eps = clip_eps,
      maxit = maxit,
      tol = tol
    )
    fit_0U <- fit_theta_nu_models(
      train_arm = train_t0,
      x_cols = x_cols,
      kappa_neg = Gamma,
      seed = if (is.null(seed_base)) NULL else seed_base + 3L,
      clip_eps = clip_eps,
      maxit = maxit,
      tol = tol
    )

    theta1_L[test_idx] <- predict_asymmetric_ls(fit_1L$theta_fit, test_df)
    theta0_L[test_idx] <- predict_asymmetric_ls(fit_0L$theta_fit, test_df)
    theta1_U[test_idx] <- predict_asymmetric_ls(fit_1U$theta_fit, test_df)
    theta0_U[test_idx] <- predict_asymmetric_ls(fit_0U$theta_fit, test_df)

    nu1_L[test_idx] <- pmax(predict_linear_reg(fit_1L$nu_fit, test_df), clip_eps)
    nu0_L[test_idx] <- pmax(predict_linear_reg(fit_0L$nu_fit, test_df), clip_eps)
    nu1_U[test_idx] <- pmax(predict_linear_reg(fit_1U$nu_fit, test_df), clip_eps)
    nu0_U[test_idx] <- pmax(predict_linear_reg(fit_0U$nu_fit, test_df), clip_eps)
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
  p_left  <- stats::pnorm(zU)
  p_two   <- min(1, 2 * min(p_right, p_left))

  p_value <- switch(
    side,
    "right" = p_right,
    "left" = p_left,
    "two.sided" = p_two
  )

  list(
    p_value = p_value,
    p_right = p_right,
    p_left = p_left,
    p_two = p_two,
    tauLhat = tauLhat,
    tauUhat = tauUhat,
    sigmaLhat = sigmaLhat,
    sigmaUhat = sigmaUhat,
    phiL = phiL,
    phiU = phiU,
    ehat = ehat,
    Gamma = Gamma,
    side = side,
    K = K,
    x_cols = x_cols
  )
}

combined_partial_missing_test <- function(
  complete_pairs,
  incomplete_data,
  alpha = 0.05,
  Gamma_M = 1,
  Gamma_A = 1,
  match_side = c("two.sided", "right", "left"),
  incomplete_side = c("two.sided", "right", "left"),
  K = 5,
  x_cols = NULL,
  score = c("huber", "sign", "identity"),
  huber_c = 1.345,
  mc_reps_match = 20000,
  seed = NULL
) {
  match_side <- match.arg(match_side)
  incomplete_side <- match.arg(incomplete_side)
  score <- match.arg(score)
  if (is.null(x_cols)) {
    x_cols <- setdiff(names(incomplete_data), c("Y", "T"))
  }

  pM <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = Gamma_M,
    side = match_side,
    score = score,
    huber_c = huber_c,
    mc_reps = mc_reps_match,
    seed = seed
  )

  pA <- crossfit_aipw_sensitivity(
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
    alpha = alpha,
    Gamma_M = Gamma_M,
    Gamma_A = Gamma_A,
    p_match = pM,
    p_incomplete = pA,
    match_side = match_side,
    incomplete_side = incomplete_side
  )
}

simulate_population <- function(N = 2000, eta = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z1 <- rep(1, N)
  Z2 <- stats::rnorm(N, mean = 0, sd = 1)
  Z3 <- stats::runif(N, min = -0.5, max = 0.5)
  Z4 <- stats::rbinom(N, size = 1, prob = 0.5)
  Z5 <- stats::rnorm(N, mean = 0, sd = 4)
  Z6 <- stats::rbinom(N, size = 1, prob = 0.2)
  Z7 <- stats::rnorm(N, mean = 0, sd = 2)
  Z8 <- stats::rnorm(N, mean = 0, sd = 3)
  Z9 <- stats::rnorm(N, mean = 0, sd = 0.5)
  Z10 <- stats::rbinom(N, size = 1, prob = 0.8)
  Z11 <- stats::runif(N, min = -1, max = 1)

  Z <- data.frame(
    Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5, Z6 = Z6,
    Z7 = Z7, Z8 = Z8, Z9 = Z9, Z10 = Z10, Z11 = Z11
  )
  gamma <- c(-2, 0.15, 0.2, 0.1, 0.05, 0.3, 0.25, 0.08, 0.1, 0.75, 1.75)
  lin_ps <- as.numeric(as.matrix(Z) %*% gamma)
  p_treat <- expit(lin_ps)
  T <- stats::rbinom(N, size = 1, prob = p_treat)

  mean_y0 <- (2 / 11) * rowSums(Z) - 1
  Y0 <- stats::rnorm(N, mean = mean_y0, sd = 1)
  Y1 <- Y0 + eta
  Y_obs <- ifelse(T == 1, Y1, Y0)

  dat <- cbind.data.frame(id = seq_len(N), Z, T = T, Y0 = Y0, Y1 = Y1, Y = Y_obs, p_treat = p_treat)
  dat
}

build_optmatch_pairs <- function(pop_df, z_cols = paste0("Z", 2:11)) {
  assert_required_packages()
  dat <- pop_df
  dat$logit_ps <- stats::qlogis(clip_prob(dat$p_treat, eps = 1e-6))

  # Distances from logit propensity scores via match_on() and pairmatch().
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
    for (nm in c("Z1", z_cols)) {
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

induce_pair_level_missingness <- function(
  pair_df,
  gamma0_missing = -0.075,
  gamma_missing_rest = c(0.05, 0.02, 0.01, 0.05, 0.03, 0.25, 0.08, 0.01, 0.075, 0.175),
  seed = NULL,
  balance_categories = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  gamma_missing <- c(gamma0_missing, gamma_missing_rest)
  avg_cols <- paste0("avg_Z", 1:11)
  Zbar <- as.matrix(pair_df[, avg_cols, drop = FALSE])
  p_missing <- expit(as.numeric(Zbar %*% gamma_missing))
  OM <- stats::rbinom(nrow(pair_df), size = 1, prob = p_missing)
  B <- stats::rbinom(nrow(pair_df), size = 1, prob = 0.5)

  category <- ifelse(OM == 0, "complete", ifelse(B == 0, "treated_only", "control_only"))
  pair_df$OM <- OM
  pair_df$B <- B
  pair_df$category <- category

  if (balance_categories) {
    m <- min(sum(category == "complete"), sum(category == "treated_only"), sum(category == "control_only"))
    # Keep strict inequality n + n1 + n2 < N_init by design.
    max_balanced_m <- floor((nrow(pair_df) - 1) / 3)
    m <- min(m, max_balanced_m)
    if (m < 1) {
      stop("Could not form at least one balanced triplet of pair types. Try a different missingness intercept.", call. = FALSE)
    }
    idx_complete <- sample(which(category == "complete"), m)
    idx_tonly <- sample(which(category == "treated_only"), m)
    idx_conly <- sample(which(category == "control_only"), m)
  } else {
    idx_complete <- which(category == "complete")
    idx_tonly <- which(category == "treated_only")
    idx_conly <- which(category == "control_only")
    m <- length(idx_complete)
  }

  complete_pairs <- pair_df[idx_complete, c("pair_id", "Y_t", "Y_c"), drop = FALSE]

  treated_rows <- pair_df[idx_tonly, , drop = FALSE]
  control_rows <- pair_df[idx_conly, , drop = FALSE]

  treated_only <- data.frame(
    pair_id = treated_rows$pair_id,
    T = 1,
    Y = treated_rows$Y_t,
    Z1 = treated_rows$t_Z1,
    Z2 = treated_rows$t_Z2,
    Z3 = treated_rows$t_Z3,
    Z4 = treated_rows$t_Z4,
    Z5 = treated_rows$t_Z5,
    Z6 = treated_rows$t_Z6,
    Z7 = treated_rows$t_Z7,
    Z8 = treated_rows$t_Z8,
    Z9 = treated_rows$t_Z9,
    Z10 = treated_rows$t_Z10,
    Z11 = treated_rows$t_Z11
  )

  control_only <- data.frame(
    pair_id = control_rows$pair_id,
    T = 0,
    Y = control_rows$Y_c,
    Z1 = control_rows$c_Z1,
    Z2 = control_rows$c_Z2,
    Z3 = control_rows$c_Z3,
    Z4 = control_rows$c_Z4,
    Z5 = control_rows$c_Z5,
    Z6 = control_rows$c_Z6,
    Z7 = control_rows$c_Z7,
    Z8 = control_rows$c_Z8,
    Z9 = control_rows$c_Z9,
    Z10 = control_rows$c_Z10,
    Z11 = control_rows$c_Z11
  )

  incomplete_data <- rbind(treated_only, control_only)
  rownames(incomplete_data) <- NULL

  list(
    complete_pairs = complete_pairs,
    incomplete_data = incomplete_data,
    selected_pair_ids = c(pair_df$pair_id[idx_complete], pair_df$pair_id[idx_tonly], pair_df$pair_id[idx_conly]),
    counts = c(n = nrow(complete_pairs), n1 = nrow(treated_only), n2 = nrow(control_only)),
    missingness_prob = p_missing,
    categorized_pairs = pair_df
  )
}

run_one_experiment <- function(
  eta,
  N = 2000,
  alpha = 0.05,
  Gamma_M = 1,
  Gamma_A = 1,
  gamma0_missing = -0.075,
  match_side = "two.sided",
  incomplete_side = "two.sided",
  K = 5,
  score = "huber",
  huber_c = 1.345,
  mc_reps_match = 20000,
  seed = NULL
) {
  pop <- simulate_population(N = N, eta = eta, seed = seed)
  pairs <- build_optmatch_pairs(pop)
  missing_obj <- induce_pair_level_missingness(pairs, gamma0_missing = gamma0_missing, seed = if (is.null(seed)) NULL else seed + 1L)

  complete_pairs <- missing_obj$complete_pairs
  incomplete_data <- missing_obj$incomplete_data
  x_cols <- paste0("Z", 2:11)

  proposed <- combined_partial_missing_test(
    complete_pairs = complete_pairs,
    incomplete_data = incomplete_data,
    alpha = alpha,
    Gamma_M = Gamma_M,
    Gamma_A = Gamma_A,
    match_side = match_side,
    incomplete_side = incomplete_side,
    K = K,
    x_cols = x_cols,
    score = score,
    huber_c = huber_c,
    mc_reps_match = mc_reps_match,
    seed = if (is.null(seed)) NULL else seed + 2L
  )

  match_only <- matched_sensitivity_pvalue(
    complete_pairs = complete_pairs,
    Gamma = Gamma_M,
    side = match_side,
    score = score,
    huber_c = huber_c,
    mc_reps = mc_reps_match,
    seed = if (is.null(seed)) NULL else seed + 3L
  )

  list(
    eta = eta,
    proposed_p = proposed$p_value,
    proposed_reject = as.numeric(proposed$p_value <= alpha),
    match_p = match_only$p_value,
    match_reject = as.numeric(match_only$p_value <= alpha),
    counts = missing_obj$counts,
    proposed = proposed,
    match_only = match_only
  )
}

summarize_rejection_rates <- function(trial_df, conf_level = 0.90) {
  agg <- stats::aggregate(
    reject ~ eta + method,
    data = trial_df,
    FUN = function(x) c(rate = mean(x), n_rep = length(x))
  )
  summary_df <- data.frame(
    eta = agg$eta,
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

run_eta_simulation <- function(
  eta_grid = seq(-1, 1, by = 0.25),
  n_rep = 200,
  N = 2000,
  alpha = 0.05,
  Gamma_M = 1,
  Gamma_A = 1,
  gamma0_missing = -0.075,
  match_side = "two.sided",
  incomplete_side = "two.sided",
  K = 5,
  score = "huber",
  huber_c = 1.345,
  mc_reps_match = 20000,
  base_seed = 12345,
  n_cores = detect_available_cores(),
  show_progress = TRUE,
  save_path = NULL
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
    pct <- 100 * done_tasks / n_tasks
    cat(sprintf("Progress: %d/%d (%.1f%%)\n", done_tasks, n_tasks, pct))
    flush.console()
  }

  run_task <- function(task_id) {
    j <- grid_df$eta_idx[task_id]
    r <- grid_df$replication[task_id]
    eta <- eta_grid[j]
    seed_here <- base_seed + 100000L * j + r
    one <- run_one_experiment(
      eta = eta,
      N = N,
      alpha = alpha,
      Gamma_M = Gamma_M,
      Gamma_A = Gamma_A,
      gamma0_missing = gamma0_missing,
      match_side = match_side,
      incomplete_side = incomplete_side,
      K = K,
      score = score,
      huber_c = huber_c,
      mc_reps_match = mc_reps_match,
      seed = seed_here
    )
    data.frame(
      eta = eta,
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
        "grid_df", "eta_grid", "base_seed", "N", "alpha",
        "Gamma_M", "Gamma_A", "gamma0_missing", "match_side",
        "incomplete_side", "K", "score", "huber_c", "mc_reps_match",
        "run_one_experiment"
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
  summary_df <- summarize_rejection_rates(trial_df, conf_level = 0.90)

  out <- list(
    trials = trial_df,
    summary = summary_df,
    settings = list(
      eta_grid = eta_grid,
      n_rep = n_rep,
      N = N,
      alpha = alpha,
      Gamma_M = Gamma_M,
      Gamma_A = Gamma_A,
      gamma0_missing = gamma0_missing,
      match_side = match_side,
      incomplete_side = incomplete_side,
      K = K,
      score = score,
      huber_c = huber_c,
      mc_reps_match = mc_reps_match,
      base_seed = base_seed,
      n_cores = n_cores,
      show_progress = show_progress
    )
  )

  if (!is.null(save_path)) saveRDS(out, file = save_path)
  out
}

plot_rejection_rates <- function(
  results_obj,
  x_lab = expression(eta),
  y_lab = "Rejection rate",
  color_lab = "Test",
  palette = c("Proposed test" = "#1b9e77", "Matching test" = "#d95f02"),
  line_width = 1.4,
  point_size = 2.4,
  legend_position = "bottom",
  axis_title_size = 18,
  axis_text_size = 15,
  legend_title_size = 16,
  legend_text_size = 14,
  show_ci_band = TRUE,
  ci_alpha = 0.18,
  alpha_line = NULL,
  alpha_line_width = 1.2
) {
  assert_required_packages()
  if (is.list(results_obj) && !is.null(results_obj$summary)) {
    summary_df <- results_obj$summary
    if (show_ci_band && (!all(c("lower_ci", "upper_ci") %in% names(summary_df))) && !is.null(results_obj$trials)) {
      summary_df <- summarize_rejection_rates(results_obj$trials, conf_level = 0.90)
    }
  } else {
    summary_df <- results_obj
  }

  y_breaks <- NULL
  if (!is.null(alpha_line)) {
    y_vals <- summary_df$rejection_rate
    y_breaks <- sort(unique(c(stats::na.omit(y_vals), alpha_line)))
    y_breaks <- pretty(y_breaks, n = 6)
    y_breaks <- sort(unique(c(y_breaks, alpha_line)))
  }

  p <- ggplot2::ggplot(summary_df, ggplot2::aes_string(x = "eta", y = "rejection_rate", color = "method"))

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
    ggplot2::labs(x = x_lab, y = y_lab, color = color_lab) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
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

# Convenience wrapper for the requested simulation setup -------------------
run_requested_study <- function(
  eta_grid = seq(-1, 1, by = 0.25),
  n_rep = 200,
  N = 2000,
  alpha = 0.05,
  gamma0_missing = -0.075,
  K = 5,
  mc_reps_match = 20000,
  base_seed = 12345,
  n_cores = detect_available_cores(),
  show_progress = TRUE,
  save_path = "R_files/combined_test_results.rds"
) {
  out <- run_eta_simulation(
    eta_grid = eta_grid,
    n_rep = n_rep,
    N = N,
    alpha = alpha,
    Gamma_M = 1,
    Gamma_A = 1,
    gamma0_missing = gamma0_missing,
    match_side = "two.sided",
    incomplete_side = "two.sided",
    K = K,
    score = "huber",
    huber_c = 1.345,
    mc_reps_match = mc_reps_match,
    base_seed = base_seed,
    n_cores = n_cores,
    show_progress = show_progress,
    save_path = save_path
  )
  out
}

# Example usage ------------------------------------------------------------
#
# results <- run_requested_study(
#   eta_grid = seq(-1, 1, by = 0.25),
#   n_rep = 200,
#   N = 2000,
#   alpha = 0.05,
#   gamma0_missing = -0.075,
#   K = 5,
#   mc_reps_match = 20000,
#   base_seed = 12345,
#   n_cores = parallel::detectCores(logical = TRUE),
#   save_path = "R_files/combined_test_results.rds"
# )
#
# p <- plot_rejection_rates(results)
# print(p)
#
# ## Diagnostics can be changed later without rerunning the simulations:
# saved <- readRDS("R_files/combined_test_results.rds")
# p2 <- plot_rejection_rates(
#   saved,
#   x_lab = expression(eta),
#   y_lab = "Average rejection rate",
#   line_width = 1,
#   point_size = 2.2
# )
# print(p2)
