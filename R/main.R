#' dymo
#'
#' @description DYMO is an R package for multi-feature time-series forecasting using Dynamic Mode Decomposition (DMD) on scaled returns. It reconstructs future levels by compounding predicted dynamics and applies conformal predictive sampling to generate calibrated uncertainty distributions. The result is a fast, interpretable, and distribution-free forecasting framework with fully empirical predictive intervals and ggplot visualizations.
#'
#' @param df A data frame or matrix of numeric time features, where each column represents a different time series of equal length. Must contain at least two features.
#' @param horizon Positive integer. Number of forecast steps ahead (horizon) to predict.
#' @param n_windows Positive integer. Number of validation windows used for residual and error estimation. Default: 10.
#' @param ci Numeric in (0, 1). Confidence level for the conformal predictive intervals. Default: 0.8.
#' @param smoother Logical. Whether to apply automatic LOESS smoothing to the input time series before modeling. Default: FALSE.
#' @param min_feats,max_feats Positive integers. Minimum and maximum number of features to combine in the multifeature DMD. By default both are set to the total number of columns in df.
#' @param dates Optional vector of `Date` objects (same length as the rows of df). If supplied, forecast plots will use date labels and forecast dates will be extrapolated automatically.
#' @param error_scale Character string specifying the scaling used in normalized error metrics. Options: "naive" (average one-step absolute error of the historical series) or "deviation" (standard deviation of the series). Default: `"naive"`.
#' @param error_benchmark Character string specifying the benchmark used for relative error metrics. Options: "naive" (last value persistence) or "average" (mean value). Default: "naive".
#' @param seed Positive integer. Random seed for reproducibility. Default: 42.
#' @param n_samp Positive integer. Number of conformal samples to draw for each feature and horizon. Default: 1000.
#' @param eig_max_mod Numeric in (0,1]. Maximum allowed eigenvalue modulus for DMD eigenmodes (values greater than this are clipped to prevent divergence). Default: 0.995.
#'
#' @return A list with three top-level components:
#'
#' \describe{
#'   \item{`comb_metrics`}{A data frame containing averaged forecast metrics and
#'     scores for all tested feature combinations and DMD ranks. Columns include
#'     rank, prediction scores, and standard error metrics such as MAE, MAPE,
#'     RMSSE, MASE, etc.}
#'
#'   \item{`best_model`}{A list summarizing the best-performing feature
#'     combination and DMD rank, containing:
#'       \describe{
#'         \item{`best_combination`}{Character string describing which features
#'           and rank were selected.}
#'         \item{`quant_preds`}{List of matrices (one per feature) summarizing
#'           forecast quantiles, mean, and standard deviation for each horizon.}
#'         \item{`point_forecast`}{Matrix of point (mean) forecasts on the
#'           original level scale. Rows = features, columns = forecast steps.}
#'         \item{`samples`}{List of matrices (one per feature) containing
#'           conformal resampled level trajectories (n_samp x horizon).}
#'         \item{`empfuns`}{List of empirical distribution functions per feature
#'           and horizon. Each element contains:
#'           \code{rfun(n)}, \code{pfun(x)}, \code{qfun(p)}, \code{dfun(x)}.}
#'         \item{`testing_errors`}{Vector of average forecast error metrics for
#'           the best rank and feature combination.}
#'         \item{`plots`}{List of ggplot objects visualizing the historical data,
#'           forecast mean, and predictive intervals for each feature.}
#'       }}
#'
#'   \item{`time_log`}{Character string giving the elapsed computation time in
#'     "hh:mm:ss" format.}
#' }
#'
#'
#' @importFrom stats approx density ecdf quantile loess fitted sd setNames
#' @importFrom utils combn head tail
#' @importFrom ggplot2 ggplot geom_line geom_ribbon aes xlab ylab theme_bw
#' @importFrom stats approx density ecdf quantile loess fitted sd
#' @importFrom utils head tail globalVariables
#' @importFrom ggplot2 ggplot geom_line geom_ribbon aes xlab ylab theme_bw
#' @importFrom rlang .data
#'
#' @examples
#'dymo(time_features[,c(2, 3, 4)], horizon = 10, dates = time_features$dates)
#'
#' @export
dymo <- function(df, horizon,
                 n_windows = 10, ci = 0.8, smoother = FALSE,
                 min_feats = NULL, max_feats = NULL, dates = NULL,
                 error_scale = "naive", error_benchmark = "naive",
                 seed = 42, n_samp = 1000, eig_max_mod = 0.995) {

  t0 <- proc.time()[["elapsed"]]
  set.seed(seed)

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n_feats <- ncol(df)
  if (n_feats < 2) stop("You need at least two time features.")

  if (is.null(min_feats)) min_feats <- n_feats
  if (is.null(max_feats)) max_feats <- n_feats
  min_feats <- max(2L, min(min_feats, n_feats))
  max_feats <- max(2L, min(max_feats, n_feats))

  if (anyNA(df)) df <- impute_df(df)
  if (isTRUE(smoother)) df <- smoother_fun(df)

  # Build global scaled transformers ONCE (consistent basis)
  global_transformers <- lapply(seq_len(ncol(df)), function(j) make_transformer_scaled(df[, j]))

  comb_idx <- all_combs(n_feats, min_feats, max_feats)

  comb_models <- lapply(comb_idx, function(ix)
    windower(df[, ix, drop = FALSE], horizon, n_windows, ci,
             error_scale, error_benchmark, dates,
             n_samp = n_samp, transformers = global_transformers[ix],
             eig_max_mod = eig_max_mod))

  comb_metrics <- do.call(rbind, lapply(seq_along(comb_models), function(i) {
    h <- comb_models[[i]]$history
    h$combined_features <- paste0(comb_idx[[i]], collapse = ",")
    h
  }))
  rownames(comb_metrics) <- NULL

  best_row <- which.max(comb_metrics$pred_scores)
  best_combination <- comb_metrics$combined_features[best_row]
  best_model <- comb_models[[which(vapply(comb_idx, function(ix) paste0(ix, collapse = ",") == best_combination, logical(1)))]]

  best_out <- list(
    best_combination = sprintf("combination %s, rank %s", best_combination, best_model$best_rank),
    quant_preds      = best_model$quant_preds,
    testing_errors   = best_model$testing_errors,
    point_forecast   = best_model$point_forecast,
    samples          = best_model$samples,
    empfuns          = best_model$empfuns,
    plots            = best_model$plots
  )

  elapsed <- proc.time()[["elapsed"]] - t0
  list(comb_metrics = comb_metrics,
       best_model   = best_out,
       time_log     = fmt_elapsed(elapsed))
}

####
# -------------- Windower + Conformal ---------

#' @keywords internal

windower <- function(df, horizon, n_windows = 10, ci = 0.8,
                     error_scale, error_benchmark, dates,
                     n_samp = 1000, transformers, eig_max_mod = 0.995) {

  feat_names <- colnames(df)
  n_length <- nrow(df); n_feats <- ncol(df)

  if (floor(n_length/(n_windows + 1)) <= horizon) stop("Not enough data for validation windows.")
  idx <- c(rep(1, n_length%%(n_windows + 1)),
           rep(seq_len(n_windows + 1), each = n_length/(n_windows + 1)))

  testing_errors_for_rank <- vector("list", n_feats)
  # per-feature residual banks: list(feature j) -> list(horizon h) -> numeric residuals
  horizon_resids_pf_allranks <- vector("list", n_feats)

  for (r in seq_len(n_feats)) {
    per_win_metrics <- vector("list", n_windows)
    horizon_resids_pf <- lapply(seq_len(ncol(df)), function(.) replicate(horizon, numeric(0), simplify = FALSE))
    pred_scores_acc <- c()

    for (w in seq_len(n_windows)) {
      df_train <- df[idx <= w, , drop = FALSE]
      df_hold  <- utils::head(df[idx == (w + 1), , drop = FALSE], horizon)
      if (nrow(df_hold) < horizon) next

      last_levels <- vapply(seq_len(ncol(df_train)), function(j) tail(df_train[, j], 1), numeric(1))

      dmd <- dynamic_mode_decomposition(df_train, rank = r, transformers = transformers,
                                        eig_max_mod = eig_max_mod)
      pf  <- dmd(horizon, last_levels)
      state_pf <- pf$state_forecast

      # Metrics on ORIGINAL scale
      split_forecasts <- split_along(pf$forecast, along = 1)
      split_holdouts  <- split_along(t(df_hold),   along = 1)
      x <- mapply(function(hh, ff, aa) my_metrics(holdout = hh, forecast = ff, actuals = aa,
                                                  error_scale = error_scale, error_benchmark = error_benchmark),
                  split_holdouts, split_forecasts, split_along(df_train, along = 2),
                  SIMPLIFY = FALSE)
      per_win_metrics[[w]] <- colMeans(do.call(rbind, x), na.rm = TRUE)

      # Build transform holdout state using the SAME transformers (scaled returns)
      state_hold <- matrix(NA_real_, ncol(df_train), horizon)
      for (j in seq_len(ncol(df_train))) {
        state_hold[j, ] <- transformers[[j]]$holdout_state(last_levels[j], df_hold[, j])
      }
      nr <- min(nrow(state_hold), nrow(state_pf)); nc <- min(ncol(state_hold), ncol(state_pf))
      if (nrow(state_hold) != nrow(state_pf) || ncol(state_hold) != ncol(state_pf)) {
        state_hold <- state_hold[seq_len(nr), seq_len(nc), drop = FALSE]
        state_pf   <- state_pf[seq_len(nr), seq_len(nc),   drop = FALSE]
      }
      state_resid <- state_hold - state_pf

      # Per-feature residual banks
      for (j in seq_len(ncol(df_train))) {
        for (h in seq_len(nc)) {
          horizon_resids_pf[[j]][[h]] <- c(horizon_resids_pf[[j]][[h]], state_resid[j, h])
        }
      }

      # Prediction score (original scale)
      ps <- mean(mapply(function(pred, hld) {
        prediction_score(matrix(pred, ncol = length(pred)), hld)
      }, split_forecasts, split_holdouts))
      pred_scores_acc <- c(pred_scores_acc, ps)
    }

    testing_errors_for_rank[[r]] <- c(pred_scores = mean(pred_scores_acc, na.rm = TRUE),
                                      colMeans(do.call(rbind, per_win_metrics), na.rm = TRUE))
    horizon_resids_pf_allranks[[r]] <- horizon_resids_pf
  }

  history <- as.data.frame(cbind(rank = seq_len(n_feats),
                                 do.call(rbind, lapply(testing_errors_for_rank, as.numeric))),
                           stringsAsFactors = FALSE)
  colnames(history) <- c("rank", names(testing_errors_for_rank[[1]]))
  best_rank <- history$rank[which.max(history$pred_scores)]
  resid_bank <- horizon_resids_pf_allranks[[best_rank]]   # list(feature j)[[horizon]]

  # Final model on full history, best rank (same transformers)
  last_levels_full <- vapply(seq_len(ncol(df)), function(j) tail(df[, j], 1), numeric(1))
  dmd_full <- dynamic_mode_decomposition(df, rank = best_rank, transformers = transformers,
                                         eig_max_mod = eig_max_mod)
  pf_full  <- dmd_full(horizon, last_levels_full)
  state_point  <- pf_full$state_forecast
  point_levels <- pf_full$forecast

  # Conformal sampling: joint (path-wise) per feature
  samples_by_feat <- lapply(seq_len(ncol(df)), function(j) {
    S <- matrix(NA_real_, n_samp, horizon)
    bank_j <- resid_bank[[j]]
    for (s in seq_len(n_samp)) {
      state_path <- numeric(horizon)
      for (h in seq_len(horizon)) {
        rb <- bank_j[[h]]; rb <- rb[is.finite(rb)]
        eps <- if (length(rb)) sample(rb, 1) else 0
        state_path[h] <- state_point[j, h] + eps
      }
      # returns → compound from last level to get levels path
      S[s, ] <- transformers[[j]]$inv_forecast(state_path, last_levels_full[j])
    }
    S
  })
  names(samples_by_feat) <- feat_names

  # Empirical distribution family per (feature, horizon)
  empfuns <- lapply(seq_len(ncol(df)), function(j) {
    lapply(seq_len(horizon), function(h) empfun(samples_by_feat[[j]][, h]))
  })
  names(empfuns) <- feat_names

  # Summary quantiles
  quants <- sort(unique(c((1 - ci)/2, 0.25, 0.5, 0.75, ci + (1 - ci)/2)))
  p_stats <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(rep(NA_real_, 10))
    qv <- stats::quantile(x, probs = quants, na.rm = TRUE)
    c(min = min(x), qv, max = max(x),
      mean = mean(x), sd = stats::sd(x))
  }
  quant_preds <- setNames(vector("list", ncol(df)), feat_names)
  for (j in seq_len(ncol(df))) {
    S <- samples_by_feat[[j]]
    Summ <- t(apply(S, 2, p_stats))
    rownames(Summ) <- paste0("t", seq_len(ncol(S)))
    quant_preds[[j]] <- Summ
  }

  # Dates & plotting
  if (is.null(dates)) {
    hist_dates   <- seq_len(n_length)
    forcat_dates <- (n_length + 1):(n_length + horizon)
  } else if (inherits(dates, "Date")) {
    hist_dates <- tail(dates, n_length)
    step <- as.integer(round(mean(diff(dates))))
    if (!is.finite(step) || step <= 0) step <- 1L
    forcat_dates <- seq(from = tail(dates, 1) + step, by = step, length.out = horizon)
  } else {
    hist_dates   <- seq_len(n_length)
    forcat_dates <- (n_length + 1):(n_length + horizon)
  }

  for (j in seq_len(ncol(df))) {
    nr <- nrow(quant_preds[[j]])
    if (!is.null(dates) && inherits(dates, "Date")) {
      rn <- if (length(forcat_dates) == nr) {
        as.character(forcat_dates)
      } else if (length(forcat_dates) > nr) {
        as.character(forcat_dates[seq_len(nr)])
      } else {
        as.character(rep(forcat_dates, length.out = nr))
      }
      rownames(quant_preds[[j]]) <- rn
    } else {
      rownames(quant_preds[[j]]) <- paste0("t", seq_len(nr))
    }
  }

  x_lab <- sprintf("Forecasting Horizon (n = %d)", horizon)
  y_lab <- sprintf("Forecasting Values for %s", feat_names)

  plots <- mapply(function(y_hist, Q, yl) {
    lower <- Q[, 2]; upper <- Q[, 6]
    ts_graph(x_hist = hist_dates, y_hist = y_hist,
             x_forcat = forcat_dates, y_forcat = Q[, "50%"],
             lower = lower, upper = upper,
             label_x = x_lab, label_y = yl)
  }, as.list(df), quant_preds, y_lab, SIMPLIFY = FALSE)

  list(best_rank = best_rank,
       history = history,
       testing_errors = testing_errors_for_rank[[best_rank]],
       quant_preds = quant_preds,         # summaries from conformal samples
       point_forecast = point_levels,     # features x horizon (original scale)
       samples = samples_by_feat,         # list of (n_samp x horizon) matrices
       empfuns = empfuns,                 # per feature -> per horizon -> {r,p,q,d}
       plots = plots)
}



# ---------------- Utilities -----------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Robust pseudo-inverse via SVD
pinv <- function(A, tol = NULL) {
  s <- svd(A)
  if (is.null(tol)) tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
  s$v %*% (diag(d_inv, nrow = length(d_inv))) %*% t(s$u)
}

# Split helper
split_along <- function(x, along = 1) {
  if (is.vector(x)) return(as.list(x))
  if (!is.matrix(x)) x <- as.matrix(x)
  if (along == 1) lapply(seq_len(nrow(x)), function(i) x[i, , drop = TRUE])
  else            lapply(seq_len(ncol(x)), function(j) x[, j, drop = TRUE])
}

# Time formatter
fmt_elapsed <- function(sec) {
  sec <- as.numeric(sec)
  h <- floor(sec/3600); m <- floor((sec - 3600*h)/60); s <- round(sec - 3600*h - 60*m)
  sprintf("%02dh:%02dm:%02ds", h, m, s)
}

# Simple imputer: linear interpolation + LOCF/BOCF
impute_lin_locf <- function(x) {
  x <- as.numeric(x)
  if (!anyNA(x)) return(x)
  idx <- seq_along(x)
  if (all(is.na(x))) return(rep(0, length(x)))
  approx(x = idx[!is.na(x)], y = x[!is.na(x)], xout = idx, method = "linear", rule = 2)$y
}
impute_df <- function(df) as.data.frame(lapply(df, impute_lin_locf), stringsAsFactors = FALSE)

# Loess smoother with AIC grid search (optional pre-smoothing of levels)
loess_smoother <- function(y, spans = c(0.15, 0.25, 0.35, 0.5, 0.65, 0.8)) {
  x <- seq_along(y); best <- NULL; best_aic <- Inf
  for (sp in spans) {
    fit <- try(loess(y ~ x, span = sp, family = "gaussian"), silent = TRUE)
    if (inherits(fit, "try-error")) next
    rss <- sum((y - fitted(fit))^2, na.rm = TRUE)
    k <- length(fit$coef %||% 2)
    aic <- length(y) * log(rss/length(y)) + 2 * k
    if (is.finite(aic) && aic < best_aic) { best_aic <- aic; best <- fit }
  }
  if (is.null(best)) y else fitted(best)
}
smoother_fun <- function(data) {
  if (is.vector(data)) data <- data.frame(v = data)
  as.data.frame(lapply(data, loess_smoother), stringsAsFactors = FALSE)
}

# ---------------- Metrics --------------------

ME   <- function(y, f) mean(f - y, na.rm = TRUE)
MAE  <- function(y, f) mean(abs(f - y), na.rm = TRUE)
MSE  <- function(y, f) mean((f - y)^2, na.rm = TRUE)
RMSE <- function(y, f) sqrt(MSE(y, f))

RMSSE <- function(y, f, scale) sqrt(mean(((f - y)/scale)^2, na.rm = TRUE))
MRE   <- function(y, f) mean((f - y)/y, na.rm = TRUE)
MPE   <- function(y, f) mean((y - f)/y * 100, na.rm = TRUE)
MAPE  <- function(y, f) mean(abs((y - f)/y) * 100, na.rm = TRUE)

rMAE <- function(y, f, b) mean(abs(f - y) / pmax(abs(b - y), .Machine$double.eps), na.rm = TRUE)
rRMSE <- function(y, f, b) {
  num <- mean((f - y)^2, na.rm = TRUE)
  den <- mean((b - y)^2, na.rm = TRUE)
  sqrt(num / pmax(den, .Machine$double.eps))
}
rAME <- function(y, f, b) mean((f - y) / pmax(abs(b - y), .Machine$double.eps), na.rm = TRUE)

MASE <- function(y, f, scale) MAE(y, f) / pmax(scale, .Machine$double.eps)
sMSE <- function(y, f, scale) MSE(y, f) / pmax(scale^2, .Machine$double.eps)
sCE  <- function(y, f, scale) sum(f - y, na.rm = TRUE) / pmax(scale * length(y), .Machine$double.eps)

GMRAE <- function(y, f, b) {
  r <- abs(f - y) / pmax(abs(b - y), .Machine$double.eps)
  exp(mean(log(pmax(r, .Machine$double.eps)), na.rm = TRUE))
}

my_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive") {
  scale <- switch(error_scale,
                  deviation = stats::sd(actuals, na.rm = TRUE),
                  naive     = mean(abs(diff(actuals)), na.rm = TRUE))
  bmk <- switch(error_benchmark,
                average = rep(mean(forecast, na.rm = TRUE), length(forecast)),
                naive   = rep(tail(actuals, 1), length(forecast)))
  out <- c(
    me    = ME(holdout, forecast),
    mae   = MAE(holdout, forecast),
    mse   = MSE(holdout, forecast),
    rmsse = RMSSE(holdout, forecast, scale),
    mpe   = MPE(holdout, forecast),
    mape  = MAPE(holdout, forecast),
    rmae  = rMAE(holdout, forecast, bmk),
    rrmse = rRMSE(holdout, forecast, bmk),
    rame  = rAME(holdout, forecast, bmk),
    mase  = MASE(holdout, forecast, scale),
    smse  = sMSE(holdout, forecast, scale),
    sce   = sCE(holdout, forecast, scale),
    gmrae = GMRAE(holdout, forecast, bmk)
  )
  round(out, 3)
}


# -------------- Forecast scoring -------------

prediction_score <- function(integrated_preds, ground_truth) {
  pfuns <- apply(integrated_preds, 2, stats::ecdf)
  pv <- mapply(function(Fn, y) Fn(y), pfuns, ground_truth)
  mean(1 - 2 * abs(pv - 0.5))
}

# -------------- Combinations -----------------

all_combs <- function(n, min_range = NULL, max_range = NULL) {
  if (is.null(min_range)) min_range <- 1L
  if (is.null(max_range)) max_range <- n
  out <- list()
  for (m in 1:n) out <- c(out, utils::combn(n, m, simplify = FALSE))
  out[vapply(out, function(ix) length(ix) >= min_range && length(ix) <= max_range, logical(1))]
}

# -------------- Plotting (ggplot2 only) ------

ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat,
                     lower = NULL, upper = NULL,
                     label_x = "Horizon", label_y= "Forecasted Var",
                     line_size = 1.3) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required. Please install.packages('ggplot2').")
  }
  all_data <- data.frame(x_all = c(x_hist, x_forcat),
                         y_all = c(y_hist, y_forcat))
  fc <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)
  if (!is.null(lower) && !is.null(upper)) { fc$lower <- lower; fc$upper <- upper }

  p <- ggplot() +
    geom_line(data = all_data, aes(x = .data$x_all, y = .data$y_all),
                 color = "gray40", linewidth = line_size)
  if (!is.null(lower) && !is.null(upper)) {
    p <- p + geom_ribbon(data = fc, aes(x = .data$x_forcat, ymin = .data$lower, ymax = .data$upper),
                            alpha = 0.3, fill = "darkseagreen1")
  }
  p + geom_line(data = fc, aes(x = .data$x_forcat, y = .data$y_forcat),
                   color = "darkseagreen4", linewidth = line_size) +
    xlab(label_x) + ylab(label_y) + theme_bw()
}

# --------- Scaled-returns transformer --------

# r_t = y_t / y_{t-1} - 1  (length T-1)
make_transformer_scaled <- function(x) {
  x <- as.numeric(x)
  fwd <- function(levels) {
    v <- as.numeric(levels)
    if (length(v) < 2) return(numeric(0))
    v[-1] / v[-length(v)] - 1
  }
  inv_state <- function(state, first_level) {
    r <- as.numeric(state)
    cumprod(c(1, 1 + r)) * first_level
  }
  inv_forecast <- function(state_forecast, last_level) {
    r <- as.numeric(state_forecast)    # horizon returns
    last_level * cumprod(1 + r)        # levels for horizons 1..H
  }
  holdout_state <- function(last_level, hold_levels) {
    v <- c(last_level, as.numeric(hold_levels))
    if (length(v) < 2) return(numeric(0))
    v[-1] / v[-length(v)] - 1
  }
  list(kind = "scaled", fwd = fwd, inv_state = inv_state,
       inv_forecast = inv_forecast, holdout_state = holdout_state, shift = 0)
}

# -------------- DMD on state (centered) ------

dynamic_mode_decomposition <- function(df, rank, transformers, eig_max_mod = 0.995) {
  # Build scaled-return states using provided transformers (consistent basis)
  state_rows <- lapply(seq_len(ncol(df)), function(j) transformers[[j]]$fwd(df[, j]))
  len <- unique(vapply(state_rows, length, integer(1)))
  if (length(len) != 1L) stop("Scaled transform produced unequal state lengths among features.")
  X <- t(do.call(rbind, state_rows))  # Tstate x features
  X <- t(X)                           # features x time
  if (ncol(X) < 2L) stop("Not enough state observations for DMD.")

  # Always center per-feature
  mu <- rowMeans(X)
  Xc <- X - mu

  X1 <- Xc[, -ncol(Xc), drop = FALSE]
  X2 <- Xc[, -1,        drop = FALSE]

  s <- svd(X1)
  r <- min(rank, length(s$d), nrow(X1), ncol(X1))
  if (r < 1L) stop("Rank too small in DMD.")

  Sr_inv <- diag(1 / s$d[1:r], r, r)
  A_tilde <- t(s$u[, 1:r, drop = FALSE]) %*% X2 %*% s$v[, 1:r, drop = FALSE] %*% Sr_inv
  eig <- eigen(A_tilde)

  # eigenvalue clipping
  lam <- eig$values
  mod <- Mod(lam)
  if (any(mod > eig_max_mod)) lam <- lam * (eig_max_mod / pmax(mod, .Machine$double.eps))

  psi <- X2 %*% s$v[, 1:r, drop = FALSE] %*% Sr_inv %*% eig$vectors
  Ahat <- Re(psi %*% diag(lam, r, r) %*% pinv(psi))

  function(h, last_levels_vec) {
    M <- matrix(0, nrow(X), h)
    m <- Xc[, ncol(Xc)]  # centered state
    for (t in seq_len(h)) { m <- Re(Ahat %*% m); M[, t] <- m }
    # de-center back to state space
    M <- M + mu
    # inverse to levels (compound from last observed level), feature by feature
    L <- matrix(NA_real_, nrow = nrow(M), ncol = ncol(M))
    for (j in seq_len(nrow(M))) {
      L[j, ] <- transformers[[j]]$inv_forecast(M[j, ], last_levels_vec[j])
    }
    list(state_forecast = M, forecast = L, transformers = transformers, mu = mu)
  }
}

# ------------- Empirical distribution --------

empfun <- function(sample_vec) {
  s <- as.numeric(sample_vec)
  s <- s[is.finite(s)]
  if (!length(s)) {
    rfun <- function(n) rep(NA_real_, n)
    pfun <- function(q) rep(NA_real_, length(q))
    qfun <- function(p) rep(NA_real_, length(p))
    dfun <- function(x) rep(NA_real_, length(x))
    return(list(rfun = rfun, pfun = pfun, qfun = qfun, dfun = dfun))
  }
  ec <- stats::ecdf(s)
  rfun <- function(n) sample(s, size = n, replace = TRUE)
  pfun <- function(q) ec(q)
  qfun <- function(p) stats::quantile(s, probs = p, type = 8, names = FALSE)
  kd  <- stats::density(s, n = max(512, min(2048, length(s) * 4)))
  dfun <- function(x) approx(kd$x, kd$y, xout = x, rule = 2, ties = "ordered")$y
  list(rfun = rfun, pfun = pfun, qfun = qfun, dfun = dfun)
}
