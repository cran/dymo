#' dymo
#'
#' @param df A data frame with time features on columns. You need at least two time features. In case of missing values, automatic missing imputation through kalman filter will be performed.
#' @param seq_len Positive integer. Time-step number of the forecasting sequence. Default: NULL (automatic selection between 1 and the square root of full length).
#' @param n_windows Positive integer. Number of validation windows to test prediction error. Default: 10.
#' @param ci Confidence interval for prediction. Default: 0.8
#' @param smoother Logical. Flag to TRUE for loess smoothing. Default: FALSE.
#' @param min_feats Positive integer. Minimum number of time features to combine. Default: NULL (set equal to the total number of features)
#' @param max_feats Positive integer. Maximum number of time features to combine. Default: NULL (set equal to the total number of features)
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item comb_metrics: error metrics for all possible combinations of time features (for each combination, pred_score, me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae, are averaged across features, ranks and validation windows)
#' \item best_model: best combination resulting from the average prediction score across different ranks and features, including:
#' \itemize{
#' \item best_combination: combination of indexes and rank for the best model
#' \item testing_errors: testing errors for each time feature averaged across validation windows
#' \item quant_preds: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, median range ratio, upside probability and divergence for each point fo predicted sequences
#' \item plots: standard plot with confidence interval for each time feature
#' }
#' \item time_log
#' }
#'
#' @export
#'
#' @import purrr
#' @import tictoc
#' @import MASS
#' @import matlib
#' @import greybox
#' @import ggplot2
#' @importFrom scales number
#' @importFrom narray split
#' @importFrom readr parse_number
#' @importFrom lubridate seconds_to_period is.Date as.duration
#' @importFrom modeest mlv1
#' @importFrom moments kurtosis skewness
#' @import stats
#' @importFrom imputeTS na_kalman
#' @importFrom fANCOVA loess.as
#' @importFrom utils combn head tail


#'@examples
#'dymo(time_features[,c(2, 3, 4)], seq_len = 10, dates = time_features$dates)
#'


###
dymo <- function(df, seq_len, n_windows = 10, ci = 0.8, smoother = FALSE, min_feats = NULL, max_feats = NULL, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  tic("time")
  set.seed(seed)

  n_feats <- ncol(df)
  if(n_feats == 1){stop("You need at least two time features")}

  if(is.null(min_feats)){min_feats <- n_feats}
  if(is.null(max_feats)){max_feats <- n_feats}
  if(min_feats < 2){min_feats <- 2}
  if(min_feats > n_feats){min_feats <- n_feats}
  if(max_feats < 2){max_feats <- 2}
  if(max_feats > n_feats){max_feats <- n_feats}

  if(anyNA(df)){df <- na_kalman(df)}
  if(smoother == TRUE){df <- smoother(df)}

  comb_idx <- all_combs(n_feats, min_feats, max_feats)
  comb_models <- map(comb_idx, ~ windower(df[, .x, drop = FALSE], seq_len, n_windows = 10, ci = 0.8, error_scale, error_benchmark, dates))
  comb_metrics <- map(comb_models, ~ .x$history)
  comb_metrics <- Reduce(rbind, map2(comb_idx, comb_metrics, ~ cbind(combined_features = t(list(.x)), .y)))
  names(comb_models) <- as.character(unique(comb_metrics$combined_features))
  best_idx <- which.max(comb_metrics$pred_scores)
  best_combination <- as.character(comb_metrics$combined_features)[best_idx]

  rank_models <- flatten(comb_models[best_combination])
  quant_preds <- rank_models$quant_preds
  best_rank <- rank_models$best_rank
  testing_errors <- rank_models$testing_errors
  plots <- rank_models$plots
  best_combination <- paste0("combination ", best_combination, ", rank ", best_rank)
  best_model <- list(best_combination = best_combination, quant_preds = quant_preds, testing_errors = testing_errors, plots = plots)

  toc(log = TRUE)
  time_log<-seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(comb_metrics = comb_metrics, best_model = best_model, time_log = time_log)
  return(outcome)
}

###
windower <- function(df, seq_len, n_windows = 10, ci = 0.8, error_scale, error_benchmark, dates)
{
  deriv <- apply(df, 2, best_deriv)
  feat_names <- colnames(df)
  n_length <- nrow(df)
  n_feats <- ncol(df)

  too_short <- floor(n_length/(n_windows + 1)) <= seq_len
  if(too_short){stop("not enough data for the validation windows")}
  idx <- c(rep(1, n_length%%(n_windows + 1)), rep(1:(n_windows + 1), each = n_length/(n_windows + 1)))

  raw_errors_for_rank <- list()
  testing_errors_for_rank <- list()

  for(r in 1:n_feats)
  {
    pfuns <- map(1:n_windows, ~ dynamic_mode_decomposition(df[idx <= .x, , drop = FALSE], rank = r, deriv))
    seed_preds <- map(pfuns, ~ .x(seq_len))

    split_forecasts <- map(seed_preds, ~ split(.x, along = 1))
    split_holdouts <- map(1:n_windows, ~ split(head(df[idx == .x + 1,, drop = FALSE], seq_len), along = 2))
    split_actuals <- map(1:n_windows, ~ split(df[idx <= .x,, drop = FALSE], along = 2))

    testing_errors <- mapply(function(w) pmap(list(split_holdouts[[w]], split_forecasts[[w]], split_actuals[[w]]), ~ my_metrics(holdout = ..1, forecast = ..2, actuals = ..3, error_scale, error_benchmark)), w = 1:n_windows, SIMPLIFY = FALSE)
    testing_errors <- map(transpose(testing_errors), ~ colMeans(Reduce(rbind, .x)))
    testing_errors <- Reduce(rbind, testing_errors)
    raw_errors <- map2(1:n_windows, seed_preds, ~ t(head(df[idx == .x + 1,, drop = FALSE], seq_len)) - .y)
    collected_former_pred <- map(tail(seed_preds, -1), ~ as.data.frame(t(.x)))
    collected_former_errors <- transpose(map(head(raw_errors, -1), ~ as.data.frame(t(.x))))
    collected_former_errors <- map(collected_former_errors, ~ matrix(.x, 1))
    collected_former_errors <- map_depth(collected_former_errors, 2, ~ matrix(.x, 1, seq_len))
    collected_former_errors <- transpose(map(collected_former_errors, ~ accumulate(.x, rbind)))
    collected_former_holdouts <-  tail(map(1:n_windows, ~ head(df[idx == .x + 1,, drop = FALSE], seq_len)), -1)
    former_preds <- map(1:(n_windows - 1), ~ map2(collected_former_pred[[.x]], collected_former_errors[[.x]], ~ mapply(function(t) .x[t] + sample(.y[, t], size = 1000, replace = TRUE), t = 1:seq_len)))

    pred_scores <- colMeans(Reduce(rbind, map2(former_preds, collected_former_holdouts, ~ map2_dbl(.x, .y, ~ prediction_score(.x, .y)))))
    testing_errors <- cbind(pred_scores = pred_scores, testing_errors)
    rownames(testing_errors) <- feat_names

    testing_errors_for_rank[[r]] <- testing_errors
    raw_errors_for_rank[[r]] <- raw_errors
  }

  history <- as.data.frame(cbind(rank = 1:n_feats, Reduce(rbind, map(testing_errors_for_rank, ~ colMeans(.x)))))
  rownames(history) <- NULL
  best_rank <- history$rank[which.max(history$pred_scores)]
  best_raw <- raw_errors_for_rank[[best_rank]]

  raw_errors <- transpose(map(best_raw, ~ split(.x, along = 1)))
  collected_errors <- map(raw_errors, ~ Reduce(rbind, .x))
  preds <- dynamic_mode_decomposition(df, rank = best_rank, deriv)(seq_len)
  preds <- map2(split(preds, along = 1), collected_errors, ~ mapply(function(t) .x[t] + sample(.y[, t], size = 1000, replace = TRUE), t = 1:seq_len))
  preds <- map2(df, preds, ~ doxa_filter(.x, .y, n_class = NULL))

  quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))
  p_stats <- function(x){stats <- c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = tryCatch(suppressWarnings(modeest::mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(moments::kurtosis(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA), skewness = tryCatch(suppressWarnings(moments::skewness(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA)); return(stats)}
  quant_preds <- map(preds, ~ t(apply(.x, 2, p_stats)))
  names(quant_preds) <- feat_names

  iqr_to_range <- map(quant_preds, ~ tryCatch((.x[, "75%"] - .x[, "25%"])/(.x[, "max"] - .x[, "min"]), error = function(e) NA))
  median_range_ratio <- map(quant_preds, ~ tryCatch((.x[, "max"] - .x[, "50%"])/(.x[, "50%"] - .x[, "min"]), error = function(e) NA))
  upside_prob <- map(preds, ~ tryCatch(c(NA, colMeans(apply(.x[,-1]/.x[,-ncol(.x)], 2, function(x) x > 1))), error = function(e) NA))
  pvalues <- map(preds, ~ apply(.x, 2, function(x) ecdf(x)(seq(min(.x), max(.x), length.out = 100))), error = function(e) NA)
  divergence <- map(pvalues, ~ tryCatch(c(NA, apply(.x[,-1] - .x[,-ncol(.x)], 2, function(x) abs(max(x, na.rm = TRUE)))), error = function(e) NA))
  quant_preds <- pmap(list(quant_preds, iqr_to_range, median_range_ratio, upside_prob, divergence), ~ cbind(..1, iqr_to_range = ..2, median_range_ratio = ..3, upside_prob = ..4, divergence = ..5))

  if(is.null(dates)){hist_dates <- 1:n_length; forcat_dates <- (n_length + 1):(n_length + seq_len); quant_preds <- map(quant_preds, ~ {rownames(.x) <- paste0("t", 1:seq_len); return(.x)})}
  if(!is.null(dates) & is.Date(dates)){hist_dates <- tail(dates, n_length); forcat_dates <- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len); quant_preds <- map(quant_preds, ~ {rownames(.x) <- as.character(forcat_dates); return(.x)})}
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_lab <- paste0("Forecasting Values for ", feat_names)

  plots <- pmap(list(df, quant_preds, y_lab), ~ ts_graph(x_hist = hist_dates, y_hist = ..1, x_forcat = forcat_dates, y_forcat = ..2[,"50%"], lower = ..2[,2], upper = ..2[,6], label_x = x_lab, label_y = ..3))

  outcome <- list(best_rank = best_rank, history = history, testing_errors = testing_errors, quant_preds = quant_preds, plots = plots)
  return(outcome)
}

###
dynamic_mode_decomposition <- function(df, rank, deriv)
{
  df <- map2(df, max(deriv) - deriv, ~ smart_tail(.x, -.y))
  diff_models <- map2(df, deriv, ~ recursive_diff(.x, .y))

  if(sd(deriv) == 0){ddf <- as.data.frame(map(diff_models, ~ .x$vector))}
  if(sd(deriv) > 0){ddf <- as.data.frame(map(diff_models, ~ .x$vector))}
  tails <- map(diff_models, ~ .x$tail_value)
  mat <- as.matrix(t(ddf))

  mat1 <- mat[,-ncol(mat)]
  mat2 <- mat[,-1]
  svd_model <- svd(mat1)
  prox_A <- t(svd_model$u[, 1:rank]) %*% mat2 %*% svd_model$v[, 1:rank] %*% ginv(diag(svd_model$d[1:rank], rank, rank))
  eigen_A <- eigen(prox_A)
  psi <- mat2 %*% svd_model$v[, 1:rank] %*% ginv(diag(svd_model$d[1:rank], rank, rank)) %*% eigen_A$vectors
  A <- psi %*% diag(eigen_A$values, rank, rank) %*% ginv(psi)

  pred_fun <- function(n_steps)
  {
    p <- matrix(0, nrow(mat), n_steps)
    m <- mat[, ncol(mat)]
    for(n in 1:n_steps)
    {
      m <- Re(A %*% m)
      p[, n] <- m
    }

    p <- as.matrix(t(as.data.frame(map2(as.data.frame(t(p)), tails, ~ invdiff(.x, .y)))))
    return(p)
  }

  return(pred_fun)
}


###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "seagreen2", forcat_line = "seagreen3", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = c(y_hist, y_forcat))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- lower; forcat_data$upper <- upper}

  plot <- ggplot()+geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}

###
my_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive")
{
  scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
  benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))

  me <- ME(holdout, forecast, na.rm = TRUE)
  mae <- MAE(holdout, forecast, na.rm = TRUE)
  mse <- MSE(holdout, forecast, na.rm = TRUE)
  rmsse <- RMSSE(holdout, forecast, scale, na.rm = TRUE)
  mre <- MRE(holdout, forecast, na.rm = TRUE)
  mpe <- MPE(holdout, forecast, na.rm = TRUE)
  mape <- MAPE(holdout, forecast, na.rm = TRUE)
  rmae <- rMAE(holdout, forecast, benchmark, na.rm = TRUE)
  rrmse <- rRMSE(holdout, forecast, benchmark, na.rm = TRUE)
  rame <- rAME(holdout, forecast, benchmark, na.rm = TRUE)
  mase <- MASE(holdout, forecast, scale, na.rm = TRUE)
  smse <- sMSE(holdout, forecast, scale, na.rm = TRUE)
  sce <- sCE(holdout, forecast, scale, na.rm = TRUE)
  gmrae <- GMRAE(holdout, forecast, benchmark, na.rm = TRUE)

  out <- round(c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae), 3)
  return(out)
}

###
doxa_filter <- function(orig, mat, n_class = NULL)
{
  discrete_check <- all(orig%%1 == 0)
  all_positive_check <- all(orig >= 0)
  all_negative_check <- all(orig <= 0)
  monotonic_increase_check <- all(diff(orig) >= 0)
  monotonic_decrease_check <- all(diff(orig) <= 0)
  class_check <- FALSE
  if(is.integer(n_class)){class_check <- length(unique(orig)) <= n_class}

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(all_positive_check){mat[mat < 0] <- 0}
  if(all_negative_check){mat[mat > 0] <- 0}
  if(discrete_check){mat <- floor(mat)}
  if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
  if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
  if(class_check){mat[!(mat %in% unique(orig))] <- ((mat[!(mat %in% unique(orig))] > max(unique(orig))) * max(unique(orig))) + ((mat[!(mat %in% unique(orig))] < min(unique(orig))) * min(unique(orig)))}

  return(mat)
}

###
recursive_diff <- function(vector, deriv)
{
  vector <- unlist(vector)
  head_value <- vector("numeric", deriv)
  tail_value <- vector("numeric", deriv)
  if(deriv==0){head_value = NULL; tail_value = NULL}
  if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
  outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
  return(outcome)
}

###
invdiff <- function(vector, heads, add = FALSE)
{
  vector <- unlist(vector)
  if(is.null(heads)){return(vector)}
  for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
  if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
}

###
best_deriv <- function(ts, max_diff = 3, thresh = 0.001)
{
  pvalues <- vector(mode = "double", length = as.integer(max_diff))

  for(d in 1:(max_diff + 1))
  {
    model <- lm(ts ~ t, data.frame(ts, t = 1:length(ts)))
    pvalues[d] <- with(summary(model), pf(fstatistic[1], fstatistic[2], fstatistic[3],lower.tail=FALSE))
    ts <- diff(ts)
  }

  best <- tail(cumsum(pvalues < thresh), 1)

  return(best)
}

###
prediction_score <- function(integrated_preds, ground_truth)
{
  pfuns <- apply(integrated_preds, 2, ecdf)
  pvalues <- map2_dbl(pfuns, ground_truth, ~ .x(.y))
  scores <- mean(1 - 2 * abs(pvalues - 0.5))
  return(scores)
}

###
all_combs <- function(n, min_range = NULL, max_range = NULL)
{
  if(is.null(min_range)){min_range <- 1}
  if(is.null(max_range)){max_range <- n}
  comb_list <- mapply(function(m) split(combn(n,m), along = 2), m=1:n, SIMPLIFY = FALSE)
  comb_list <- flatten(comb_list)
  filtered_comb <- keep(comb_list, ~ length(.x) >= min_range)
  filtered_comb <- keep(filtered_comb, ~ length(.x) <= max_range)
  return(filtered_comb)
}

###
smoother <- function(data)
{
  if(is.vector(data)){data <- as.data.frame(data)}
  smoothed <- as.data.frame(purrr::map(data, ~ loess.as(x=1:length(.x), y=.x)$fitted))
  return(smoothed)
}

###
smart_head <- function(x, n)
{
  if(n != 0){return(head(x, n))}
  if(n == 0){return(x)}
}

###
smart_tail <- function(x, n)
{
  if(n != 0){return(tail(x, n))}
  if(n == 0){return(x)}
}
