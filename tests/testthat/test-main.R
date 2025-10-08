# tests/testthat/test-temper.R
library(testthat)
library(dymo)

# Keep test runs deterministic and fast
set.seed(42)

# tiny toy generator used by multiple tests
make_toy_df <- function(n = 60) {
  t <- seq_len(n)
  x1 <- exp(0.01 * t) * (1 + 0.1 * sin(2*pi*t/12)) + rnorm(n, 0, 0.02)
  x2 <- 0.7 * x1 + 0.3 * exp(0.006 * t) + rnorm(n, 0, 0.03)
  x3 <- 0.5 * x1 + 0.4 * x2 + rnorm(n, 0, 0.02)
  data.frame(series_A = x1, series_B = x2, series_C = x3)
}


test_that("dymo basic run returns expected structure and sizes", {
  set.seed(123)
  df <- make_toy_df(60)
  dates <- as.Date("2024-01-01") + 0:(nrow(df)-1)

  res <- dymo(df,
              horizon    = 8,
              n_windows  = 3,
              ci         = 0.8,
              dates      = dates,
              n_samp     = 200,
              eig_max_mod = 0.99,
              min_feats  = 3,
              max_feats  = 3)

  # top-level names
  expect_true(all(c("comb_metrics", "best_model", "time_log") %in% names(res)))

  bm <- res$best_model
  expect_true(all(c("best_combination","quant_preds","point_forecast",
                    "samples","empfuns","testing_errors","plots") %in% names(bm)))

  # dims & types
  expect_equal(nrow(bm$point_forecast), 3L)
  expect_equal(ncol(bm$point_forecast), 8L)

  expect_length(bm$quant_preds, 3L)
  for (m in bm$quant_preds) {
    expect_true(is.matrix(m))
    expect_equal(nrow(m), 8L)
  }

  # samples size
  expect_length(bm$samples, 3L)
  for (S in bm$samples) {
    expect_equal(dim(S), c(200L, 8L))
    expect_true(all(is.finite(S)))
  }

  # empfuns behave sensibly
  emp_A_h3 <- bm$empfuns$series_A[[3]]
  qs <- emp_A_h3$qfun(c(0.1, 0.5, 0.9))
  expect_true(is.numeric(qs) && length(qs) == 3L)
  expect_true(qs[1] <= qs[2] && qs[2] <= qs[3])

  ps <- emp_A_h3$pfun(qs)
  expect_true(all(ps >= 0 & ps <= 1))

  # plots are ggplot objects
  expect_s3_class(bm$plots[[1]], "ggplot")
})

test_that("dymo results are reproducible with fixed seed", {
  df <- make_toy_df(50)
  dates <- as.Date("2024-06-01") + 0:(nrow(df)-1)

  res1 <- dymo(df, horizon = 6, n_windows = 3, ci = 0.8, dates = dates,
               n_samp = 150, eig_max_mod = 0.995, min_feats = 3, max_feats = 3, seed = 123)
  res2 <- dymo(df, horizon = 6, n_windows = 3, ci = 0.8, dates = dates,
               n_samp = 150, eig_max_mod = 0.995, min_feats = 3, max_feats = 3, seed = 123)

  # identical because sampling seed used inside
  expect_identical(res1$best_model$point_forecast, res2$best_model$point_forecast)
  # quant summaries should match exactly with same samples
  qp1 <- res1$best_model$quant_preds$series_A
  qp2 <- res2$best_model$quant_preds$series_A
  expect_identical(qp1, qp2)
})

test_that("errors: need at least two features", {
  df <- data.frame(only_one = rnorm(30))
  expect_error(dymo(df, horizon = 5), "at least two time features", ignore.case = TRUE)
})

test_that("errors: not enough data for validation windows", {
  # With n=20, n_windows=10, horizon=5 => floor(20/11) = 1 <= 5 -> error
  set.seed(1)
  df <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_error(dymo(df, horizon = 5, n_windows = 10, min_feats = 2, max_feats = 2),
               "not enough data", ignore.case = TRUE)
})


test_that("eigenvalue clipping prevents explosive forecasts", {
  set.seed(7)
  n <- 80
  t <- seq_len(n)
  x1 <- exp(0.02 * t) + rnorm(n, 0, 0.05)
  x2 <- 0.9 * x1 + rnorm(n, 0, 0.05)
  df <- data.frame(a = x1, b = x2, c = 0.5 * x1 + 0.3 * x2 + rnorm(n, 0, 0.03))

  res <- dymo(df, horizon = 5, n_windows = 4, ci = 0.8,
              n_samp = 100, eig_max_mod = 0.9, min_feats = 3, max_feats = 3, seed = 99)

  pf <- res$best_model$point_forecast
  expect_true(all(is.finite(pf)))

  # sanity bound: forecasts shouldn't exceed a large multiple of last levels
  last_levels <- vapply(df, function(x) tail(x, 1), numeric(1))
  bound <- 100 * max(last_levels)  # very loose
  expect_true(max(pf, na.rm = TRUE) < bound)
})


test_that("ggplot objects are built without error", {
  set.seed(2)
  df <- make_toy_df(40)
  res <- dymo(df, horizon = 6, n_windows = 3, ci = 0.8, n_samp = 100,
              min_feats = 3, max_feats = 3)
  p <- res$best_model$plots[[1]]
  expect_s3_class(p, "ggplot")
})

