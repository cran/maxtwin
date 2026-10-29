test_that("maxtwin returns correct number of test points", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- sin(X[, 1]) + cos(X[, 2])

  for (gam in c(0.1, 0.2, 0.25, 0.5)) {
    idx <- maxtwin(X, Y, splitratio = gam)
    expect_equal(length(idx), round(nrow(X) * gam))
  }
})

test_that("initialmaxtwin returns correct number of test points", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- sin(X[, 1]) + cos(X[, 2])

  idx <- initialmaxtwin(X, Y, splitratio = 0.2)
  expect_equal(length(idx), 20L)
})

test_that("returned indices are valid row indices", {
  set.seed(42)
  X <- matrix(runif(300), ncol = 3)
  Y <- rowSums(X^2)

  idx <- maxtwin(X, Y, splitratio = 0.2)
  expect_true(all(idx >= 1L & idx <= nrow(X)))
  expect_equal(length(unique(idx)), length(idx))  # no duplicates
})

test_that("output is sorted integer vector", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(sin(X))

  idx <- maxtwin(X, Y, splitratio = 0.2)
  expect_true(is.numeric(idx))
  expect_true(!is.unsorted(idx))
})

test_that("deterministic output with same seed", {
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(X)

  set.seed(123)
  idx1 <- maxtwin(X, Y, splitratio = 0.2)
  set.seed(123)
  idx2 <- maxtwin(X, Y, splitratio = 0.2)
  expect_identical(idx1, idx2)
})

test_that("works with univariate input", {
  set.seed(42)
  X <- matrix(runif(100), ncol = 1)
  Y <- sin(X[, 1])

  idx <- maxtwin(X, Y, splitratio = 0.2)
  expect_equal(length(idx), 20L)
  expect_true(all(idx >= 1L & idx <= 100L))
})

test_that("works with high-dimensional input", {
  set.seed(42)
  X <- matrix(runif(160 * 8), ncol = 8)
  Y <- apply(X, 1, function(x) sum(sin(x * pi)))

  idx <- maxtwin(X, Y, splitratio = 0.2)
  expect_equal(length(idx), 32L)
})

test_that("factor response is supported", {
  idx <- maxtwin(
    X = iris[, 1:4],
    Y = iris[, 5],
    splitratio = 0.2
  )
  idx_initial <- initialmaxtwin(
    X = iris[, 1:4],
    Y = iris[, 5],
    splitratio = 0.2
  )

  expect_length(idx, 30L)
  expect_length(idx_initial, 30L)
  expect_equal(length(unique(idx)), 30L)
  expect_true(all(idx >= 1L & idx <= nrow(iris)))
})

test_that("factor predictors are supported", {
  X <- iris[, c("Sepal.Length", "Species")]
  Y <- iris$Sepal.Width

  idx <- maxtwin(X, Y, splitratio = 0.2)
  idx_initial <- initialmaxtwin(X, Y, splitratio = 0.2)

  expect_length(idx, 30L)
  expect_length(idx_initial, 30L)
  expect_equal(length(unique(idx)), 30L)
  expect_equal(length(unique(idx_initial)), 30L)
  expect_true(all(idx >= 1L & idx <= nrow(iris)))
})

test_that("a separate response is optional", {
  idx <- maxtwin(iris, splitratio = 0.2)
  idx_initial <- initialmaxtwin(iris, splitratio = 0.2)

  expect_length(idx, 30L)
  expect_length(idx_initial, 30L)
  expect_equal(length(unique(idx)), 30L)
  expect_equal(length(unique(idx_initial)), 30L)
  expect_true(all(idx >= 1L & idx <= nrow(iris)))
})

test_that("initialization option changes result", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(X)

  set.seed(99)
  idx1 <- maxtwin(X, Y, splitratio = 0.2, initialization = TRUE)
  set.seed(99)
  idx2 <- maxtwin(X, Y, splitratio = 0.2, initialization = FALSE)

  expect_equal(length(idx1), length(idx2))
  expect_true(all(idx1 >= 1L & idx1 <= nrow(X)))
  expect_true(all(idx2 >= 1L & idx2 <= nrow(X)))
})

test_that("separation distance is positive", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(sin(X))

  idx <- maxtwin(X, Y, splitratio = 0.2)
  train_idx <- setdiff(seq_len(nrow(X)), idx)

  sep_dists <- as.numeric(FNN::knnx.dist(
    X[train_idx, , drop = FALSE],
    X[idx, , drop = FALSE], k = 1
  ))
  expect_true(all(sep_dists > 0))
})

test_that("handles paired/clustered data without error", {
  set.seed(42)
  X_half <- matrix(runif(100), ncol = 2)
  X <- rbind(X_half, X_half + matrix(runif(100, -0.01, 0.01), ncol = 2))
  Y <- rowSums(sin(X))

  expect_no_error(idx <- maxtwin(X, Y, splitratio = 0.2))
  expect_equal(length(idx), 20L)
})

test_that("gam=0.5 produces half the data as test", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(X)

  idx <- maxtwin(X, Y, splitratio = 0.5)
  expect_equal(length(idx), 50L)
})

test_that("separation improves over random split for paired data", {
  set.seed(42)
  X_half <- matrix(runif(200), ncol = 2)
  X <- rbind(X_half, X_half + matrix(runif(200, -0.005, 0.005), ncol = 2))
  Y <- rowSums(sin(X))
  N <- nrow(X)

  idx_mt <- maxtwin(X, Y, splitratio = 0.2)

  set.seed(42)
  idx_rnd <- sort(sample(N, round(N * 0.2)))

  sep_mt <- mean(as.numeric(FNN::knnx.dist(
    X[-idx_mt, , drop = FALSE], X[idx_mt, , drop = FALSE], k = 1
  )))
  sep_rnd <- mean(as.numeric(FNN::knnx.dist(
    X[-idx_rnd, , drop = FALSE], X[idx_rnd, , drop = FALSE], k = 1
  )))

  expect_true(sep_mt > sep_rnd)
})

test_that("maxtwin matches original code results", {
  orig_path <- Sys.getenv("MAXTWIN_ORIGINAL_IMPL", unset = "")
  skip_if_not(nzchar(orig_path), "MAXTWIN_ORIGINAL_IMPL is not set")
  skip_if_not(file.exists(orig_path), "Original code not available")

  orig_env <- new.env(parent = globalenv())
  library(twinning)
  library(FNN)
  source(orig_path, local = orig_env)

  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(sin(X * pi))
  gam <- 0.2
  had_gam <- exists("gam", envir = globalenv(), inherits = FALSE)
  old_gam <- if (had_gam) get("gam", envir = globalenv()) else NULL
  on.exit({
    if (had_gam) {
      assign("gam", old_gam, envir = globalenv())
    } else if (exists("gam", envir = globalenv(), inherits = FALSE)) {
      rm("gam", envir = globalenv())
    }
  }, add = TRUE)
  assign("gam", gam, envir = globalenv())

  set.seed(99)
  idx_pkg <- maxtwin(X, Y, splitratio = gam, initialization = FALSE)
  set.seed(99)
  idx_orig <- orig_env$maxtwin(X, Y, splitratio = gam, initialization = FALSE)

  expect_equal(length(idx_pkg), length(idx_orig))

  sep_pkg <- mean(as.numeric(FNN::knnx.dist(
    X[-idx_pkg, ], X[idx_pkg, ], k = 1
  )))
  sep_orig <- mean(as.numeric(FNN::knnx.dist(
    X[-idx_orig, ], X[idx_orig, ], k = 1
  )))
  expect_true(abs(sep_pkg - sep_orig) / sep_orig < 0.1)
})

# ---- Input validation tests ----

test_that("invalid splitratio is rejected", {
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(X)

  expect_error(maxtwin(X, Y, splitratio = 0.3), "reciprocal")
  expect_error(maxtwin(X, Y, splitratio = 0.15), "reciprocal")
  expect_error(maxtwin(X, Y, splitratio = 0.6), "reciprocal")
  expect_error(maxtwin(X, Y, splitratio = 0.07), "reciprocal")
  expect_error(maxtwin(X, Y, splitratio = 0.4), "1/splitratio = 2.5")
})

test_that("Y length mismatch is rejected", {
  X <- matrix(runif(200), ncol = 2)
  Y <- runif(50)

  expect_error(maxtwin(X, Y, splitratio = 0.2), "length")
})

test_that("constant column in X is rejected", {
  X <- cbind(runif(100), rep(1, 100))
  Y <- runif(100)

  expect_error(maxtwin(X, Y, splitratio = 0.2), "constant")
})

test_that("constant Y is rejected", {
  X <- matrix(runif(200), ncol = 2)
  Y <- rep(5, 100)

  expect_error(maxtwin(X, Y, splitratio = 0.2), "constant")
})

test_that("too small N is rejected", {
  X <- matrix(runif(6), ncol = 2)
  Y <- runif(3)

  expect_error(maxtwin(X, Y, splitratio = 0.5), "at least")
})

test_that("result is independent of global gam variable", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(X)

  assign("gam", 0.99, envir = globalenv())
  idx <- maxtwin(X, Y, splitratio = 0.2)
  expect_equal(length(idx), 20L)
  rm("gam", envir = globalenv())
})

test_that("Ntest matches twinning output for various N", {
  for (N in c(13, 16, 17, 21, 50, 100)) {
    set.seed(42)
    X <- matrix(runif(N * 2), ncol = 2)
    Y <- rowSums(X)

    idx <- maxtwin(X, Y, splitratio = 0.2)
    # twinning returns ceiling(N / r) points
    expect_equal(length(idx), ceiling(N / 5),
                 label = paste("N =", N))
  }
})
