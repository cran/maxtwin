test_that(".validate_inputs rejects unsupported X columns", {
  expect_error(
    maxtwin:::.validate_inputs(matrix("a", 10, 2), runif(10), 0.2),
    "numeric or factors"
  )
})

test_that(".validate_inputs accepts factor predictors and response", {
  X <- data.frame(
    value = seq_len(20),
    group = factor(rep(c("a", "b"), 10))
  )
  Y <- factor(rep(c("no", "yes"), each = 10))

  expect_no_error(maxtwin:::.validate_inputs(X, Y, 0.2))
  expect_no_error(maxtwin:::.validate_inputs(X, NULL, 0.2))
})

test_that(".encode_helmert matches Helmert coding", {
  group <- factor(c("b", "a", "c", "b"))
  encoded <- maxtwin:::.encode_helmert(data.frame(group = group))
  observed <- unique(group)
  expected <- stats::contr.helmert(length(observed))[
    match(group, observed),
    ,
    drop = FALSE
  ]

  expect_equal(unname(encoded), unname(expected))
  expect_equal(colnames(encoded), c("group.H1", "group.H2"))
})

test_that(".prepare_data uses encoded predictor dimension", {
  X <- iris[, c("Sepal.Length", "Species")]
  prepared <- maxtwin:::.prepare_data(X, iris$Species, 0.2)

  expect_equal(dim(prepared$X), c(150L, 3L))
  expect_equal(dim(prepared$Z), c(150L, 5L))
  expect_true(is.numeric(prepared$X))
  expect_true(is.numeric(prepared$Z))

  prepared_without_y <- maxtwin:::.prepare_data(X, NULL, 0.2)
  expect_equal(dim(prepared_without_y$X), c(150L, 3L))
  expect_equal(prepared_without_y$Z, prepared_without_y$X)
})

test_that(".validate_inputs rejects Y length mismatch", {
  expect_error(
    maxtwin:::.validate_inputs(matrix(runif(20), ncol = 2), runif(5), 0.2),
    "length"
  )
})

test_that(".validate_inputs rejects non-reciprocal splitratio", {
  X <- matrix(runif(200), ncol = 2)
  expect_error(maxtwin:::.validate_inputs(X, runif(100), 0.3), "reciprocal")
  expect_error(maxtwin:::.validate_inputs(X, runif(100), 0.7), "reciprocal")
  expect_error(
    maxtwin:::.validate_inputs(X, runif(100), 0.4),
    "1/splitratio = 2.5"
  )
})

test_that(".validate_inputs rejects malformed splitratio", {
  X <- matrix(runif(200), ncol = 2)
  Y <- runif(100)

  expect_error(maxtwin:::.validate_inputs(X, Y, "0.2"), "single finite positive")
  expect_error(maxtwin:::.validate_inputs(X, Y, c(0.2, 0.25)), "single finite positive")
  expect_error(maxtwin:::.validate_inputs(X, Y, NA_real_), "single finite positive")
  expect_error(maxtwin:::.validate_inputs(X, Y, NaN), "single finite positive")
  expect_error(maxtwin:::.validate_inputs(X, Y, Inf), "single finite positive")
  expect_error(maxtwin:::.validate_inputs(X, Y, 0), "single finite positive")
  expect_error(maxtwin:::.validate_inputs(X, Y, -0.2), "single finite positive")
  expect_error(maxtwin(X, Y), "single finite positive")
})

test_that(".validate_inputs accepts valid reciprocal splitratios", {
  X <- matrix(runif(200), ncol = 2)
  Y <- runif(100)
  for (sr in c(0.5, 1/3, 0.25, 0.2, 0.1)) {
    expect_no_error(maxtwin:::.validate_inputs(X, Y, sr))
  }
})

test_that(".validate_inputs rejects constant columns", {
  X <- cbind(runif(100), rep(0, 100))
  expect_error(maxtwin:::.validate_inputs(X, runif(100), 0.2), "constant")
})

test_that(".validate_inputs rejects invalid categorical data", {
  X_character <- data.frame(
    x = seq_len(20),
    group = rep(c("a", "b"), 10)
  )
  X_constant_factor <- data.frame(
    x = seq_len(20),
    group = factor(rep("a", 20))
  )

  expect_error(
    maxtwin:::.validate_inputs(X_character, runif(20), 0.2),
    "factor"
  )
  expect_error(
    maxtwin:::.validate_inputs(X_constant_factor, runif(20), 0.2),
    "constant"
  )
  expect_error(
    maxtwin:::.validate_inputs(
      data.frame(x = seq_len(20)),
      factor(rep("a", 20)),
      0.2
    ),
    "constant"
  )
})

test_that(".validate_inputs rejects missing and non-finite values", {
  X_missing <- data.frame(x = c(seq_len(19), NA_real_))
  X_infinite <- data.frame(x = c(seq_len(19), Inf))
  X <- data.frame(x = seq_len(20))

  expect_error(
    maxtwin:::.validate_inputs(X_missing, runif(20), 0.2),
    "missing"
  )
  expect_error(
    maxtwin:::.validate_inputs(X_infinite, runif(20), 0.2),
    "non-finite"
  )
  expect_error(
    maxtwin:::.validate_inputs(X, c(runif(19), Inf), 0.2),
    "non-finite"
  )
})

test_that(".run_twinning returns correct number of indices", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- sin(X[, 1])
  Z <- cbind(X, Y)

  idx <- maxtwin:::.run_twinning(Z, splitratio = 0.2,
                                  initialization = TRUE,
                                  standardization = TRUE)
  expect_equal(length(idx), 20L)
  expect_true(all(idx >= 1L & idx <= 100L))
})

test_that(".run_twinning works without initialization", {
  set.seed(42)
  X <- matrix(runif(200), ncol = 2)
  Y <- rowSums(X)
  Z <- cbind(X, Y)

  idx <- maxtwin:::.run_twinning(Z, splitratio = 0.2,
                                  initialization = FALSE,
                                  standardization = TRUE)
  expect_equal(length(idx), 20L)
})

test_that(".run_twinning applies the requested standardization", {
  set.seed(42)
  Z <- cbind(
    x = runif(100, 0, 1),
    y = runif(100, 0, 100)
  )
  u1 <- which.min(rowSums(sweep(Z, 2, colMeans(Z), "-")^2))

  observed_raw <- maxtwin:::.run_twinning(
    Z,
    splitratio = 0.2,
    initialization = TRUE,
    standardization = FALSE
  )
  expected_raw <- twinning::twin(
    Z,
    r = 5,
    u1 = u1,
    format_data = FALSE
  )
  observed_scaled <- maxtwin:::.run_twinning(
    Z,
    splitratio = 0.2,
    initialization = TRUE,
    standardization = TRUE
  )
  expected_scaled <- twinning::twin(
    scale(Z),
    r = 5,
    u1 = u1,
    format_data = FALSE
  )

  expect_identical(observed_raw, expected_raw)
  expect_identical(observed_scaled, expected_scaled)
})

test_that(".run_twinning validates logical options", {
  Z <- matrix(runif(200), ncol = 2)

  expect_error(
    maxtwin:::.run_twinning(Z, 0.2, NA, TRUE),
    "initialization"
  )
  expect_error(
    maxtwin:::.run_twinning(Z, 0.2, TRUE, "yes"),
    "standardization"
  )
})

test_that(".compute_Smin returns positive value", {
  result <- maxtwin:::.compute_Smin(rbar = 0.5, gam = 0.2, p = 2)
  expect_true(result > 0)
  expect_true(is.finite(result))
})

test_that(".compute_Smin increases with rbar", {
  s1 <- maxtwin:::.compute_Smin(rbar = 0.3, gam = 0.2, p = 4)
  s2 <- maxtwin:::.compute_Smin(rbar = 0.6, gam = 0.2, p = 4)
  expect_true(s2 > s1)
})

test_that(".compute_Smin matches Proposition 1 formula", {
  rbar <- 0.4; gam <- 0.2; p <- 6
  expected <- 0.5 * rbar * (1 - gam) / (gam + (1 - gam) * p) *
    beta(1 / p, (1 - gam) / gam)
  result <- maxtwin:::.compute_Smin(rbar, gam, p)
  expect_equal(result, expected)
})
