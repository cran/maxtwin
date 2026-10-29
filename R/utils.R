#' @importFrom FNN knn.dist knnx.dist knnx.index
#' @importFrom twinning twin
NULL

#' Validate inputs for maxtwin functions
#' @noRd
.validate_inputs <- function(X, Y, splitratio) {
  if (!is.matrix(X) && !is.data.frame(X))
    stop("X must be a matrix or data frame.", call. = FALSE)

  X_data <- as.data.frame(X, stringsAsFactors = FALSE)
  if (ncol(X_data) < 1L)
    stop("X must contain at least one column.", call. = FALSE)

  valid_X <- vapply(
    X_data,
    function(x) is.numeric(x) || is.factor(x),
    logical(1)
  )
  if (any(!valid_X))
    stop(
      "X columns must be numeric or factors. Convert nominal character ",
      "columns with factor().",
      call. = FALSE
    )

  if (anyNA(X_data))
    stop("X contains missing values.", call. = FALSE)

  finite_X <- vapply(
    X_data,
    function(x) !is.numeric(x) || all(is.finite(x)),
    logical(1)
  )
  if (any(!finite_X))
    stop("X contains non-finite numeric values.", call. = FALSE)

  N <- nrow(X_data)

  if (!is.null(Y)) {
    if (length(Y) != N)
      stop("length(Y) must equal nrow(X).", call. = FALSE)

    if (is.matrix(Y) || is.data.frame(Y) ||
        (!is.numeric(Y) && !is.factor(Y)))
      stop("Y must be NULL or a numeric or factor vector.", call. = FALSE)

    if (anyNA(Y))
      stop("Y contains missing values.", call. = FALSE)

    if (is.numeric(Y) && any(!is.finite(Y)))
      stop("Y contains non-finite numeric values.", call. = FALSE)
  }

  if (missing(splitratio) || !is.numeric(splitratio) ||
      length(splitratio) != 1L || !is.finite(splitratio) || splitratio <= 0)
    stop("splitratio must be a single finite positive number.", call. = FALSE)

  # twinning::twin requires r = 1/splitratio to be integer >= 2
  r <- 1 / splitratio
  if (abs(r - round(r)) > 1e-8 || round(r) < 2L)
    stop(
      "splitratio must be the reciprocal of an integer r >= 2 because ",
      "twinning::twin() expects integer r = 1/splitratio. ",
      "Received splitratio = ", splitratio, ", so 1/splitratio = ", r, ". ",
      "Examples: 0.5, 1/3, 0.25, 0.2, 0.1.",
      call. = FALSE
    )

  if (N < 2 * round(r))
    stop("N must be at least 2 * (1/splitratio) = ", 2 * round(r), ".",
         call. = FALSE)

  constant_X <- vapply(
    X_data,
    function(x) length(unique(x)) < 2L,
    logical(1)
  )
  if (any(constant_X))
    stop(
      "X contains constant column(s): ",
      paste(names(X_data)[constant_X], collapse = ", "),
      ". Remove or perturb them.",
      call. = FALSE
    )

  if (!is.null(Y) && length(unique(Y)) < 2L)
    stop("Y is constant. Provide meaningful response values.", call. = FALSE)

  invisible(NULL)
}


#' Convert numeric and factor columns to a numeric matrix
#' @noRd
.encode_helmert <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  blocks <- vector("list", ncol(data))

  for (j in seq_along(data)) {
    x <- data[[j]]
    column_name <- names(data)[j]

    if (is.numeric(x)) {
      block <- matrix(as.double(x), ncol = 1L)
      colnames(block) <- column_name
    } else {
      observed_levels <- unique(x)
      contrasts <- stats::contr.helmert(length(observed_levels))
      block <- contrasts[
        match(x, observed_levels),
        ,
        drop = FALSE
      ]
      storage.mode(block) <- "double"
      colnames(block) <- paste0(
        column_name,
        ".H",
        seq_len(ncol(block))
      )
    }

    blocks[[j]] <- block
  }

  do.call(cbind, blocks)
}


#' Validate and encode predictors and response
#' @noRd
.prepare_data <- function(X, Y, splitratio) {
  .validate_inputs(X, Y, splitratio)

  X_numeric <- .encode_helmert(X)
  Z_numeric <- X_numeric
  if (!is.null(Y)) {
    Y_numeric <- .encode_helmert(
      data.frame(Y = Y, check.names = FALSE)
    )
    Z_numeric <- cbind(X_numeric, Y_numeric)
  }

  list(
    X = X_numeric,
    Z = Z_numeric
  )
}


#' Select one distinct point from each Voronoi cell
#' @noRd
.select_voronoi_points <- function(centers, data, cluster, nndis) {
  n_centers <- nrow(centers)
  selected <- integer(n_centers)
  empty <- logical(n_centers)

  for (i in seq_len(n_centers)) {
    cell <- which(cluster == i)
    if (length(cell) == 0L) {
      empty[i] <- TRUE
    } else {
      selected[i] <- cell[which.max(nndis[cell])]
    }
  }

  used <- selected[!empty]
  for (i in which(empty)) {
    available <- setdiff(seq_len(nrow(data)), used)
    squared_distance <- rowSums(
      sweep(data[available, , drop = FALSE], 2, centers[i, ], "-")^2
    )
    nearest <- available[squared_distance == min(squared_distance)]
    selected[i] <- nearest[which.max(nndis[nearest])]
    used <- c(used, selected[i])
  }

  selected
}


#' Run twinning algorithm with initialization options
#' @noRd
.run_twinning <- function(Z, splitratio, initialization, standardization) {
  if (!is.logical(initialization) || length(initialization) != 1L ||
      is.na(initialization))
    stop("initialization must be TRUE or FALSE.", call. = FALSE)

  if (!is.logical(standardization) || length(standardization) != 1L ||
      is.na(standardization))
    stop("standardization must be TRUE or FALSE.", call. = FALSE)

  r <- round(1 / splitratio)
  Zin <- if (standardization) scale(Z) else Z
  if (initialization) {
    zbar <- colMeans(Z)
    u1 <- which.min(rowSums(sweep(Z, 2, zbar, "-")^2))
    twinning::twin(Zin, r = r, u1 = u1, format_data = FALSE)
  } else {
    twinning::twin(Zin, r = r, format_data = FALSE)
  }
}


#' Compute S_min threshold from Proposition 1 -- eq (8)
#' @noRd
.compute_Smin <- function(rbar, gam, p) {
  0.5 * rbar * (1 - gam) / (gam + (1 - gam) * p) *
    beta(1 / p, (1 - gam) / gam)
}
