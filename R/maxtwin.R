#' @title Initial MaxTwin (Step 1)
#'
#' @description Selects testing points by maximizing nearest-neighbor distance
#'   within each Voronoi cell induced by the twinning set of Vakayil and Joseph
#'   (2022). This is Step 1 of the algorithm of Joseph and Cho (in press), without
#'   iterative improvement of the separation distance.
#'
#' @param X Matrix or data frame of inputs (N x p). Columns may be numeric or
#'   factors.
#' @param Y Optional numeric or factor vector of outputs (length N). If
#'   supplied, it is included in twinning. If omitted, splitting is based on
#'   \code{X} alone.
#' @param splitratio Fraction of data for the testing set. Must be \code{1/r}
#'   where \code{r} is an integer >= 2 (e.g., 0.5, 0.25, 0.2, 0.1).
#'   This constraint comes from the twinning algorithm.
#' @param initialization Logical. If \code{TRUE}, start twinning from the point
#'   closest to the data centroid.
#' @param standardization Logical. If \code{TRUE}, standardize the encoded
#'   data used for twinning before selection.
#' @return A vector of testing point indices.
#'
#' @references
#' Vakayil, A. and Joseph, V.R. (2022). Data Twinning. \emph{Statistical
#' Analysis and Data Mining: The ASA Data Science Journal}, 15(5), 598--610.
#' \doi{10.1002/sam.11574}.
#'
#' Joseph, V.R. and Cho, Y. (in press). A Robust Deterministic Data Splitting
#' Method for Model Validation. \emph{Journal of Computational and Graphical
#' Statistics}.
#'
#' @export
initialmaxtwin <- function(X, Y = NULL, splitratio,
                           initialization = TRUE,
                           standardization = TRUE) {
  prepared <- .prepare_data(X, Y, splitratio)
  X <- prepared$X
  Z <- prepared$Z

  N <- nrow(X)
  p <- ncol(X)

  twin.ind <- .run_twinning(Z, splitratio, initialization, standardization)
  Ntest <- length(twin.ind)

  s <- X[twin.ind, , drop = FALSE]
  cluster <- as.integer(FNN::knnx.index(s, X, k = 1))
  nndis <- as.numeric(FNN::knn.dist(X, k = 1))
  test_idx <- .select_voronoi_points(s, X, cluster, nndis)

  sort(test_idx)
}


#' @title MaxTwin Algorithm
#'
#' @description Implements the robust deterministic data-splitting algorithm
#'   of Joseph and Cho (in press). Step 1 selects testing points using twinning,
#'   Voronoi partitioning, and maximum nearest-neighbor distance. Step 2
#'   iteratively improves separation if needed.
#'
#' @inheritParams initialmaxtwin
#' @return A vector of testing point indices.
#'
#' @details
#' See Joseph and Cho (in press) for details.
#'
#' Note that \code{splitratio} must be \code{1/r} for integer \code{r >= 2},
#' due to the constraint of the underlying twinning algorithm.
#'
#' @references
#' Joseph, V.R. and Cho, Y. (in press). A Robust Deterministic Data Splitting
#' Method for Model Validation. \emph{Journal of Computational and Graphical
#' Statistics}.
#'
#' @examples
#' set.seed(42)
#' X <- matrix(runif(200), ncol = 2)
#' Y <- sin(X[,1]) + cos(X[,2])
#' idx <- maxtwin(X, Y, splitratio = 0.2)
#' X_test  <- X[idx, ]
#' X_train <- X[-idx, ]
#'
#' @export
maxtwin <- function(X, Y = NULL, splitratio,
                    initialization = TRUE,
                    standardization = TRUE) {
  prepared <- .prepare_data(X, Y, splitratio)
  X <- prepared$X
  Z <- prepared$Z

  N <- nrow(X)
  p <- ncol(X)
  gam <- splitratio

  idx <- seq_len(N)

  # ---- Step 1: Twinning ----
  twin.ind <- .run_twinning(Z, splitratio, initialization, standardization)
  Ntest <- length(twin.ind)

  ind.sel <- twin.ind
  s <- X[ind.sel, , drop = FALSE]

  # Rearrange: test candidates first, training after
  X   <- rbind(s, X[-ind.sel, , drop = FALSE])
  idx <- c(idx[ind.sel], idx[-ind.sel])
  # Voronoi assignment + select max-d in each cell
  cluster <- as.integer(FNN::knnx.index(s, X, k = 1))
  d <- as.numeric(FNN::knn.dist(X, k = 1))
  ind.sel <- .select_voronoi_points(s, X, cluster, d)

  # Compute separation distance
  dis <- as.numeric(FNN::knnx.dist(X[-ind.sel, , drop = FALSE],
                                    X[ind.sel, , drop = FALSE], k = 1))
  dbar <- mean(dis)

  # Reorder by separation (worst first)
  ord <- order(dis)
  ind.sel <- ind.sel[ord]
  s <- X[ind.sel, , drop = FALSE]

  # Rearrange with updated selection
  X   <- rbind(s, X[-ind.sel, , drop = FALSE])
  idx <- c(idx[ind.sel], idx[-ind.sel])

  # ---- Step 2: Ensuring Separation (Algorithm 1, lines for ell=2..Ntest) ----
  inds <- seq_len(Ntest)
  cluster <- as.integer(FNN::knnx.index(s, X, k = 1))
  d <- as.numeric(FNN::knn.dist(X, k = 1))

  K <- max(min(ceiling(2 * p), floor(Ntest / 2)), 1L)
  r <- FNN::knn.dist(s, k = K)[, K]
  rbar <- mean(r)

  # Pre-compute constant part of S_min (does not depend on rbar)
  Smin_const <- 0.5 * (1 - gam) / (gam + (1 - gam) * p) *
    beta(1 / p, (1 - gam) / gam)

  if (dbar < Smin_const * rbar) {
    # Step 2 loop: ell = 2, ..., Ntest per Algorithm 1
    for (i in 2:Ntest) {
      # Voronoi for remaining active centers (i:Ntest), excluding t_1..t_{i-1}
      cluster <- as.integer(FNN::knnx.index(
        matrix(s[-(1:(i - 1)), ], ncol = p), X, k = 1
      )) + (i - 1L)
      d <- c(rep(0, i - 1L), as.numeric(
        FNN::knn.dist(X[-(1:(i - 1)), , drop = FALSE], k = 1)
      ))

      # A replacement cannot reuse a point already selected elsewhere.
      occupied <- setdiff(seq_len(Ntest), i)
      ind <- setdiff(which(cluster == i), occupied)
      cell_size <- length(ind)

      if (cell_size > 1L) {
        # Ball constraint: penalize points outside r_i radius
        rad <- sqrt(rowSums(sweep(
          matrix(X[ind, ], ncol = p), 2, s[i, ], "-"
        )^2))
        d[ind[rad > r[i]]] <- -1

        inds[i] <- ind[which.max(d[ind])]
      } else {
        inds[i] <- i
      }

      # Rearrange with updated selection
      s   <- X[inds, , drop = FALSE]
      X   <- rbind(s, X[-inds, , drop = FALSE])
      idx <- c(idx[inds], idx[-inds])
      inds <- seq_len(Ntest)

      # Recheck separation threshold
      dbar <- mean(as.numeric(
        FNN::knnx.dist(X[-inds, , drop = FALSE], s, k = 1)
      ))
      r <- FNN::knn.dist(s, k = K)[, K]
      rbar <- mean(r)

      if (dbar > Smin_const * rbar) break
    }
  }

  sort(idx[seq_len(Ntest)])
}
