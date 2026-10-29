# maxtwin

Robust deterministic data splitting for model validation.

`maxtwin` splits data into training and testing sets by balancing two objectives:

1. **Distribution preservation** -- the test set has the same distribution as the full data (measured by energy distance)
2. **Separation** -- test points are far from training points, ensuring non-trivial predictions

The method is robust to paired and clustered data, where methods like random splitting or DUPLEX can fail.

## Installation

```r
# After CRAN release
install.packages("maxtwin")
```

## Usage

```r
library(maxtwin)

set.seed(42)
X <- matrix(runif(500), ncol = 5)
Y <- sin(X[,1]) + cos(X[,2]) + X[,3]^2

# Split with 20% test set
test_idx <- maxtwin(X, Y, splitratio = 0.2)

X_train <- X[-test_idx, ]
X_test  <- X[test_idx, ]
Y_train <- Y[-test_idx]
Y_test  <- Y[test_idx]
```

Factor responses and factor predictor columns are encoded using Helmert
contrasts. Character columns must first be converted to factors.

```r
# Classification with a factor response
test_idx <- maxtwin(
  X = iris[, 1:4],
  Y = iris$Species,
  splitratio = 0.2
)
```

When there is no separate response, splitting can be based on all columns in a
data frame:

```r
test_idx <- maxtwin(iris, splitratio = 0.2)
```

### Step 1 only (without iterative separation improvement)

```r
test_idx <- initialmaxtwin(X, Y, splitratio = 0.2)
```

### Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `splitratio` | -- | Fraction of data for test set. Must be `1/r` for integer `r >= 2` (0.5, 0.25, 0.2, 0.1, ...) |
| `initialization` | `TRUE` | Start twinning from the point closest to the centroid |
| `standardization` | `TRUE` | Standardize the augmented data before twinning |

## Algorithm

MaxTwin (Algorithm 1 in the paper) proceeds in two steps:

**Step 1**: Apply twinning to select candidate test points, then within each Voronoi cell of the twinning set, pick the point with the largest nearest-neighbor distance.

**Step 2**: If the average separation distance falls below the threshold from Proposition 1, iteratively replace test points to improve separation while staying within a ball constraint that preserves distribution quality.


## Dependencies

- [FNN](https://cran.r-project.org/package=FNN) -- KD-tree based nearest neighbor search
- [twinning](https://cran.r-project.org/package=twinning) -- Energy distance minimization

## Testing

```r
devtools::test()
```

Tests cover:
- Output correctness (size, bounds, uniqueness, sorting)
- Determinism (same seed = same result)
- Edge cases (univariate, high-dimensional, paired, and categorical data)
- Separation quality vs random splitting
- Consistency with original implementation

## Reference

Joseph, V.R. and Cho, Y. (in press). A Robust Deterministic Data Splitting Method for Model Validation. *Journal of Computational and Graphical Statistics*.

## License

GPL (>= 3)
