# maxtwin 0.1.0

* Initial release.
* Added `maxtwin()` for deterministic train/test splitting with distribution
  preservation and train/test separation.
* Added `initialmaxtwin()` for the Step 1 split without the iterative separation
  improvement.
* Added input validation, examples, and test coverage for deterministic output,
  split sizes, paired data, and high-dimensional inputs.
* Added Helmert coding for factor predictors and factor responses.
* Added predictor-only splitting when `Y` is omitted.
