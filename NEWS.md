# ipd 0.1.3

* Added a `NEWS.md` file to track changes to the package.

* Added a `pkgdown` site for the package.

* `ipd()` now allows for regression through the origin with `intercept = FALSE` argument.
  
* `ipd()` now takes an additional argument, `na_action`, to handle missing covariate data.

  * Currently supports `"na.fail"` and `"na.omit"`. Defaults to `na.fail`.
  
  * Provides a more informative error message and lists which covariates are missing observations.
  
* `ipd()` now takes an additional argument, `n_t`, which denotes the (optional) size of the training set used to generate the prediction rule. Defaults to `Inf` but is necessary for the `postpi_X` methods if `n_t` < `n`, `N`, the number of labeled and unlabeled observations, respectively, in the data being analyzed.

# ipd 0.1.4

* Added a help topic for the package itself (`man/ipd-package.Rd`) via `R/ipd-package.R` and `roxygen2`

* Updated the documentation for `ipd()`:

  * Provided a more explicit description of the `model` argument, which is meant to specify the downstream inferential model or parameter to be estimated.
  
  * Clarified that not all columns in data are used in prediction unless explicitly referenced in the `formula` argument or in the `label` argument if the data are passed as one stacked data frame. 

* Updated the documentation for `simdat()` to include a more thorough explanation of how to simulate data with this function. 

* `simdat()` now outputs a `data.frame` with a column named `"set_label"` instead of `"set"` to denote the labeled/unlabeled observation indicator.
