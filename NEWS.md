# ipd 0.1.3

* Added a `NEWS.md` file to track changes to the package.

* Added a `pkgdown` site for the package.

* `ipd()` now allows for regression through the origin with `intercept = FALSE` argument.
  
* `ipd()` now takes an additional argument, `na_action`, to handle missing covariate data.

  * Currently supports `"na.fail"` and `"na.omit"`. Defaults to `na.fail`.
  
  * Provides a more informative error message and lists which covariates are missing observations.
  
* `ipd()` now takes an additional argument, `n_t`, which denotes the (optional) size of the training set used to generate the prediction rule. Defaults to `Inf` but is necessary for the `postpi_X` methods if `n_t` < `n`, `N`, the number of labeled and unlabeled observations, respectively, in the data being analyzed.
