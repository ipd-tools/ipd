## Summary of new changes

* Added a help topic for the package itself (`man/ipd-package.Rd`) via `R/ipd-package.R` and `roxygen2`

* Updated the documentation for `ipd()`:

  * Provided a more explicit description of the `model` argument, which is meant to specify the downstream inferential model or parameter to be estimated.
  
  * Clarified that not all columns in data are used in prediction unless explicitly referenced in the `formula` argument or in the `label` argument if the data are passed as one stacked data frame. 

* Updated the documentation for `simdat()` to include a more thorough explanation of how to simulate data with this function. 

* `simdat()` now outputs a `data.frame` with a column named `"set_label"` instead of `"set"` to denote the labeled/unlabeled observation indicator.


## R CMD check results

── R CMD check results ────────────────────────────────────────────────────────────────────────────────── ipd 0.1.4 ────
Duration: 3m 51.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

character(0)
