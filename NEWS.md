# ipd 0.4.0

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

# ipd 0.2.0

## Summary:

* Preparations to archive on CRAN and move to Bioconductor.

* Slight formatting changes to conform to `styler` and `lintr` suggestions.

* Added PPIa, Chen and Chen methods from Gronsbell et al. (2025) "Another look at inference after prediction."

## Specific Changes: 

- **Version bump**  
  - Pre‑release version set to **0.99.0** for Bioconductor devel.

- **DESCRIPTION updates**  
  - `Depends: R (>= 4.4.0)`  
  - Added `biocViews: Software`  
  - Added `Suggests: BiocStyle, BiocManager`

- **Vignettes**  
  - Converted existing R Markdown vignettes to Bioconductor style with `BiocStyle::html_document` and proper `VignetteIndexEntry` headers.
  - Added examples for Chen and Chen, PPI "All" methods.

- **NEWS & CITATION**  
  - Added `NEWS.md` entry (this file) and a `CITATION` file for package citation metadata.

- **Testing & QA**  
  - Passed `BiocCheck` with no errors or warnings.  
  - Updated `testthat` suite as needed for Bioc compliance.

- **Continuous Integration**  
  - Added GitHub Actions via `usethis::use_github_action("bioc-workflow")` to run Bioconductor checks on Linux, macOS, and Windows.

- **README**  
  - Installation instructions updated to:
    ```r
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install("ipd")
    ```
  - Replaced CRAN build badge with:
    ```markdown
    [![Bioc build status](https://bioconductor.org/shields/build/release/bioc/ipd.svg)](https://bioconductor.org/packages/ipd)
    ```
    
- **New functions**  
  - `ppi_a_ols()`: implements the PPIa estimator for prediction‑powered inference.  
  - `chen_ols()`: implements the Chen & Chen estimator for inference on predicted data.
  - `.parse_inputs()` - helper to validate and split input data.
  - `.drop_unused_levels()` - helper to drop unused factor levels and report which were removed.
  - `.warn_differing_levels()` - helper to warn on differing factor levels between labeled and unlabeled data.
  - `.build_design()` - helper to build design matrices and outcome vectors.
  - `show()` - implements S4 method for `ipd` class.
  
- **Updates to functions**
  - `ipd()` - registered new methods in wrapper.
  - `ipd()` - helper functions for parsing inputs and additional warnings to users when parsing formulas.
  - `methods.R` - Updated with new S4 class for `ipd`.

- **Bioconductor submission prep**  
  - Branch created, tag `v0.99.0` applied.  
  - Will request CRAN archiving of the CRAN version upon successful Bioconductor acceptance.

# ipd 0.3.0

## Summary:

* Refactored internal Chen-Chen and PDC functions to centralize matrix-based estimators and reduce duplicated logic across models.

* Added thin ipd-compatible wrappers for OLS, logistic, and Poisson Chen-Chen and PDC methods.

* Improved internal GLM handling and documentation to support ongoing methodological development.

## Specific Changes:

* **Chen / PDC Refactor**

  * Added new internal utilities file `utils_chen.R` housing shared score/Hessian computation and matrix-based estimators.
  * Centralized core augmented estimators:

    * `compute_min_mse_est()`: Chen-Chen min-MSE form
    * `compute_pdc_est()`: PDC augmentation
    * `compute_min_mse_est_new()`: experimental Chen-Chen variant with separate auxiliary design/family
  * Added shared internal helpers:

    * `compute_residuals()`: GLM residual and Hessian computation
    * `.fit_glm_from_matrix()`: GLM fitting directly from design matrices
    * `.resolve_family()`: robust family resolution for Gaussian, Binomial, and Poisson models
    * `expit()` and `expit_derivative()`

* **Thin Wrapper Functions (ipd helpers)**

  * Rewrote `chen_ols()` to use centralized matrix-based estimators.
  * Added thin wrappers preserving existing `ipd()` method lookup:

    * `chen_logistic()`
    * `chen_poisson()`
    * `pdc_ols()`
    * `pdc_logistic()`
    * `pdc_poisson()`
  * All wrappers maintain backward compatibility with the existing helper contract (`list(est, se)`).

* **Internal Documentation**

  * Added Roxygen2 documentation to all Chen/PDC utilities (including non-exported helpers).
  * Introduced consistent section headers for estimator components and wrappers.
  * Adding PDC method examples to documentation and vignette.
  * Updating urls throughout to use doi.org links instead of journal-specific links.

* **Development Notes**

  * `compute_min_mse_est_new()` is currently internal-only and not yet wired into the `ipd()` method registry; future work will determine exposure via a new method (e.g., `chen_new`) or argument toggle.
