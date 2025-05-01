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

# ipd 0.99.0

## Summary:

* Preparations to archive on CRAN and move to Bioconductor.

* Slight formatting changes to conform to `styler` and `lintr` suggestions.

* Added PPIa, Chen and Chen methods from Gronsbell et al. (2024) "Another look at inference after prediction."

## Specific Changes: 

- **Version bump**  
  - Pre‑release version set to **0.99.0** for Bioconductor devel.

- **DESCRIPTION updates**  
  - `Depends: R (>= 4.4.0)`  
  - Added `biocViews: Software`  
  - Added `Suggests: BiocStyle, BiocManager`

- **Vignettes**  
  - Converted existing R Markdown vignettes to Bioconductor style with `BiocStyle::html_document` and proper `VignetteIndexEntry` headers.

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
  - `ipd_ppia()` — implements the PPIa estimator for prediction‑powered inference.  
  - `ipd_chenchen()` — implements the Chen & Chen estimator for inference on predicted data.

- **Bioconductor submission prep**  
  - Branch `RELEASE_3_18` created, tag `v0.99.0` applied.  
  - Will request CRAN archiving of the CRAN version upon successful Bioconductor acceptance.

