# R CMD check results

── R CMD check results ────────────────────────────────────────────────────────────────────────── ipd 0.3.0 ────
Duration: 2m 29s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

# Downstream dependencies

character(0)

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
