#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' ipd: Inference on Predicted Data
#'
#' The `ipd` package provides tools for statistical modeling and inference when
#' a significant portion of the outcome data is predicted by AI/ML algorithms.
#' It implements several state-of-the-art methods for inference on predicted
#' data (IPD), offering a user-friendly interface to facilitate their use in
#' real-world applications.
#'
#' This package is particularly useful in scenarios where predicted values
#' (e.g., from machine learning models) are used as proxies for unobserved
#' outcomes, which can introduce biases in estimation and inference. The `ipd`
#' package integrates methods designed to address these challenges.
#'
#' @section Features:
#' - Multiple IPD methods: `PostPI`, `PPI`, `PPI++`, and `PSPA` currently.
#' - Flexible wrapper functions for ease of use.
#' - Custom methods for model inspection and evaluation.
#' - Seamless integration with common data structures in R.
#' - Comprehensive documentation and examples.
#'
#' @section Key Functions:
#' - \code{\link{ipd}}: Main wrapper function which implements various methods for inference on predicted data for a specified model/outcome type (e.g., mean estimation, linear regression).
#' - \code{\link{simdat}}: Simulates data for demonstrating the use of the various IPD methods.
#' - \code{\link{print.ipd}}: Prints a brief summary of the IPD method/model combination.
#' - \code{\link{summary.ipd}}: Summarizes the results of fitted IPD models.
#' - \code{\link{tidy.ipd}}: Tidies the IPD method/model fit into a data frame.
#' - \code{\link{glance.ipd}}: Glances at the IPD method/model fit, returning a one-row summary.
#' - \code{\link{augment.ipd}}: Augments the data used for an IPD method/model fit with additional information about each observation.
#'
#' @section Documentation:
#' The package includes detailed documentation for each function, including
#' usage examples. A vignette is also provided to guide users through common
#' workflows and applications of the package.
#'
#' @section References:
#' For details on the statistical methods implemented in this package, please
#' refer to the associated manuscripts at the following references:
#' - \strong{PostPI}: Wang, S., McCormick, T. H., & Leek, J. T. (2020). Methods for correcting inference based on outcomes predicted by machine learning. Proceedings of the National Academy of Sciences, 117(48), 30266-30275.
#' - \strong{PPI}: Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic, T. (2023). Prediction-powered inference. Science, 382(6671), 669-674.
#' - \strong{PPI++}: Angelopoulos, A. N., Duchi, J. C., & Zrnic, T. (2023). PPI++: Efficient prediction-powered inference. arXiv preprint arXiv:2311.01453.
#' - \strong{PSPA}: Miao, J., Miao, X., Wu, Y., Zhao, J., & Lu, Q. (2023). Assumption-lean and data-adaptive post-prediction inference. arXiv preprint arXiv:2311.14220.
#'
#' @name ipd-package
#'
#' @keywords package
#'
#' @examples
#' #-- Generate Example Data
#'
#' set.seed(12345)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- PostPI Analytic Correction (Wang et al., 2020)
#'
#' fit_postpi1 <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' #-- PostPI Bootstrap Correction (Wang et al., 2020)
#'
#' nboot <- 200
#'
#' fit_postpi2 <- ipd(formula, method = "postpi_boot", model = "ols",
#'
#'     data = dat, label = "set_label", nboot = nboot)
#'
#' #-- PPI (Angelopoulos et al., 2023)
#'
#' fit_ppi <- ipd(formula, method = "ppi", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' #-- PPI++ (Angelopoulos et al., 2023)
#'
#' fit_plusplus <- ipd(formula, method = "ppi_plusplus", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' #-- PSPA (Miao et al., 2023)
#'
#' fit_pspa <- ipd(formula, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' #-- Print the Model
#'
#' print(fit_postpi1)
#'
#' #-- Summarize the Model
#'
#' summ_fit_postpi1 <- summary(fit_postpi1)
#'
#' #-- Print the Model Summary
#'
#' print(summ_fit_postpi1)
#'
#' #-- Tidy the Model Output
#'
#' tidy(fit_postpi1)
#'
#' #-- Get a One-Row Summary of the Model
#'
#' glance(fit_postpi1)
#'
#' #-- Augment the Original Data with Fitted Values and Residuals
#'
#' augmented_df <- augment(fit_postpi1)
#'
#' head(augmented_df)
## usethis namespace: end
NULL
