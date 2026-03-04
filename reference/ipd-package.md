# ipd: Inference on Predicted Data

Performs valid statistical inference on predicted data (IPD) using
recent methods, where for a subset of the data, the outcomes have been
predicted by an algorithm. Provides a wrapper function with specified
defaults for the type of model and method to be used for estimation and
inference. Further provides methods for tidying and summarizing results.
Salerno et al., (2025)
[doi:10.1093/bioinformatics/btaf055](https://doi.org/10.1093/bioinformatics/btaf055)
.

The `ipd` package provides tools for statistical modeling and inference
when a significant portion of the outcome data is predicted by AI/ML
algorithms. It implements several state-of-the-art methods for inference
on predicted data (IPD), offering a user-friendly interface to
facilitate their use in real-world applications.

## Details

This package is particularly useful in scenarios where predicted values
(e.g., from machine learning models) are used as proxies for unobserved
outcomes, which can introduce biases in estimation and inference. The
`ipd` package integrates methods designed to address these challenges.

## Features

- Multiple IPD methods: `PostPI`, `PPI`, `PPI++`, and `PSPA` currently.

- Flexible wrapper functions for ease of use.

- Custom methods for model inspection and evaluation.

- Seamless integration with common data structures in R.

- Comprehensive documentation and examples.

## Key Functions

- [`ipd`](https://ipd-tools.github.io/ipd/reference/ipd.md): Main
  wrapper function which implements various methods for inference on
  predicted data for a specified model/outcome type (e.g., mean
  estimation, linear regression).

- [`simdat`](https://ipd-tools.github.io/ipd/reference/simdat.md):
  Simulates data for demonstrating the use of the various IPD methods.

- [`print.ipd`](https://ipd-tools.github.io/ipd/reference/print.ipd.md):
  Prints a brief summary of the IPD method/model combination.

- [`summary.ipd`](https://ipd-tools.github.io/ipd/reference/summary.ipd.md):
  Summarizes the results of fitted IPD models.

- [`tidy.ipd`](https://ipd-tools.github.io/ipd/reference/tidy.ipd.md):
  Tidies the IPD method/model fit into a data frame.

- [`glance.ipd`](https://ipd-tools.github.io/ipd/reference/glance.ipd.md):
  Glances at the IPD method/model fit, returning a one-row summary.

- [`augment.ipd`](https://ipd-tools.github.io/ipd/reference/augment.ipd.md):
  Augments the data used for an IPD method/model fit with additional
  information about each observation.

## Documentation

The package includes detailed documentation for each function, including
usage examples. A vignette is also provided to guide users through
common workflows and applications of the package.

## References

For details on the statistical methods implemented in this package,
please refer to the associated manuscripts at the following references:

- **PostPI**: Wang, S., McCormick, T. H., & Leek, J. T. (2020). Methods
  for correcting inference based on outcomes predicted by machine
  learning. Proceedings of the National Academy of Sciences, 117(48),
  30266-30275.

- **PPI**: Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I.,
  & Zrnic, T. (2023). Prediction-powered inference. Science, 382(6671),
  669-674.

- **PPI++**: Angelopoulos, A. N., Duchi, J. C., & Zrnic, T. (2023).
  PPI++: Efficient prediction-powered inference. arXiv preprint
  arXiv:2311.01453.

- **PSPA**: Miao, J., Miao, X., Wu, Y., Zhao, J., & Lu, Q. (2023).
  Assumption-lean and data-adaptive post-prediction inference. arXiv
  preprint arXiv:2311.14220.

## See also

Useful links:

- <https://github.com/ipd-tools/ipd>

- <https://ipd-tools.github.io/ipd/>

- Report bugs at <https://github.com/ipd-tools/ipd/issues>

## Author

**Maintainer**: Stephen Salerno <ssalerno@fredhutch.org>
([ORCID](https://orcid.org/0000-0003-2763-0494)) \[copyright holder\]

Authors:

- Jiacheng Miao <jmiao24@wisc.edu>

- Awan Afiaz <aafiaz@uw.edu>

- Kentaro Hoffman <khoffm3@uw.edu>

- Jesse Gronsbell <j.gronsbell@utoronto.ca>

- Jianhui Gao <jianhui.gao@mail.utoronto.ca>

- David Cheng <dcheng@mgh.harvard.edu>

- Anna Neufeld <acn2@williams.edu>

- Qiongshi Lu <qlu@biostat.wisc.edu>

- Tyler H McCormick <tylermc@uw.edu>

- Jeffrey T Leek <jtleek@fredhutch.org>

## Examples

``` r
#-- Generate Example Data

set.seed(12345)

dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)

head(dat)
#>           X1          X2         X3          X4         Y  f set_label
#> 1  0.5855288 -0.78486098  1.1872102  1.05076285 1.4008570 NA  training
#> 2  0.7094660 -2.56005244 -0.3567140 -0.07179733 4.1079201 NA  training
#> 3 -0.1093033  0.07280078  1.2122385  0.11673662 1.4501726 NA  training
#> 4 -0.4534972  0.75024358 -0.6939527  0.97786651 1.2987926 NA  training
#> 5  0.6058875 -0.12824888  1.3560616 -1.03154201 2.5256490 NA  training
#> 6 -1.8179560 -0.48786673  0.9057313  2.19912933 0.2889297 NA  training

formula <- Y - f ~ X1

#-- PostPI Analytic Correction (Wang et al., 2020)

fit_postpi1 <- ipd(formula,
  method = "postpi_analytic", model = "ols",
  data = dat, label = "set_label"
)

#-- PostPI Bootstrap Correction (Wang et al., 2020)

nboot <- 200

fit_postpi2 <- ipd(formula,
  method = "postpi_boot", model = "ols",
  data = dat, label = "set_label", nboot = nboot
)

#-- PPI (Angelopoulos et al., 2023)

fit_ppi <- ipd(formula,
  method = "ppi", model = "ols",
  data = dat, label = "set_label"
)

#-- PPI++ (Angelopoulos et al., 2023)

fit_plusplus <- ipd(formula,
  method = "ppi_plusplus", model = "ols",
  data = dat, label = "set_label"
)

#-- PSPA (Miao et al., 2023)

fit_pspa <- ipd(formula,
  method = "pspa", model = "ols",
  data = dat, label = "set_label"
)

#-- Print the Model

print(fit_postpi1)
#> IPD inference summary
#>   Method:   postpi_analytic 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)  0.62451    0.10803  5.7811 7.420e-09 ***
#> X1           0.79449    0.12965  6.1279 8.904e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- Summarize the Model

summ_fit_postpi1 <- summary(fit_postpi1)

#-- Print the Model Summary

print(summ_fit_postpi1)
#> 
#> Call:
#>   Y - f ~ X1 
#> 
#> Method:    postpi_analytic 
#> Model:     ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)  0.62451    0.10803  5.7811 7.420e-09 ***
#> X1           0.79449    0.12965  6.1279 8.904e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- Tidy the Model Output

tidy(fit_postpi1)
#> # A tibble: 2 × 5
#>   term        estimate std.error conf.low conf.high
#>   <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.625     0.108    0.413     0.836
#> 2 X1             0.794     0.130    0.540     1.05 

#-- Get a One-Row Summary of the Model

glance(fit_postpi1)
#> # A tibble: 1 × 6
#>   method          model intercept nobs_labeled nobs_unlabeled call      
#>   <chr>           <chr> <lgl>            <int>          <int> <chr>     
#> 1 postpi_analytic ols   TRUE               300            300 Y - f ~ X1

#-- Augment the Original Data with Fitted Values and Residuals

augmented_df <- augment(fit_postpi1)

head(augmented_df)
#>             X1         X2          X3         X4          Y           f
#> 601 -1.6366291 -0.1746226 -0.18264862 -0.1995761 -1.1215870 -1.27598909
#> 602  0.2115626 -0.6706167 -1.39196798  0.7013954 -0.3900187 -0.27576946
#> 603 -0.4648317  0.5074258  0.70824781  0.8477244  1.2041955  0.80200786
#> 604 -0.6623572  1.2474343  0.18896582  0.2288555  0.1193828 -0.04566093
#> 605 -0.1329536 -1.2482755 -0.21736688 -0.1678426 -0.8563329  0.54360381
#> 606 -1.3217017 -1.9347187  0.07163463 -0.6305805 -0.2815313 -0.40867991
#>     set_label     .fitted      .resid
#> 601 unlabeled -0.67577706 -0.44580990
#> 602 unlabeled  0.79259713 -1.18261583
#> 603 unlabeled  0.25520700  0.94898851
#> 604 unlabeled  0.09827451  0.02110826
#> 605 unlabeled  0.51888164 -1.37521458
#> 606 unlabeled -0.42556966  0.14403841
```
