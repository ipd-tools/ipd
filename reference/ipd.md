# Inference on Predicted Data (ipd)

The main wrapper function to conduct ipd using various methods and
models, and returns a list of fitted model components.

## Usage

``` r
ipd(
  formula,
  method,
  model,
  data,
  label = NULL,
  unlabeled_data = NULL,
  intercept = TRUE,
  alpha = 0.05,
  alternative = "two-sided",
  na_action = "na.fail",
  ...
)
```

## Arguments

- formula:

  An object of class `formula`: a symbolic description of the model to
  be fitted. Must be of the form `Y - f ~ X`, where `Y` is the name of
  the column corresponding to the observed outcome in the labeled data,
  `f` is the name of the column corresponding to the predicted outcome
  in both labeled and unlabeled data, and `X` corresponds to the
  features of interest (i.e., `X = X1 + ... + Xp`). See **1. Formula**
  in the **Details** below for more information.

- method:

  The IPD method to be used for fitting the model. Must be one of
  `"chen"`, `"pdc"`, `"postpi_analytic"`, `"postpi_boot"`, `"ppi"`,
  `"ppi_a"`, `"ppi_plusplus"`, or `"pspa"`. See **3. Method** in the
  **Details** below for more information.

- model:

  The type of downstream inferential model to be fitted, or the
  parameter being estimated. Must be one of `"mean"`, `"quantile"`,
  `"ols"`, `"logistic"`, or `"poisson"`. See **4. Model** in the
  **Details** below for more information.

- data:

  A `data.frame` containing the variables in the model, either a stacked
  data frame with a specific column identifying the labeled versus
  unlabeled observations (`label`), or only the labeled data set. Must
  contain columns for the observed outcomes (`Y`), the predicted
  outcomes (`f`), and the features (`X`) needed to specify the
  `formula`. See **2. Data** in the **Details** below for more
  information.

- label:

  A `string`, `int`, or `logical` specifying the column in the data that
  distinguishes between the labeled and unlabeled observations. See the
  `Details` section for more information. If NULL, `unlabeled_data` must
  be specified. See **2. Data** in the **Details** below for more
  information.

- unlabeled_data:

  (optional) A `data.frame` of unlabeled data. If NULL, `label` must be
  specified. Specifying both the `label` and `unlabeled_data` arguments
  will result in an error message. If specified, must contain columns
  for the predicted outcomes (`f`), and the features (`X`) needed to
  specify the `formula`. See **2. Data** in the **Details** below for
  more information.

- intercept:

  `Logical`. Should an intercept be included in the model? Default is
  `TRUE`.

- alpha:

  The significance level for confidence intervals. Default is `0.05`.

- alternative:

  A string specifying the alternative hypothesis. Must be one of
  `"two-sided"`, `"less"`, or `"greater"`.

- na_action:

  (string, optional) How missing covariate data should be handled.
  Currently `"na.fail"` and `"na.omit"` are accommodated. Defaults to
  `"na.fail"`.

- ...:

  Additional arguments to be passed to the fitting function. See the
  `Details` section for more information. See **5. Auxiliary Arguments**
  and **6. Other Arguments** in the **Details** below for more
  information.

## Value

a summary of model output.

An S4 object of class `IPD` with the following slots:

- `coefficients`:

  Named [`numeric`](https://rdrr.io/r/base/numeric.html) vector of
  estimated parameters.

- `se`:

  Named [`numeric`](https://rdrr.io/r/base/numeric.html) vector of
  standard errors.

- `ci`:

  A [`matrix`](https://rdrr.io/r/base/matrix.html) of confidence
  intervals, with columns `lower` and `upper`.

- `coefTable`:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) summarizing
  Estimate, Std. Error, z-value, and Pr(\>\|z\|) (glm-style).

- `fit`:

  The raw output [`list`](https://rdrr.io/r/base/list.html) returned by
  the method-specific helper function.

- `formula`:

  The [`formula`](https://rdrr.io/r/stats/formula.html) used for fitting
  the IPD model.

- `data_l`:

  The labeled [`data.frame`](https://rdrr.io/r/base/data.frame.html)
  used in the analysis.

- `data_u`:

  The unlabeled [`data.frame`](https://rdrr.io/r/base/data.frame.html)
  used in the analysis.

- `method`:

  A [`character`](https://rdrr.io/r/base/character.html) string
  indicating which IPD method was applied.

- `model`:

  A [`character`](https://rdrr.io/r/base/character.html) string
  indicating the downstream inferential model.

- `intercept`:

  A [`logical`](https://rdrr.io/r/base/logical.html) indicating whether
  an intercept was included.

## Details

**1. Formula:**

The `ipd` function uses one formula argument that specifies both the
calibrating model (e.g., PostPI "relationship model", PPI "rectifier"
model) and the inferential model. These separate models will be created
internally based on the specific `method` called.

**2. Data:**

The data can be specified in two ways:

1.  Single data argument (`data`) containing a stacked `data.frame` and
    a label identifier (`label`).

2.  Two data arguments, one for the labeled data (`data`) and one for
    the unlabeled data (`unlabeled_data`).

For option (1), provide one data argument (`data`) which contains a
stacked `data.frame` with both the unlabeled and labeled data and a
`label` argument that specifies the column identifying the labeled
versus the unlabeled observations in the stacked `data.frame` (e.g.,
`label = "set_label"` if the column "set_label" in the stacked data
denotes which set an observation belongs to).

NOTE: Labeled data identifiers can be:

- String:

  "l", "lab", "label", "labeled", "labelled", "tst", "test", "true"

- Logical:

  TRUE

- Factor:

  Non-reference category (i.e., binary 1)

Unlabeled data identifiers can be:

- String:

  "u", "unlab", "unlabeled", "unlabelled", "val", "validation", "false"

- Logical:

  FALSE

- Factor:

  Non-reference category (i.e., binary 0)

For option (2), provide separate data arguments for the labeled data set
(`data`) and the unlabeled data set (`unlabeled_data`). If the second
argument is provided, the function ignores the `label` identifier and
assumes the data provided are not stacked.

NOTE: Not all columns in `data` or `unlabeled_data` may be used unless
explicitly referenced in the `formula` argument or in the `label`
argument (if the data are passed as one stacked data frame).

**3. Method:**

Use the `method` argument to specify the fitting method:

- "chen":

  Gronsbell et al. (2026) Chen and Chen Correction

- "pdc":

  Gan et al. (2024) Prediction Decorrelated Inference

- "postpi_analytic":

  Wang et al. (2020) Post-Prediction Inference (PostPI) Analytic
  Correction

- "postpi_boot":

  Wang et al. (2020) Post-Prediction Inference (PostPI) Bootstrap
  Correction

- "ppi":

  Angelopoulos et al. (2023) Prediction-Powered Inference (PPI)

- "ppi_a":

  Gronsbell et al. (2025) PPI "All" Correction

- "ppi_plusplus":

  Angelopoulos et al. (2023) PPI++

- "pspa":

  Miao et al. (2023) Assumption-Lean and Data-Adaptive Post-Prediction
  Inference (PSPA)

**4. Model:**

Use the `model` argument to specify the type of downstream inferential
model or parameter to be estimated:

- "mean":

  Mean value of a continuous outcome

- "quantile":

  `q`th quantile of a continuous outcome

- "ols":

  Linear regression coefficients for a continuous outcome

- "logistic":

  Logistic regression coefficients for a binary outcome

- "poisson":

  Poisson regression coefficients for a count outcome

The `ipd` wrapper function will concatenate the `method` and `model`
arguments to identify the required helper function, following the naming
convention "method_model".

**5. Auxiliary Arguments:**

The wrapper function will take method-specific auxiliary arguments
(e.g., `q` for the quantile estimation models) and pass them to the
helper function through the "..." with specified defaults for
simplicity.

**6. Other Arguments:**

All other arguments that relate to all methods (e.g., alpha, ci.type),
or other method-specific arguments, will have defaults.

## Examples

``` r
#-- Generate Example Data

dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)

head(dat)
#>            X1          X2         X3         X4          Y  f set_label
#> 1  0.91295168  0.82805181  0.9800060 -1.3185450  1.9821781 NA  training
#> 2  0.34629713  0.95618711  1.5572516 -0.4150736  1.9104987 NA  training
#> 3 -0.01794417 -0.02357563 -0.6602489  1.7181628 -0.5284219 NA  training
#> 4  0.42511788 -0.20367165 -0.4274529  0.3139744  0.6832558 NA  training
#> 5  0.98476232 -0.66591755 -1.4923011 -0.2940388  0.3307442 NA  training
#> 6 -0.57141674 -0.34011487  0.8388489  0.7944130 -2.2485908 NA  training

formula <- Y - f ~ X1

#-- Chen and Chen Correction (Gronsbell et al., 2026)

ipd(formula,
  method = "chen", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   chen 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.653932   0.084634  7.7266 1.105e-14 ***
#> X1          1.041767   0.098024 10.6277 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- Prediction Decorrelated Inference (Gan et al., 2024)

ipd(formula,
  method = "chen", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   chen 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.653932   0.084634  7.7266 1.105e-14 ***
#> X1          1.041767   0.098024 10.6277 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- PostPI Analytic Correction (Wang et al., 2020)

ipd(formula,
  method = "postpi_analytic", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   postpi_analytic 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)  0.69068    0.11352  6.0843 1.170e-09 ***
#> X1           1.03082    0.13585  7.5877 3.255e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- PostPI Bootstrap Correction (Wang et al., 2020)

nboot <- 200

ipd(formula,
  method = "postpi_boot", model = "ols",
  data = dat, label = "set_label", nboot = nboot
)
#> IPD inference summary
#>   Method:   postpi_boot 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.635274   0.095072   6.682 2.357e-11 ***
#> X1          1.024740   0.094036  10.897 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- PPI (Angelopoulos et al., 2023)

ipd(formula,
  method = "ppi", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   ppi 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.636905   0.093397  6.8194 9.145e-12 ***
#> X1          1.032763   0.100322 10.2945 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- PPI "All" (Gronsbell et al., 2025)

ipd(formula,
  method = "ppi_a", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   ppi_a 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.655324   0.084638  7.7426 9.737e-15 ***
#> X1          1.046635   0.093822 11.1555 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- PPI++ (Angelopoulos et al., 2023)

ipd(formula,
  method = "ppi_plusplus", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   ppi_plusplus 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.652305   0.085908  7.5931 3.123e-14 ***
#> X1          1.044527   0.096407 10.8346 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- PSPA (Miao et al., 2023)

ipd(formula,
  method = "pspa", model = "ols",
  data = dat, label = "set_label"
)
#> IPD inference summary
#>   Method:   pspa 
#>   Model:    ols 
#>   Formula:  Y - f ~ X1 
#> 
#> Coefficients:
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.656443   0.085606  7.6682 1.744e-14 ***
#> X1          1.039811   0.096154 10.8141 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
