# Data generation function for various underlying models

Data generation function for various underlying models

## Usage

``` r
simdat(
  n = c(300, 300, 300),
  effect = 1,
  sigma_Y = 1,
  model = "ols",
  shift = 0,
  scale = 1
)
```

## Arguments

- n:

  Integer vector of size 3 indicating the sample sizes in the training,
  labeled, and unlabeled data sets, respectively

- effect:

  Regression coefficient for the first variable of interest for
  inference. Defaults is 1.

- sigma_Y:

  Residual variance for the generated outcome. Defaults is 1.

- model:

  The type of model to be generated. Must be one of `"mean"`,
  `"quantile"`, `"ols"`, `"logistic"`, or `"poisson"`. Default is
  `"ols"`.

- shift:

  Scalar shift of the predictions for continuous outcomes (i.e., "mean",
  "quantile", and "ols"). Defaults to 0.

- scale:

  Scaling factor for the predictions for continuous outcomes (i.e.,
  "mean", "quantile", and "ols"). Defaults to 1.

## Value

A data.frame containing n rows and columns corresponding to the labeled
outcome (Y), the predicted outcome (f), a character variable (set_label)
indicating which data set the observation belongs to (training, labeled,
or unlabeled), and four independent, normally distributed predictors
(X1, X2, X3, and X4), where applicable.

## Details

The `simdat` function generates three datasets consisting of independent
realizations of \\Y\\ (for `model` = `"mean"` or `"quantile"`), or
\\\\Y, \boldsymbol{X}\\\\ (for `model` = `"ols"`, `"logistic"`, or
`"poisson"`): a *training* dataset of size \\n_t\\, a *labeled* dataset
of size \\n_l\\, and an *unlabeled* dataset of size \\n_u\\. These sizes
are specified by the argument `n`.

NOTE: In the *unlabeled* data subset, outcome data are still generated
to facilitate a benchmark for comparison with an "oracle" model that
uses the true \\Y^{\mathcal{U}}\\ values for estimation and inference.

**Generating Data**

For `"mean"` and `"quantile"`, we simulate a continuous outcome, \\Y \in
\mathbb{R}\\, with mean given by the `effect` argument and error
variance given by the `sigma_y` argument.

For `"ols"`, `"logistic"`, or `"poisson"` models, predictor data,
\\\boldsymbol{X} \in \mathbb{R}^4\\ are simulated such that the \\i\\th
observation follows a standard multivariate normal distribution with a
zero mean vector and identity covariance matrix:

\$\$ \boldsymbol{X_i} = (X\_{i1}, X\_{i2}, X\_{i3}, X\_{i4}) \sim
\mathcal{N}\_4(\boldsymbol{0}, \boldsymbol{I}). \$\$

For `"ols"`, a continuous outcome \\Y \in \mathbb{R}\\ is simulated to
depend on \\X_1\\ through a linear term with the effect size specified
by the `effect` argument, while the other predictors, \\\boldsymbol{X}
\setminus X_1\\, have nonlinear effects:

\$\$ Y_i = effect \times Z\_{i1} + \frac{1}{2} Z\_{i2}^2 + \frac{1}{3}
Z\_{i3}^3 + \frac{1}{4} Z\_{i4}^2 + \varepsilon_y, \$\$

and \\\varepsilon_y \sim \mathcal{N}(0, sigma_y)\\, where the `sigma_y`
argument specifies the error variance.

For `"logistic"`, we simulate:

\$\$ \Pr(Y_i = 1 \mid \boldsymbol{X}) = logit^{-1}(effect \times
Z\_{i1} + \frac{1}{2} Z\_{i2}^2 + \frac{1}{3} Z\_{i3}^3 + \frac{1}{4}
Z\_{i4}^2 + \varepsilon_y) \$\$

and generate:

\$\$ Y_i \sim Bern\[1, \Pr(Y_i = 1 \mid \boldsymbol{X})\] \$\$

where \\\varepsilon_y \sim \mathcal{N}(0, sigma\\y)\\.

For `"poisson"`, we simulate:

\$\$ \lambda_Y = exp(effect \times Z\_{i1} + \frac{1}{2} Z\_{i2}^2 +
\frac{1}{3} Z\_{i3}^3 + \frac{1}{4} Z\_{i4}^2 + \varepsilon_y) \$\$

and generate:

\$\$ Y_i \sim Poisson(\lambda_Y) \$\$

**Generating Predictions**

To generate predicted outcomes for `"mean"` and `"quantile"`, we
simulate a continuous variable with mean given by the empirical mean of
the training data and error variance given by the `sigma_y` argument.

For `"ols"`, we fit a generalized additive model (GAM) on the simulated
*training* dataset and calculate predictions for the *labeled* and
*unlabeled* datasets as deterministic functions of \\\boldsymbol{X}\\.
Specifically, we fit the following GAM:

\$\$ Y^{\mathcal{T}} = s_0 + s_1(X_1^{\mathcal{T}}) +
s_2(X_2^{\mathcal{T}}) + s_3(X_3^{\mathcal{T}}) +
s_4(X_4^{\mathcal{T}}) + \varepsilon_p, \$\$

where \\\mathcal{T}\\ denotes the *training* dataset, \\s_0\\ is an
intercept term, and \\s_1(\cdot)\\, \\s_2(\cdot)\\, \\s_3(\cdot)\\, and
\\s_4(\cdot)\\ are smoothing spline functions for \\X_1\\, \\X_2\\,
\\X_3\\, and \\X_4\\, respectively, with three target equivalent degrees
of freedom. Residual error is modeled as \\\varepsilon_p\\.

Predictions for *labeled* and *unlabeled* datasets are calculated as:

\$\$ f(\boldsymbol{X}^{\mathcal{L}\cup\mathcal{U}}) = \hat{s}\_0 +
\hat{s}\_1(X_1^{\mathcal{L}\cup\mathcal{U}}) +
\hat{s}\_2(X_2^{\mathcal{L}\cup\mathcal{U}}) +
\hat{s}\_3(X_3^{\mathcal{L}\cup\mathcal{U}}) +
\hat{s}\_4(X_4^{\mathcal{L}\cup\mathcal{U}}), \$\$

where \\\hat{s}\_0, \hat{s}\_1, \hat{s}\_2, \hat{s}\_3\\, and
\\\hat{s}\_4\\ are estimates of \\s_0, s_1, s_2, s_3\\, and \\s_4\\,
respectively.

NOTE: For continuous outcomes, we provide optional arguments `shift` and
`scale` to further apply a location shift and scaling factor,
respectively, to the predicted outcomes. These default to `shift = 0`
and `scale = 1`, i.e., no location shift or scaling.

For `"logistic"`, we train k-nearest neighbors (k-NN) classifiers on the
simulated *training* dataset for values of \\k\\ ranging from 1 to 10.
The optimal \\k\\ is chosen via cross-validation, minimizing the
misclassification error on the validation folds. Predictions for the
*labeled* and *unlabeled* datasets are obtained by applying the k-NN
classifier with the optimal \\k\\ to \\\boldsymbol{X}\\.

Specifically, for each observation in the *labeled* and *unlabeled*
datasets:

\$\$ \hat{Y} = \operatorname{argmax}\_c \sum\_{i \in \mathcal{N}\_k}
I(Y_i = c), \$\$

where \\\mathcal{N}\_k\\ represents the set of \\k\\ nearest neighbors
in the training dataset, \\c\\ indexes the possible classes (0 or 1),
and \\I(\cdot)\\ is an indicator function.

For `"poisson"`, we fit a generalized linear model (GLM) with a log link
function to the simulated *training* dataset. The model is of the form:

\$\$ \log(\mu^{\mathcal{T}}) = \gamma_0 + \gamma_1 X_1^{\mathcal{T}} +
\gamma_2 X_2^{\mathcal{T}} + \gamma_3 X_3^{\mathcal{T}} + \gamma_4
X_4^{\mathcal{T}}, \$\$

where \\\mu^{\mathcal{T}}\\ is the expected count for the response
variable in the *training* dataset, \\\gamma_0\\ is the intercept, and
\\\gamma_1\\, \\\gamma_2\\, \\\gamma_3\\, and \\\gamma_4\\ are the
regression coefficients for the predictors \\X_1\\, \\X_2\\, \\X_3\\,
and \\X_4\\, respectively.

Predictions for the *labeled* and *unlabeled* datasets are calculated
as:

\$\$ \hat{\mu}^{\mathcal{L} \cup \mathcal{U}} = \exp(\hat{\gamma}\_0 +
\hat{\gamma}\_1 X_1^{\mathcal{L} \cup \mathcal{U}} + \hat{\gamma}\_2
X_2^{\mathcal{L} \cup \mathcal{U}} + \hat{\gamma}\_3 X_3^{\mathcal{L}
\cup \mathcal{U}} + \hat{\gamma}\_4 X_4^{\mathcal{L} \cup \mathcal{U}}),
\$\$

where \\\hat{\gamma}\_0\\, \\\hat{\gamma}\_1\\, \\\hat{\gamma}\_2\\,
\\\hat{\gamma}\_3\\, and \\\hat{\gamma}\_4\\ are the estimated
coefficients.

## Examples

``` r
#-- Mean

dat_mean <- simdat(c(100, 100, 100),
  effect = 1, sigma_Y = 1,
  model = "mean"
)

head(dat_mean)
#>            Y  f set_label
#> 1 -0.1074374 NA  training
#> 2  0.7298870 NA  training
#> 3  3.0675835 NA  training
#> 4  1.0686494 NA  training
#> 5  1.6759636 NA  training
#> 6  1.1935525 NA  training

#-- Linear Regression

dat_ols <- simdat(c(100, 100, 100),
  effect = 1, sigma_Y = 1,
  model = "ols"
)

head(dat_ols)
#>          X1          X2          X3         X4         Y  f set_label
#> 1 1.5910561 -0.59940831  0.30908385 -0.2626245 0.5273091 NA  training
#> 2 0.1205511  0.60279555 -0.28111856 -0.0494309 0.3698371 NA  training
#> 3 1.8115029 -3.62378708  1.94845629  0.2513857 8.8433709 NA  training
#> 4 1.3911349 -1.02943895  0.01860701 -0.2364061 0.2271772 NA  training
#> 5 1.5462584 -0.29433239  0.20530142 -1.9968885 2.3138664 NA  training
#> 6 0.6071025  0.09565687 -0.91779026  1.4275819 2.1756725 NA  training
```
