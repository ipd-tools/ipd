#===============================================================================
#  DATA GENERATION FUNCTION
#===============================================================================

#--- DATA GENERATION FOR VARIOUS MODELS ----------------------------------------

#' Data generation function for various underlying models
#'
#' @param n Integer vector of size 3 indicating the sample sizes in the
#' training, labeled, and unlabeled data sets, respectively
#'
#' @param effect Regression coefficient for the first variable of interest for
#' inference. Defaults is 1.
#'
#' @param sigma_Y Residual variance for the generated outcome. Defaults is 1.
#'
#' @param model The type of model to be generated. Must be one of
#' \code{"mean"}, \code{"quantile"}, \code{"ols"}, \code{"logistic"}, or
#' \code{"poisson"}. Default is \code{"ols"}.
#'
#' @param shift Scalar shift of the predictions for continuous outcomes
#' (i.e., "mean", "quantile", and "ols"). Defaults to 0.
#'
#' @param scale Scaling factor for the predictions for continuous outcomes
#' (i.e., "mean", "quantile", and "ols"). Defaults to 1.
#'
#  @inheritParams gam::gam
#'
#' @return A data.frame containing n rows and columns corresponding to
#' the labeled outcome (Y), the predicted outcome (f), a character variable
#' (set_label) indicating which data set the observation belongs to (training,
#' labeled, or unlabeled), and four independent, normally distributed
#' predictors (X1, X2, X3, and X4), where applicable.
#'
#' @details
#'
#' The `simdat` function generates three datasets consisting of independent
#' realizations of \eqn{Y} (for \code{model} = \code{"mean"} or
#' \code{"quantile"}), or \eqn{\{Y, \boldsymbol{X}\}} (for \code{model} =
#' \code{"ols"}, \code{"logistic"}, or \code{"poisson"}): a \emph{training}
#' dataset of size \eqn{n_t}, a \emph{labeled} dataset of size \eqn{n_l}, and
#' an \emph{unlabeled} dataset of size \eqn{n_u}. These sizes are specified by
#' the argument \code{n}.
#'
#' NOTE: In the \emph{unlabeled} data subset, outcome data are still generated
#' to facilitate a benchmark for comparison with an "oracle" model that uses
#' the true \eqn{Y^{\mathcal{U}}} values for estimation and inference.
#'
#' \strong{Generating Data}
#'
#' For \code{"mean"} and \code{"quantile"}, we simulate a continuous outcome,
#' \eqn{Y \in \mathbb{R}}, with mean given by the \code{effect} argument and
#' error variance given by the \code{sigma_y} argument.
#'
#' For \code{"ols"}, \code{"logistic"}, or \code{"poisson"} models, predictor
#' data, \eqn{\boldsymbol{X} \in \mathbb{R}^4} are simulated such that the
#' \eqn{i}th observation follows a standard multivariate normal distribution
#' with a zero mean vector and identity covariance matrix:
#'
#' \deqn{
#'   \boldsymbol{X_i} = (X_{i1}, X_{i2}, X_{i3}, X_{i4}) \sim \mathcal{N}_4(\boldsymbol{0}, \boldsymbol{I}).
#' }
#'
#' For \code{"ols"}, a continuous outcome \eqn{Y \in \mathbb{R}} is simulated
#' to depend on \eqn{X_1} through a linear term with the effect size specified
#' by the \code{effect} argument, while the other predictors,
#' \eqn{\boldsymbol{X} \setminus X_1}, have nonlinear effects:
#'
#' \deqn{
#'   Y_i = effect \times Z_{i1} + \frac{1}{2} Z_{i2}^2 + \frac{1}{3} Z_{i3}^3 + \frac{1}{4} Z_{i4}^2 + \varepsilon_y,
#' }
#'
#' and \eqn{\varepsilon_y \sim \mathcal{N}(0, sigma_y)}, where the
#' \code{sigma_y} argument specifies the error variance.
#'
#' For \code{"logistic"}, we simulate:
#'
#' \deqn{
#'   \Pr(Y_i = 1 \mid \boldsymbol{X}) = logit^{-1}(effect \times Z_{i1} + \frac{1}{2} Z_{i2}^2 + \frac{1}{3} Z_{i3}^3 + \frac{1}{4} Z_{i4}^2 + \varepsilon_y)
#' }
#'
#' and generate:
#'
#' \deqn{
#'   Y_i \sim Bern[1, \Pr(Y_i = 1 \mid \boldsymbol{X})]
#' }
#'
#' where \eqn{\varepsilon_y \sim \mathcal{N}(0, sigma\_y)}.
#'
#' For \code{"poisson"}, we simulate:
#'
#' \deqn{
#'   \lambda_Y = exp(effect \times Z_{i1} + \frac{1}{2} Z_{i2}^2 + \frac{1}{3} Z_{i3}^3 + \frac{1}{4} Z_{i4}^2 + \varepsilon_y)
#' }
#'
#' and generate:
#'
#' \deqn{
#'   Y_i \sim Poisson(\lambda_Y)
#' }
#'
#' \strong{Generating Predictions}
#'
#' To generate predicted outcomes for \code{"mean"} and \code{"quantile"}, we
#' simulate a continuous variable with mean given by the empirical mean of the
#' training data and error variance given by the \code{sigma_y} argument.
#'
#' For \code{"ols"}, we fit a generalized additive model (GAM) on the
#' simulated \emph{training} dataset and calculate predictions for the
#' \emph{labeled} and \emph{unlabeled} datasets as deterministic functions of
#' \eqn{\boldsymbol{X}}. Specifically, we fit the following GAM:
#'
#' \deqn{
#'   Y^{\mathcal{T}} = s_0 + s_1(X_1^{\mathcal{T}}) + s_2(X_2^{\mathcal{T}}) +
#'   s_3(X_3^{\mathcal{T}}) + s_4(X_4^{\mathcal{T}}) + \varepsilon_p,
#' }
#'
#' where \eqn{\mathcal{T}} denotes the \emph{training} dataset, \eqn{s_0} is an
#' intercept term, and \eqn{s_1(\cdot)}, \eqn{s_2(\cdot)}, \eqn{s_3(\cdot)},
#' and \eqn{s_4(\cdot)} are smoothing spline functions for \eqn{X_1}, \eqn{X_2},
#' \eqn{X_3}, and \eqn{X_4}, respectively, with three target equivalent degrees
#' of freedom. Residual error is modeled as \eqn{\varepsilon_p}.
#'
#' Predictions for \emph{labeled} and \emph{unlabeled} datasets are calculated
#' as:
#'
#' \deqn{
#'  f(\boldsymbol{X}^{\mathcal{L}\cup\mathcal{U}}) = \hat{s}_0 + \hat{s}_1(X_1^{\mathcal{L}\cup\mathcal{U}}) +
#' \hat{s}_2(X_2^{\mathcal{L}\cup\mathcal{U}}) + \hat{s}_3(X_3^{\mathcal{L}\cup\mathcal{U}}) +
#' \hat{s}_4(X_4^{\mathcal{L}\cup\mathcal{U}}),
#' }
#'
#' where \eqn{\hat{s}_0, \hat{s}_1, \hat{s}_2, \hat{s}_3}, and \eqn{\hat{s}_4}
#' are estimates of \eqn{s_0, s_1, s_2, s_3}, and \eqn{s_4}, respectively.
#'
#' NOTE: For continuous outcomes, we provide optional arguments \code{shift} and
#' \code{scale} to further apply a location shift and scaling factor,
#' respectively, to the predicted outcomes. These default to \code{shift = 0}
#' and \code{scale = 1}, i.e., no location shift or scaling.
#'
#' For \code{"logistic"}, we train k-nearest neighbors (k-NN) classifiers on
#' the simulated \emph{training} dataset for values of \eqn{k} ranging from 1
#' to 10. The optimal \eqn{k} is chosen via cross-validation, minimizing the
#' misclassification error on the validation folds. Predictions for the
#' \emph{labeled} and \emph{unlabeled} datasets are obtained by applying the
#' k-NN classifier with the optimal \eqn{k} to \eqn{\boldsymbol{X}}.
#'
#' Specifically, for each observation in the \emph{labeled} and \emph{unlabeled}
#' datasets:
#'
#' \deqn{
#'   \hat{Y} = \operatorname{argmax}_c \sum_{i \in \mathcal{N}_k} I(Y_i = c),
#' }
#'
#' where \eqn{\mathcal{N}_k} represents the set of \eqn{k} nearest neighbors in
#' the training dataset, \eqn{c} indexes the possible classes (0 or 1), and
#' \eqn{I(\cdot)} is an indicator function.
#'
#' For \code{"poisson"}, we fit a generalized linear model (GLM) with a log link
#' function to the simulated \emph{training} dataset. The model is of the form:
#'
#' \deqn{
#'   \log(\mu^{\mathcal{T}}) = \gamma_0 + \gamma_1 X_1^{\mathcal{T}} + \gamma_2 X_2^{\mathcal{T}} +
#'   \gamma_3 X_3^{\mathcal{T}} + \gamma_4 X_4^{\mathcal{T}},
#' }
#'
#' where \eqn{\mu^{\mathcal{T}}} is the expected count for the response variable
#' in the \emph{training} dataset, \eqn{\gamma_0} is the intercept, and
#' \eqn{\gamma_1}, \eqn{\gamma_2}, \eqn{\gamma_3}, and \eqn{\gamma_4} are the
#' regression coefficients for the predictors \eqn{X_1}, \eqn{X_2}, \eqn{X_3},
#' and \eqn{X_4}, respectively.
#'
#' Predictions for the \emph{labeled} and \emph{unlabeled} datasets are
#' calculated as:
#'
#' \deqn{
#'   \hat{\mu}^{\mathcal{L} \cup \mathcal{U}} = \exp(\hat{\gamma}_0 + \hat{\gamma}_1 X_1^{\mathcal{L} \cup \mathcal{U}} +
#'   \hat{\gamma}_2 X_2^{\mathcal{L} \cup \mathcal{U}} + \hat{\gamma}_3 X_3^{\mathcal{L} \cup \mathcal{U}} +
#'   \hat{\gamma}_4 X_4^{\mathcal{L} \cup \mathcal{U}}),
#' }
#'
#' where \eqn{\hat{\gamma}_0}, \eqn{\hat{\gamma}_1}, \eqn{\hat{\gamma}_2}, \eqn{\hat{\gamma}_3},
#' and \eqn{\hat{\gamma}_4} are the estimated coefficients.
#'
#' @examples
#'
#' #-- Mean
#'
#' dat_mean <- simdat(c(100, 100, 100), effect = 1, sigma_Y = 1,
#'
#'   model = "mean")
#'
#' head(dat_mean)
#'
#' #-- Linear Regression
#'
#' dat_ols <- simdat(c(100, 100, 100), effect = 1, sigma_Y = 1,
#'
#'   model = "ols")
#'
#' head(dat_ols)
#'
#' @import stats gam caret
#'
#' @export

simdat <- function(n = c(300, 300, 300), effect = 1, sigma_Y = 1,

  model = "ols", shift = 0, scale = 1) {

  #-- CHECK FOR VALID MODEL

  if (!(model %in% c("mean", "quantile", "ols", "logistic", "poisson"))) {

    stop(paste("'model' must be one of c('mean', 'quantile', 'ols',",

      "'logistic', 'poisson')."))
  }

  #-- GENERATE SYSTEMATIC COMPONENT

  if (model %in% c("mean", "quantile")) {

    X <- 1

    mu <- effect * 1

  } else if (model %in% c("ols", "logistic", "poisson")) {

    X <- matrix(rnorm(sum(n) * 4), ncol = 4, nrow = sum(n))

    mu <- effect * X[,1] + (1/2) * X[,2]^2 + (1/3) * X[,3]^3 + (1/4) * X[,4]^2
  }

  #-- GENERATE ERROR COMPONENT

  eps <- rnorm(sum(n), 0, sigma_Y)

  #-- GENERATE OUTCOMES

  if (model %in% c("mean", "quantile", "ols")) {

    Y <- mu + eps

  } else if (model == "logistic") {

    p_Y <- plogis(mu + eps)

    Y <- rbinom(sum(n), 1, p_Y)

  } else if (model == "poisson") {

    lam_Y <- exp(mu + eps)

    Y <- rpois(sum(n), lam_Y)
  }

  #-- CREATE DATA.FRAME

  set_label <- rep(c("training", "labeled", "unlabeled"), n)

  if (model %in% c("mean", "quantile")) {

    dat <- data.frame(Y, f = NA, set_label)

  } else if (model %in% c("ols", "logistic", "poisson")) {

    dat <- data.frame(X, Y, f = NA, set_label)
  }

  #-- GENERATE PREDICTIONS

  if (model %in% c("mean", "quantile")) {

    dat[set_label == "labeled", "f"] <- (

      mean(dat[set_label == "training", "Y"]) +

        rnorm(n[2], 0, sigma_Y) - shift) / scale

    dat[set_label == "unlabeled", "f"] <- (

      mean(dat[set_label == "training", "Y"]) +

        rnorm(n[3], 0, sigma_Y) - shift) / scale

  } else if (model == "ols") {

    fit_gam <- gam::gam(Y ~ gam::s(X1) + gam::s(X2) + gam::s(X3) + gam::s(X4),

      data = dat[set_label == "training",])

    dat[set_label == "labeled", "f"] <- (predict(

      fit_gam, newdat = dat[set_label == "labeled",]) - shift) / scale

    dat[set_label == "unlabeled", "f"] <- (predict(

      fit_gam, newdat = dat[set_label == "unlabeled",]) - shift) / scale

  } else if (model == "logistic") {

    knn_tune <- caret::train(

      factor(Y) ~ X1 + X2 + X3 + X4, data = dat[set_label == "training",],

      method = "knn", trControl = trainControl(method = "cv"),

      tuneGrid = data.frame(k = c(1:10)))

    fit_knn <- caret::knn3(factor(Y) ~ X1 + X2 + X3 + X4,

      data = dat[set_label == "training",], k = knn_tune$bestTune$k)

    dat[set_label == "labeled", "f"] <- predict(

      fit_knn, dat[set_label == "labeled",], type = "class") |>

      as.numeric() - 1

    dat[set_label == "unlabeled", "f"] <- predict(

      fit_knn, dat[set_label == "unlabeled",], type = "class") |>

      as.numeric() - 1

  } else if (model == "poisson") {

    fit_poisson <- glm(Y ~ X1 + X2 + X3 + X4, family = poisson,

      data = dat[set_label == "training",])

    dat[set_label == "labeled", "f"] <- predict(fit_poisson,

      newdata = dat[set_label == "labeled",], type = "response")

    dat[set_label == "unlabeled", "f"] <- predict(fit_poisson,

      newdata = dat[set_label == "unlabeled",], type = "response")
  }

  return(dat)
}

#=== END =======================================================================
