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
#' \code{"mean"}, \code{"quantile"}, \code{"ols"}, or \code{"logistic"}.
#' Default is \code{"ols"}.
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
#' (set) indicating which data set the observation belongs to (training,
#' labeled, or unlabeled), and four independent, normally distributed
#' predictors (X1, X2, X3, and X4), where applicable.
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

  set <- rep(c("training", "labeled", "unlabeled"), n)

  if (model %in% c("mean", "quantile")) {

    dat <- data.frame(Y, f = NA, set)

  } else if (model %in% c("ols", "logistic", "poisson")) {

    dat <- data.frame(X, Y, f = NA, set)
  }

  #-- GENERATE PREDICTIONS

  if (model %in% c("mean", "quantile")) {

    dat[set == "labeled", "f"] <- (mean(dat[set == "training", "Y"]) +

      rnorm(n[2], 0, sigma_Y) - shift) / scale

    dat[set == "unlabeled", "f"] <- (mean(dat[set == "training", "Y"]) +

      rnorm(n[3], 0, sigma_Y) - shift) / scale

  } else if (model == "ols") {

    fit_gam <- gam::gam(Y ~ gam::s(X1) + gam::s(X2) + gam::s(X3) + gam::s(X4),

      data = dat[set == "training",])

    dat[set == "labeled", "f"] <- (predict(

      fit_gam, newdat = dat[set == "labeled",]) - shift) / scale

    dat[set == "unlabeled", "f"] <- (predict(

      fit_gam, newdat = dat[set == "unlabeled",]) - shift) / scale

  } else if (model == "logistic") {

    knn_tune <- caret::train(

      factor(Y) ~ X1 + X2 + X3 + X4, data = dat[set == "training",],

      method = "knn", trControl = trainControl(method = "cv"),

      tuneGrid = data.frame(k = c(1:10)))

    fit_knn <- caret::knn3(factor(Y) ~ X1 + X2 + X3 + X4,

      data = dat[set == "training",], k = knn_tune$bestTune$k)

    dat[set == "labeled", "f"] <- predict(

      fit_knn, dat[set == "labeled",], type = "class") |> as.numeric() - 1

    dat[set == "unlabeled", "f"] <- predict(

      fit_knn, dat[set == "unlabeled",], type = "class") |> as.numeric() - 1

  } else if (model == "poisson") {

    fit_poisson <- glm(Y ~ X1 + X2 + X3 + X4, family = poisson,

      data = dat[set == "training",])

    dat[set == "labeled", "f"] <- predict(fit_poisson,

      newdata = dat[set == "labeled",], type = "response")

    dat[set == "unlabeled", "f"] <- predict(fit_poisson,

      newdata = dat[set == "unlabeled",], type = "response")
  }

  return(dat)
}

#=== END =======================================================================
