#' Data generation functions for OLS regression example
#'
#' @param n vector of size 3 indicating the sample size in the training,
#' labeled/test, and unlabeled/validation data sets
#'
#' @param beta1 first regression coefficient (or, regression coefficient of
#' variable of interest for inference)
#'
#  @inheritParams gam::gam
#'
#' @return A data frame containing 4 columns: labeled outcome, predicted
#' outcome and a character variable indicating which data set the observation
#' belongs to (training, test/labeled, validation/unlabeled).
#'
#' @examples
#'
#' # Return a stacked data set with 100 observations for each individual sets
#'
#' simdat(c(100, 100, 100), 1)
#'
#' @import stats gam
#'
#' @export

simdat <- function(n = c(300, 300, 300), beta1 = 1) {

  X1 <- rnorm(sum(n), 1)
  X2 <- rnorm(sum(n), 1)
  X3 <- rnorm(sum(n), 1)
  X4 <- rnorm(sum(n), 2)

  Y <- c(beta1*X1 + 0.5*X2 + 3*smooth(X3) + 4*smooth(X4) + rnorm(sum(n)))

  set <- rep(c("training", "labeled", "unlabeled"), n)

  dat <- data.frame(X1, X2, X3, X4, Y, Yhat = NA, set)

  fit_gam <- gam::gam(

    Y ~ s(X1) + s(X2) + s(X3) + s(X4), data = dat[set == "training",])

  dat[set == "labeled", "Yhat"] <- predict(

    fit_gam, newdat = dat[set == "labeled",])

  dat[set == "unlabeled", "Yhat"] <- predict(

    fit_gam, newdat = dat[set == "unlabeled",])

  return(dat)
}

#' Data generation for logistic regression examples
#'
#' @param n vector of size 3 indicating the sample size in the training,
#' labeled/test, and unlabeled/validation data sets
#'
#' @param betac first regression coefficient (or, regression coefficient of
#' variable of interest for inference)
#'
#' @return A data frame containing 4 columns: labeled outcome, predicted
#' outcome and a character variable indicating which data set the observation
#' belongs to (training, test/labeled, validation/unlabeled).
#'
#' @examples
#'
#' # Return a stacked data set with 100 observations for each individual sets
#'
#' simdat(c(100, 100, 100), 1)
#'
#' @import stats caret
#'
#' @export

simdat_logistic <- function(n = c(300, 300, 300), betac = 1) {

  X1 <- rnorm(sum(n), 1)
  X2 <- rnorm(sum(n), 2)
  Xc <- sample(0:2, sum(n), replace = T)

  p_Y <- plogis(c(betac*as.numeric(Xc == 2) +

    1*as.numeric(Xc == 1) + 1*smooth(X1) - 2*smooth(X2) + rnorm(sum(n))))

  Y <- rbinom(sum(n), 1, p_Y)

  set <- rep(c("training", "labeled", "unlabeled"), n)

  dat <- data.frame(X1, X2, Xc, Y, Yhat = NA, set)

  knn_tune <- caret::train(

    factor(Y) ~ X1 + X2 + Xc, data = dat[set == "training",],

    method = "knn", trControl = trainControl(method = "cv"),

    tuneGrid = data.frame(k = c(1:10)))

  fit_knn <- caret::knn3(factor(Y) ~ X1 + X2 + Xc,

    data = dat[set == "training",], k = knn_tune$bestTune$k)

  dat[set == "labeled", "Yhat"] <- predict(

    fit_knn, dat[set == "labeled",], type = "class")

  dat[set == "unlabeled", "Yhat"] <- predict(

    fit_knn, dat[set == "unlabeled",], type = "class")

  return(dat)
}

#=== END =======================================================================
