#' Data generation for OLS examples
#'
#' @param n vector of size 3 indicating the sample size in the training,
#' labelled, and unlabelled data sets
#'
#' @param beta1 first regression coefficient (or, regression coefficient of
#' variable of interest for inference)
#'
#  @inheritParams gam::gam
#'
#' @return A data frame containing 4 regressors, labelled outcome, predicted
#' outcome and a character variable indicating which dat set the observartion
#' belongs to (training, test, validation).
#'
#' @examples
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

  set <- rep(c("trn", "tst", "val"), n)

  dat <- data.frame(X1, X2, X3, X4, Y, Yhat = NA, set)

  fit_gam <- gam::gam(

    Y ~ s(X1) + s(X2) + s(X3) + s(X4), data = dat[set == "trn",])

  dat[set == "tst", "Yhat"] <- predict(fit_gam, newdat = dat[set == "tst",])

  dat[set == "val", "Yhat"] <- predict(fit_gam, newdat = dat[set == "val",])

  return(dat)
}

#' Data generation for logistic regression examples
#'
#' @param n vector of size 3 indicating the sample size in the training,
#' labelled, and unlabelled data sets
#'
#' @param betac first regression coefficient (or, regression coefficient of
#' variable of interest for inference)
#'
#' @return A data frame containing 4 regressors, labelled outcome, predicted
#' outcome and a character variable indicating which dat set the observartion
#' belongs to (training, test, validation).
#'
#' @examples
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

  p_Y <- plogis(c(betac*as.numeric(Xc == 2) + 1*as.numeric(Xc == 1) +

                    1*smooth(X1) - 2*smooth(X2) + rnorm(sum(n))))

  Y <- rbinom(sum(n), 1, p_Y)    #factor()

  set <- rep(c("trn", "tst", "val"), n)

  dat <- data.frame(X1, X2, Xc, Y, Yhat = NA, set)

  knn_tune <- caret::train(factor(Y) ~ X1 + X2 + Xc, data = dat[set == "trn",],

                           method = "knn", trControl = trainControl(method = "cv"),

                           tuneGrid = data.frame(k = c(1:10)))

  fit_knn <- caret::knn3(factor(Y) ~ X1 + X2 + Xc, data = dat[set == "trn",],

                         k = knn_tune$bestTune$k)

  dat[set == "tst", "Yhat"] <- predict(

    fit_knn, dat[set == "tst",], type = "class")

  dat[set == "val", "Yhat"] <- predict(

    fit_knn, dat[set == "val",], type = "class")

  return(dat)
}

#=== END =======================================================================
