#===============================================================================
#  DATA GENERATION FUNCTIONS
#===============================================================================

#--- DATA GENERATION FOR OLS REGRESSION ----------------------------------------

#' Data generation functions for OLS regression example
#'
#' @param n (int vector) vector of size 3 indicating the sample size in the
#' training, labeled, and unlabeled data sets
#'
#' @param effect (float) regression coefficient for the first variable of
#' interest for inference (defaults to 1)
#'
#' @param sigma_Y (float) residual variance for the generated outcome
#'
#  @inheritParams gam::gam
#'
#' @return A data.frame containing n rows and seven columns corresponding to
#' the labeled outcome (Y), the predicted outcome (f), a character variable
#' (set) indicating which data set the observation belongs to (training,
#' labeled, or unlabeled), and four independent, normally distributed
#' predictors (X1, X2, X3, and X4).
#'
#' @examples
#'
#' # Return a stacked data set with 100 observations for each individual sets
#'
#' dat <- simdat(c(100, 100, 100), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' @import stats gam
#'
#' @export

simdat <- function(n = c(300, 300, 300), effect = 1, sigma_Y = 1) {

  eps <- rnorm(sum(n), 0, sigma_Y)

  X <- matrix(rnorm(sum(n) * 4), ncol = 4, nrow = sum(n))

  Y <- effect * X[,1] + (1/2) * X[,2]^2 + 3 * X[,3]^3 + 4 * X[,4]^2 + eps

  set <- rep(c("training", "labeled", "unlabeled"), n)

  dat <- data.frame(X, Y, f = NA, set)

  fit_gam <- gam::gam(Y ~ gam::s(X1) + gam::s(X2) + gam::s(X3) + gam::s(X4),

    data = dat[set == "training",])

  dat[set == "labeled", "f"] <- predict(

    fit_gam, newdat = dat[set == "labeled",])

  dat[set == "unlabeled", "f"] <- predict(

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

  dat <- data.frame(X1, X2, Xc, Y, f = NA, set)

  knn_tune <- caret::train(

    factor(Y) ~ X1 + X2 + Xc, data = dat[set == "training",],

    method = "knn", trControl = trainControl(method = "cv"),

    tuneGrid = data.frame(k = c(1:10)))

  fit_knn <- caret::knn3(factor(Y) ~ X1 + X2 + Xc,

    data = dat[set == "training",], k = knn_tune$bestTune$k)

  dat[set == "labeled", "f"] <- predict(

    fit_knn, dat[set == "labeled",], type = "class")

  dat[set == "unlabeled", "f"] <- predict(

    fit_knn, dat[set == "unlabeled",], type = "class")

  return(dat)
}

#=== END =======================================================================
