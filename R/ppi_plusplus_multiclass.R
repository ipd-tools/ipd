#===============================================================================
# PPI++ MULTICLASS LOGISTIC REGRESSION
#===============================================================================

library(nnet)                                                                   # check these
library(Matrix)
library(stats)

#=== HELPER FUNCTIONS ==========================================================

#--- ONE-HOT ENCODING ----------------------------------------------------------

one_hot <- function(Y, classes = unique(Y)) {

  n <- nrow(Y)

  K <- length(classes)

  Y_one_hot <- matrix(0, nrow = n, ncol = K, dimnames = list(NULL, classes))

  for (i in 1:n) {

    class_index <- match(Y[i], classes)

    Y_one_hot[i, class_index] <- 1
  }

  return(Y_one_hot)
}

#--- SOFTMAX -------------------------------------------------------------------

softmax <- function(theta, X, K) {

  n <- nrow(X)

  p <- ncol(X)

  numer <- matrix(0, nrow = n, ncol = K - 1)

  for (class in 1:(K - 1)) {

    idx <- 1:p + (class - 1)*p

    numer[, class] <- exp(X %*% theta[idx])
  }

  denom <- rowSums(numer) + 1

  probs <- cbind(numer, 1) / denom

  return(probs)
}

#--- W FOR MULTICLASS HESSIAN MATRIX -------------------------------------------

block_W <- function(probs) {

  k <- ncol(probs)
  n <- nrow(probs)
  N <- n*(k - 1)

  W <- lapply(1:(k - 1), function(class) {

    diag(probs[, class]*(1 - probs[, class])) }) |>

    Matrix::bdiag() |> as.matrix()

  if (k > 2) {

    for (i in 1:(k - 2)) {

      for(j in  c((i + 1):(k - 1))) {

        row_idx <- n*(i - 1) + 1

        col_idx <- n*(j - 1) + 1

        W[row_idx:(row_idx + n - 1),
          col_idx:(col_idx + n - 1)] <- diag(-probs[,i]*probs[,j])
      }
    }

    W <- W + t(W)

    diag(W) <- diag(W)/2
  }

  return(W)
}

#=== PPI++ MULTICLASS LOGISTIC REGRESSION POINT ESTIMATE =======================

#' PPI++ Multiclass Logistic Regression Point Estimate
#'
#' @description
#' Computes the prediction-powered point estimate of the multiclass logistic
#' regression coefficients.
#'
#' @param X_l (matrix): n x p matrix of covariates in the labeled data.
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param X_u (matrix): N x p matrix of covariates in the unlabeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param lhat (float, optional): Power-tuning parameter.
#' The default value `NULL` will estimate the optimal value from data.
#' Setting `lhat = 1` recovers PPI with no power tuning.
#' Setting `lhat = 0` recovers the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`.
#' If `None`, it optimizes the total variance over all coordinates.
#' Must be in {1, ..., p} where p is the shape of the estimand.
#'
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
#'
#' @returns (vector): p-vector of prediction-powered point estimates of the
#' logistic regression coefficients.
#'
#' @examples
#'
#' ## NEED EXAMPLES
#'
#' @import stats
#'
#' @export

ppi_plusplus_multiclass_est <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, opts = NULL) {                                      # IMPLEMENT WEIGHTS?

  n <- nrow(Y_l)

  p <- ncol(X_l)

  N <- nrow(f_u)

  ##
  Y_l <- factor(Y_l) |> as.numeric() |> matrix(ncol = 1)                         # Do this better?
  f_l <- factor(f_l) |> as.numeric() |> matrix(ncol = 1)
  f_u <- factor(f_u) |> as.numeric() |> matrix(ncol = 1)
  ##

  classes <- sort(unique(Y_l))

  K <- length(classes)

  if (is.null(opts) || !("factr" %in% names(opts))) {

    opts <- list(factr = 1e-15)
  }

  #-- Rectified Multiclass Logistic Regression Loss Function

  rectified_multiclass_loss <- function(theta) {

    #-- Default Reference Class: 0

    theta[1:p] <- 0                                                           # Put this back?

    EY_mat <- matrix(0, K * p, K * p)

    for (k in 1:K) {

      Ey <- matrix(0, p, K * p)

      Ey[, ((k - 1) * p + 1):(k * p)] <- diag(p)

      EY_mat[((k - 1) * p + 1):(k * p), ] <- Ey
    }

    EY <- array(EY_mat, dim = c(K, p, K * p))

    loss0 <- loss1 <- loss2 <- 0

    for (i in 1:n) {

      y_l <- Y_l[i]

      Xi_l <- X_l[i, ]

      Ey <- matrix(EY[y_l, , ], p, K * p, byrow = T)

      loss0 <- loss0 -(Xi_l %*% Ey) %*% theta +

        log(sum(exp(apply(EY, 3, function(ey) Xi_l %*% matrix(ey, p, K)) %*% theta)))

      yhat <- f_l[i]

      Eyhat <- matrix(EY[yhat, , ], p, K * p, byrow = T)

      loss1 <- loss1 -(Xi_l %*% Eyhat) %*% theta +

        log(sum(exp(apply(EY, 3, function(ey) Xi_l %*% matrix(ey, p, K)) %*% theta)))
    }

    for (i in 1:N) {

      y_unlabeled <- f_u[i]

      Ey_unlabeled <- matrix(EY[y_unlabeled, , ], p, K * p, byrow = T)

      Xi_u <- X_u[i, ]

      loss2 <- loss2 -(Xi_u %*% Ey_unlabeled) %*% theta +

        log(sum(exp(apply(EY, 3, function(ey) Xi_u %*% matrix(ey, p, K)) %*% theta)))
    }

    loss <- (1 / n) * loss0 - lhat_curr / n * loss1 + lhat_curr / N * loss2

    return(loss)
  }

  # Helper function to compute gradient of the parameter

  gradient <- function(X, Y, theta, K, classes) {

    Y_onehot <- one_hot(Y, classes)

    theta_2D <- matrix(theta, nrow = K, byrow = T)

    # by default class 0 is the reference class

    theta_2D[1, ] <- rep(0, ncol(X))

    a <- exp(X %*% t(theta_2D))

    P <- a / rowSums(a)

    gd <- - t(Y_onehot - P) %*% X

    return(as.vector(t(gd)))
  }

  #-- Rectified Multiclass Logistic Regression Gradient

  rectified_multiclass_grad <- function(theta) {

    gd_rect <- gradient(X_u, f_u, theta, K, classes) * lhat_curr / N -

    gradient(X_l, f_l, theta, K, classes) * lhat_curr / n +

    gradient(X_l, Y_l, theta, K, classes) / n

    return(gd_rect)
  }

  #-- Initialize Theta

  fit_0 <- multinom(factor(Y_l) ~ X_l - 1, trace = F)

  # theta_0 <- rbind(rep(0, p), coefficients(fit_0))

  theta_0 <- rbind(rep(0, p), matrix(1/(K*p), nrow = K - 1, ncol = p))# switch back once QC'd

  lhat_curr <- ifelse(is.null(lhat), 1, lhat)

  #-- Optimize Over (K-1)*p DF - Forcing Reference Class Parameters to Zero

  est_full <- optim(par = c(t(theta_0)), fn = rectified_multiclass_loss,         # better way to do par?

    gr = rectified_multiclass_grad, method = "L-BFGS-B",

    control = opts)$par

  #-- Remove Reference Class Parameters

  cat("est_full:", est_full, "\n")

  est <- matrix(est_full, nrow = K, byrow = TRUE)[-1, , drop = FALSE] |> t() |> c()  #better way to do this?

  cat("est:\n")

  print(est)

  if (is.null(lhat)) {

    stats <- multiclass_get_stats(est, X_l, Y_l, f_l, X_u, f_u)

    print("stats$grads:")

    print(stats$grads)

    print("stats$inv_hessian:")

    print(stats$inv_hessian)

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, clip = T)

    cat("lhat:", lhat, "\n")

    return(ppi_plusplus_multiclass_est(X_l, Y_l, f_l, X_u, f_u,

      opts = opts, lhat = lhat, coord = coord))

  } else {

    return(est)
  }
}


################################################################################
################################################################################


multiclass_get_stats <- function(est, X_l, Y_l, f_l, X_u, f_u, use_u = TRUE) {   # Implement weights?

  n <- nrow(Y_l)

  p <- ncol(X_l)

  N <- nrow(f_u)

  classes <- sort(unique(Y_l))

  K <- length(classes)

  hessian <- matrix(0, nrow = (K - 1) * p, ncol = (K - 1) * p)

  grads_hat_unlabeled <- matrix(0, nrow = N, ncol = (K - 1) * p)

  if (use_u) {

    for (i in 1:N) {

      Xi_u <- matrix(X_u[i,], nrow = 1)

      Yi_hat_unlabeled <- f_u[i]

      probs <- softmax(est, Xi_u, K)

      probs_vec <- probs[, 1:(K - 1)]

      probs_vec <- matrix(probs_vec, ncol = 1, byrow = TRUE)

      y_unlabeled <- sapply(classes[-K],

        function(class) as.numeric(Yi_hat_unlabeled == class)) |>

        matrix(ncol = 1, byrow = TRUE)

      X_bmat <- replicate(K - 1, Xi_u, simplify = F) |>

        Matrix::bdiag() |> as.matrix()

      W_bmat <- block_W(probs)

      grads_i <- t(X_bmat) %*% (y_unlabeled - probs_vec)

      grads_hat_unlabeled[i, ] <- grads_i
    }
  }

  grads <- matrix(0, nrow = n, ncol = (K - 1) * p)

  grads_hat <- matrix(0, nrow = n, ncol = (K - 1) * p)

  for (i in 1:n) {

    Xi_l <- matrix(X_l[i,], nrow = 1)

    Yi_l <- Y_l[i]

    Yi_hat <- f_l[i]

    probs <- softmax(est, Xi_l, K)

    probs_vec <- probs[, 1:(K - 1)]

    probs_vec <- matrix(probs_vec, ncol = 1, byrow = TRUE)

    y <- sapply(classes[-K], function(class) as.numeric(Yi_l == class)) |>

      matrix(ncol = 1, byrow = TRUE)

    if (i == 2) {

      print("est:")

      print(est)

      print("Xi_l:")

      print(Xi_l)

      print("Yi_l:")

      print(Yi_l)

      print("probs:")

      print(probs)

      print("y:")

      print(y)
    }

    yhat <- sapply(classes[-K], function(class) as.numeric(Yi_hat == class)) |>

      matrix(ncol = 1, byrow = TRUE)

    X_bmat <- replicate(K - 1, Xi_l, simplify = F) |>

      Matrix::bdiag() |> as.matrix()

    W_bmat <- block_W(probs)

    grads[i, ] <- t(X_bmat) %*% (y - probs_vec)

    grads_hat[i, ] <- t(X_bmat) %*% (yhat - probs_vec)
  }

  # Function to compute Gradients and Hessian
  multiclass_logistic_get_stats <- function(Y, X, classes, K, theta) {

    probs <- softmax(theta, X, K)

    probs_vec <- probs[, 1:(K - 1)]

    probs_vec <- matrix(probs_vec, ncol = 1, byrow = TRUE)

    y <- sapply(classes[1:(K - 1)], function(class) (Y == class) + 0)

    y <- matrix(y, ncol = 1, byrow = TRUE)

    X_bmat <- replicate(K - 1, X, simplify = F) |>

      Matrix::bdiag() |> as.matrix()

    W_bmat <- block_W(probs)

    grads <- t(X_bmat) %*% (y - probs_vec)

    hessian <- -t(X_bmat) %*% W_bmat %*% X_bmat

    return(list(grads = grads, hessian = hessian))
  }

  # Stats of labeled data
  stats_X_Y <- multiclass_logistic_get_stats(Y_l, X_l, classes, K, est)

  if (use_u) {
    stats_unlabeled <- multiclass_logistic_get_stats(f_u, X_u, classes, K, est)
    hessian <- 1 / (n + N) * (stats_X_Y$hessian + stats_unlabeled$hessian)
  } else {
    hessian <- stats_X_Y$hessian
  }

  inv_hessian <- solve(hessian)

  return(list(grads = grads, grads_hat = grads_hat, grads_hat_unlabeled = grads_hat_unlabeled, hessian = hessian, inv_hessian = inv_hessian))
}


##################### CHECK ABOVE ##############################################

calc_lhat_glm <- function(grads, grads_hat, grads_hat_unlabeled, inv_hessian,

  coord = NULL, clip = FALSE) {

  if (is.null(dim(grads))) {

    grads <- matrix(grads, ncol = 1)
  }

  if (is.null(dim(grads_hat))) {

    grads_hat <- matrix(grads_hat, ncol = 1)
  }

  if (is.null(dim(grads_hat_unlabeled))) {

    grads_hat_unlabeled <- matrix(grads_hat_unlabeled, ncol = 1)
  }

  n <- nrow(grads)

  N <- nrow(grads_hat_unlabeled)

  p <- ncol(inv_hessian)

  cov_grads <- matrix(0, nrow = p, ncol = p)

  for (i in 1:n) {

    cov_grads <- cov_grads + (1/n) * (outer(grads[i,] - colMeans(grads),

      grads_hat[i,] - colMeans(grads_hat)) +

        outer(grads_hat[i,] - colMeans(grads_hat),

          grads[i,] - colMeans(grads)))
  }

  var_grads_hat <- cov(rbind(grads_hat, grads_hat_unlabeled))

  if (is.null(coord)) {

    vhat <- inv_hessian

  } else {

    vhat <- inv_hessian %*% diag(p)[, coord]
  }

  if (p > 1) {

    num <- ifelse(is.null(coord),

      sum(diag(vhat %*% cov_grads %*% vhat)), vhat %*% cov_grads %*% vhat)

    denom <- ifelse(is.null(coord),

      2*(1 + n/N) * sum(diag(vhat %*% var_grads_hat %*% vhat)),

      2*(1 + n/N) * vhat %*% var_grads_hat %*% vhat)

  } else {

    num <- vhat * cov_grads * vhat

    denom <- 2*(1 + n/N) * vhat * var_grads_hat * vhat
  }

  lhat <- num / denom

  if (clip) {

    lhat <- pmax(0, pmin(lhat, 1))
  }

  return(as.numeric(lhat))
}


#=== PPI++ MULTICLASS LOGISTIC REGRESSION ======================================

#' PPI++ Multiclass Logistic Regression Estimator and Inference
#'
#' @description
#' Computes the prediction-powered confidence interval for the multiclass
#' logistic regression coefficients using the PPI++ algorithm.
#'
#' @param X_l (matrix): n x p matrix of covariates in the labeled data.
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param X_u (matrix): N x p matrix of covariates in the unlabeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param lhat (float, optional): Power-tuning parameter.
#' The default value `NULL` will estimate the optimal value from data.
#' Setting `lhat = 1` recovers PPI with no power tuning.
#' Setting `lhat = 0` recovers the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`.
#' If `None`, it optimizes the total variance over all coordinates.
#' Must be in {1, ..., p} where p is the shape of the estimand.
#'
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
#'
#' @returns (list): A list containing the following:
#'
#' \describe{
#'    \item{est}{(vector): p-vector of PPI++ logistic regression coefficient
#'    estimates.}
#'    \item{se}{(vector): p-vector of standard errors of the coefficients.}
#'    \item{lambda}{(float): estimated power-tuning parameter.}
#' }
#'
#' @examples
#'
#' ##NEED EXAMPLES                                                               # Need examples
#'
#' @import stats
#'
#' @export

ppi_plusplus_multiclass <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, opts = NULL) {                                      # Implement weights?

  n <- nrow(Y_l)

  p <- ncol(X_l)

  N <- nrow(f_u)

  classes <- sort(unique(Y_l), decreasing = TRUE)

  K <- length(classes)

  df <- (K - 1) * p

  use_u <- is.null(lhat) || lhat != 0

  est <- ppi_plusplus_multiclass_est(X_l, Y_l, f_l, X_u, f_u,

    lhat = lhat, coord = coord, opts = opts)

  stats <- multiclass_get_stats(est, X_l, Y_l, f_l, X_u, f_u, use_u = use_u)

  if (is.null(lhat)) {

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, clip = TRUE)

    return(ppi_plusplus_multiclass(X_l, Y_l, f_l, X_u, f_u,

      lhat = lhat, coord = coord, opts = opts))
  }

  var_u <- cov(lhat * t(stats$grads_hat_unlabeled))

  var_l <- cov(t(stats$grads) - lhat * t(stats$grads_hat))

  Sigma_hat <- stats$inv_hessian %*% (n/N * var_u + var_l) %*% stats$inv_hessian

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n), lambda = lhat))
}

#=== END =======================================================================


## QC

X_l <- read.csv("C:/Users/saler/Downloads/X.csv", header = F) |>

  as.matrix()

Y_l <- read.csv("C:/Users/saler/Downloads/Y.csv", header = F) |>

  as.matrix()

f_l <- read.csv("C:/Users/saler/Downloads/Yhat.csv", header = F) |>

  as.matrix()

X_u <- read.csv("C:/Users/saler/Downloads/X_unlabeled.csv", header = F) |>

  as.matrix()

f_u <- read.csv("C:/Users/saler/Downloads/Yhat_unlabeled.csv", header = F) |>

  as.matrix()

ppi_plusplus_multiclass_est(X_l, Y_l, f_l, X_u, f_u)

ppi_plusplus_multiclass(X_l, Y_l, f_l, X_u, f_u)
