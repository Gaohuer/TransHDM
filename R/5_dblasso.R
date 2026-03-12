#' Fit Debiased LASSO with Transfer Learning
#'
#' Fits a debiased LASSO regression model under transfer learning framework, supporting feature
#' selection and coefficient estimation by combining target and source data.
#'
#' @param target A list containing two elements:
#'   \itemize{
#'     \item{x: Feature matrix of target data (numeric matrix)}
#'     \item{y: Response vector of target data (numeric vector)}
#'   } Required.
#' @param source A list (optional, default: NULL) containing two elements:
#'   \itemize{
#'     \item{x: Feature matrix of source data (numeric matrix)}
#'     \item{y: Response vector of source data (numeric vector)}
#'   } Used when transfer = TRUE.
#' @param transfer A logical value (default: FALSE) indicating whether to enable transfer learning
#'   (combining source data with target data for estimation).
#' @param level A numeric value (default: 0.95) specifying confidence level for confidence intervals.
#' @param lambda A string specifying criterion for selecting regularization parameter:
#'   \itemize{
#'     \item{'lambda.min': Lambda value giving minimum cross-validation error}
#'     \item{'lambda.1se': Largest lambda within 1 standard error of minimum error}
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item{dbcoef.hat: Debiased LASSO coefficient vector (including intercept)}
#'   \item{coef.hat: Original LASSO coefficient vector}
#'   \item{CI: Data frame with confidence intervals (lb = lower bound, ub = upper bound)}
#'   \item{var.est: Variance estimates for debiased coefficients}
#'   \item{se.est: Standard errors for debiased coefficients}
#'   \item{P.value: Vector of p-values for coefficients}
#' }
#'
#' @examples
#' # Prepare target and source data
#' target <- list(x = matrix(rnorm(100 * 20), 100, 20), y = rnorm(100))
#' source <- list(x = matrix(rnorm(200 * 20), 200, 20), y = rnorm(200))
#'
#' # Non-transfer mode
#' result_no_transfer <- dblasso(target = target, transfer = FALSE,
#'                              level = 0.95, lambda = 'lambda.min')
#' summary(result_no_transfer)
#'
#' # Transfer learning mode
#' result_transfer <- dblasso(target = target, source = source, transfer = TRUE,
#'                           level = 0.95, lambda = 'lambda.min')
#' summary(result_transfer)
#' @import glmnet
#' @import foreach
#' @export
dblasso <- function(
    target,
    source = NULL,
    transfer = FALSE,
    level = 0.95,
    lambda = "lambda.1se"
) {

  ## ---------- internal: nodewise regression ----------
  gamma_est <- function(target, source = NULL, transfer = FALSE) {

    p_x <- ncol(target$x)

    if (!transfer) {

      target_data <- data.frame(Y = target$y, target$x)
      colnames(target_data) <- c("Y", paste0("X", seq_len(p_x)))

      cv_glm <- glmnet::cv.glmnet(
        x = as.matrix(target_data[, -1, drop = FALSE]),
        y = target_data$Y,
        alpha = 1,
        intercept = FALSE,
        family = "gaussian"
      )

      beta_TL <- as.vector(
        cv_glm$glmnet.fit$beta[, which(cv_glm$lambda == cv_glm[[lambda]])]
      )

    } else {

      target_data <- data.frame(Y = target$y, target$x)
      source_data <- data.frame(Y = source$y, source$x)

      colnames(target_data) <- colnames(source_data) <-
        c("Y", paste0("X", seq_len(p_x)))

      data_all <- rbind(target_data, source_data)

      cv_glm <- glmnet::cv.glmnet(
        x = as.matrix(data_all[, -1, drop = FALSE]),
        y = data_all$Y,
        alpha = 1,
        intercept = FALSE,
        family = "gaussian"
      )

      w <- as.vector(
        cv_glm$glmnet.fit$beta[, which(cv_glm$lambda == cv_glm[[lambda]])]
      )

      offset <- as.numeric(as.matrix(target_data[, -1, drop = FALSE]) %*% w)

      cv_glm <- glmnet::cv.glmnet(
        x = as.matrix(target_data[, -1, drop = FALSE]),
        y = target_data$Y,
        alpha = 1,
        intercept = FALSE,
        offset = offset,
        family = "gaussian"
      )

      delta <- as.vector(
        cv_glm$glmnet.fit$beta[, which(cv_glm$lambda == cv_glm[[lambda]])]
      )

      beta_TL <- w + delta
    }

    beta_TL
  }

  ## ---------- CI level ----------
  r.level <- level + (1 - level) / 2
  j<-NULL

  ## ---------- initial lasso ----------
  coef.hat <- as.vector(
    lasso(target = target, source = source, transfer = transfer)
  )

  ## ---------- build design matrix with FIXED names ----------
  orig_xnames <- colnames(target$x)
  if (is.null(orig_xnames)) {
    orig_xnames <- paste0("X", seq_len(ncol(target$x)))
  }

  D <- list(target = target, source = source)

  D$target$x <- cbind("(Intercept)" = 1, as.matrix(target$x))
  colnames(D$target$x) <- c("(Intercept)", orig_xnames)

  if (transfer) {
    D$source$x <- cbind("(Intercept)" = 1, as.matrix(source$x))
    colnames(D$source$x) <- colnames(D$target$x)
  }

  coef.name <- colnames(D$target$x)

  ## ---------- dimensions ----------
  p <- ncol(D$target$x)

  ## ---------- combined design ----------
  X.comb <- if (transfer) {
    rbind(D$source$x, D$target$x)
  } else {
    D$target$x
  }

  ## ---------- Sigma ----------
  Sigma.hat <- crossprod(X.comb) / nrow(X.comb)

  ## ---------- Theta ----------
  L <- foreach::foreach(j = seq_len(p), .combine = "rbind") %do% {

    D1 <- D
    D1$target$y <- D1$target$x[, j]
    D1$target$x <- D1$target$x[, -j, drop = FALSE]

    if (transfer) {
      D1$source$y <- D1$source$x[, j]
      D1$source$x <- D1$source$x[, -j, drop = FALSE]
      gamma <- gamma_est(D1$target, D1$source, TRUE)
    } else {
      gamma <- gamma_est(D1$target)
    }

    tau2 <- Sigma.hat[j, j] -
      Sigma.hat[j, -j, drop = FALSE] %*% gamma

    theta <- numeric(p)
    theta[j] <- 1
    theta[-j] <- -gamma

    c(theta, tau2)
  }

  Theta.hat <- solve(diag(L[, p + 1])) %*% L[, 1:p]

  ## ---------- debiased estimator ----------
  residuals <- D$target$y - D$target$x %*% coef.hat

  dbcoef.hat <- as.vector(
    coef.hat +
      Theta.hat %*% crossprod(D$target$x, residuals) / length(residuals)
  )

  ## ---------- inference ----------
  var.est <- diag(Theta.hat %*% Sigma.hat %*% t(Theta.hat))
  se.est  <- sqrt(var.est / length(residuals))

  CI <- data.frame(
    lb = dbcoef.hat - qnorm(r.level) * se.est,
    ub = dbcoef.hat + qnorm(r.level) * se.est
  )

  Z <- dbcoef.hat / se.est
  P.value <- 2 * pnorm(abs(Z), lower.tail = FALSE)

  ## ---------- assign names ----------
  names(coef.hat)   <- coef.name
  names(dbcoef.hat) <- coef.name
  rownames(CI)      <- coef.name
  names(var.est)    <- coef.name
  names(se.est)     <- coef.name
  names(P.value)    <- coef.name

  res <- list(
    dbcoef.hat = dbcoef.hat,
    coef.hat   = coef.hat,
    CI         = CI,
    var.est    = var.est,
    se.est     = se.est,
    P.value    = P.value
  )

  class(res) <- "dblasso"
  res
}
