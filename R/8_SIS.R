#' Sure Independence Screening for High-Dimensional Mediation Analysis
#'
#' This function performs dimension reduction for high-dimensional mediation analysis
#' using Sure Independence Screening (SIS). Mediators are ranked based on the product
#' of their marginal associations with the exposure and the outcome, and the top-ranked
#' mediators are retained for downstream analysis.
#'
#' The function supports transfer learning, allowing information from a source dataset
#' to be leveraged to improve screening stability and robustness in the target dataset.
#'
#' @param target_data A data frame containing the target dataset. All variables must be numeric.
#' @param source_data A list of data frames containing source datasets (optional, default: NULL).
#'   All variables must be numeric and have the same column names as target_data.
#' @param Y Character string specifying the outcome variable name.
#' @param D Character string specifying the exposure (treatment) variable name.
#' @param M Character vector specifying mediator variable names.
#' @param X Character vector specifying covariate variable names.
#' @param topN An integer specifying the number of mediators to retain after screening.
#' If \code{NULL}, the number is automatically determined as
#' \eqn{\lceil 2n / \log(n) \rceil}, where \eqn{n} is the target sample size.
#' @param transfer A logical value (default: FALSE) indicating whether to apply transfer
#' learning by incorporating the source dataset in the screening procedure.
#' @param verbose A logical value (default: TRUE) controlling whether progress messages
#' are printed to the console.
#' @param ncore An integer (default: 1) specifying the number of CPU cores for parallel
#' computation.
#' @param dblasso_method A logical value (default: FALSE). If TRUE, the debiased lasso
#' (dblasso) is used to estimate marginal effects. If FALSE, standard linear or generalized
#' linear models are used.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{target_SIS}: A data frame containing the outcome, exposure, selected
#'   mediators, and covariates from the target dataset.
#'   \item \code{source_SIS}: A data frame containing the same variables from the source
#'   dataset if \code{transfer = TRUE}; otherwise \code{NULL}.
#'   \item \code{M_ID_name_SIS}: A character vector of selected mediator names.
#' }
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' # Target data
#' M_target <- matrix(rnorm(200 * 50), nrow = 200)
#' colnames(M_target) <- paste0("M", 1:50)
#'
#' target_data <- data.frame(
#'   Y = rnorm(200),
#'   D = rnorm(200),
#'   M_target,
#'   X1 = rnorm(200)
#' )
#'
#' # Source data
#' M_source <- matrix(rnorm(300 * 50), nrow = 300)
#' colnames(M_source) <- paste0("M", 1:50)
#'
#' source_data <- data.frame(
#'   Y = rnorm(300),
#'   D = rnorm(300),
#'   M_source,
#'   X1 = rnorm(300)
#' )
#'
#' # Run SIS
#' result <- SIS(
#'   target_data = target_data,
#'   source_data = source_data,
#'   Y = "Y",
#'   D = "D",
#'   M = paste0("M", 1:50),
#'   X = "X1",
#'   transfer = TRUE,
#'   topN = 10
#' )
#'
#' result$M_ID_name_SIS
#' }
#'
#' @export


SIS <- function(
    target_data,
    source_data = NULL,
    Y,
    D,
    M,
    X,
    topN = NULL,
    transfer = FALSE,
    verbose = TRUE,
    ncore = 1,
    dblasso_method = FALSE
){

  if (verbose)
    message("Step 1: Sure Independence Screening ...  (",
            format(Sys.time(), "%X"), ")")

  if (ncore > 1) {
    doParallel::registerDoParallel(ncore)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }


  ## ---------- resolve columns ----------
  .resolve_cols <- function(data, cols, name) {
    if (is.null(cols)) return(NULL)
    if (is.numeric(cols)) return(colnames(data)[cols])
    if (is.character(cols)) return(cols)
    stop(paste(name, "must be names or indices"))
  }

  Y <- .resolve_cols(target_data, Y, "Y")
  D <- .resolve_cols(target_data, D, "D")
  M <- .resolve_cols(target_data, M, "M")
  X <- .resolve_cols(target_data, X, "X")

  stopifnot(length(Y) == 1, length(D) == 1, length(M) >= 1)

  p_m <- length(M)
  p_x <- ifelse(is.null(X), 0, length(X))
  n_t <- nrow(target_data)

  ## ---------- determine d_0 ----------
  if (is.null(topN)) {
    d_0 <- ceiling(2 * n_t / log(n_t))
  } else {
    d_0 <- topN
  }
  d_0 <- min(p_m, d_0)

  ## ---------- extract matrices ----------
  Y_t <- target_data[[Y]]
  D_t <- target_data[[D]]
  M_t <- as.matrix(target_data[, M, drop = FALSE])
  X_t <- if (p_x > 0) as.matrix(target_data[, X, drop = FALSE]) else NULL
  # X_t <- if (p_x > 0) .make_numeric_matrix(target_data, X) else NULL


  if (transfer) {
    Y_s <- source_data[[Y]]
    D_s <- source_data[[D]]
    M_s <- as.matrix(source_data[, M, drop = FALSE])
    X_s <- if (p_x > 0) as.matrix(source_data[, X, drop = FALSE]) else NULL
    # X_s <- if (p_x > 0) .make_numeric_matrix(source_data, X) else NULL
  }

  ## ============================================================
  ## 1. beta_SIS: mediator -> outcome
  ## ============================================================
  beta_SIS <- numeric(p_m)
  names(beta_SIS) <- M

  if (ncore == 1) {

    for (i in seq_len(p_m)) {

      x_t <- cbind(M_t[, i], D_t, X_t)
      colnames(x_t)[1] <- M[i]

      if (!transfer) {

        if (!dblasso_method) {
          # no trans, no db
          fit <- lasso(target = list(x = x_t, y = Y_t))
          beta_SIS[i] <- fit[M[i]]
          # message("mediator ",i, "\t beta",beta_SIS[i])
        } else {
          # no trans, db
          beta_SIS[i] <- dblasso(
            target = list(x = x_t, y = Y_t)
          )$dbcoef.hat[M[i]]
          # message("mediator ",i, "\t beta",beta_SIS[i])
        }

      } else {
        x_s <- cbind(M_s[, i], D_s, X_s)
        colnames(x_s)[1] <- M[i]

        if (dblasso_method) {
          # trans, db
          beta_SIS[i] <- dblasso(
            target = list(x = x_t, y = Y_t),
            source = list(x = x_s, y = Y_s),
            transfer = TRUE
          )$dbcoef.hat[M[i]]
          # message("mediator ",i, "\t beta",beta_SIS[i])
        }
        else{
          #trans, no db
          fit <- lasso(target = list(x = x_t, y = Y_t),
                       source = list(x = x_s, y = Y_s),
                       transfer = TRUE)
          beta_SIS[i] <- fit[M[i]]
          # message("mediator ",i, "\t beta",beta_SIS[i])
        }
      }
    }

  } else {

    beta_SIS <- foreach::foreach(
      i = seq_len(p_m),
      .combine = "c",
      .packages = c("stats", "glmnet","foreach")
      # .export = c("dblasso","lasso")
    ) %dopar% {

      x_t <- cbind(M_t[, i], D_t, X_t)
      colnames(x_t)[1] <- M[i]

      if (!transfer) {

        if (!dblasso_method) {
          fit <- lasso(target = list(x = x_t, y = Y_t))
          fit[M[i]]
        } else {
          dblasso(
            target = list(x = x_t, y = Y_t)
          )$dbcoef.hat[M[i]]
        }

      } else {

        x_s <- cbind(M_s[, i], D_s, X_s)
        colnames(x_s)[1] <- M[i]

        if (dblasso_method) {
          dblasso(
            target = list(x = x_t, y = Y_t),
            source = list(x = x_s, y = Y_s),
            transfer = TRUE
          )$dbcoef.hat[M[i]]
        }
        else{
          fit <- lasso(target = list(x = x_t, y = Y_t),
                       source = list(x = x_s, y = Y_s),
                       transfer = TRUE)
          fit[M[i]]
        }
      }
    }
    names(beta_SIS) <- M
  }


  ## ============================================================
  ## 2. alpha_SIS: exposure -> mediator
  ## ============================================================
  alpha_SIS <- numeric(p_m)
  names(alpha_SIS) <- M

  DX_t <- cbind(D_t, X_t)
  colnames(DX_t)[1] <- D

  if (transfer) DX_s <- cbind(D_s, X_s)

  if (ncore == 1) {

    for (i in seq_len(p_m)) {

      if (!transfer) {

        if (!dblasso_method) {
          # fit <- glm(M_t[, i] ~ ., data = as.data.frame(DX_t))
          # alpha_SIS[i] <- coef(fit)[D]
          fit <- lasso(target = list(x = DX_t, y = M_t[, i]))
          alpha_SIS[i] <- fit[D]
          # message("mediator ",i, "\t alpha",alpha_SIS[i])
        } else {
          alpha_SIS[i] <- dblasso(
            target = list(x = DX_t, y = M_t[, i])
          )$dbcoef.hat[D]
          # message("mediator ",i, "\t alpha",alpha_SIS[i])
        }

      } else {
        if(dblasso_method){
          alpha_SIS[i] <- dblasso(
            target = list(x = DX_t, y = M_t[, i]),
            source = list(x = DX_s, y = M_s[, i]),
            transfer = TRUE
          )$dbcoef.hat[D]
          # message("mediator ",i, "\t alpha",alpha_SIS[i])
        }
        else{
          fit <- lasso(target = list(x = DX_t, y = M_t[, i]),
                       source = list(x = DX_s, y = M_s[, i]),
                       transfer = TRUE)
          alpha_SIS[i] <- fit[D]
          # message("mediator ",i, "\t alpha",alpha_SIS[i])
        }
      }
    }

  } else {

    alpha_SIS <- foreach::foreach(
      i = seq_len(p_m),
      .combine = "c",
      .packages = c("stats", "glmnet","foreach")
      # .export = c("dblasso","lasso")
    ) %dopar% {

      if (!transfer) {

        if (!dblasso_method) {
          # fit <- glm(M_t[, i] ~ ., data = as.data.frame(DX_t))
          # coef(fit)[D]
          fit <- lasso(target = list(x = DX_t, y = M_t[, i]))
          fit[D]
        } else {
          dblasso(
            target = list(x = DX_t, y = M_t[, i])
          )$dbcoef.hat[D]
        }

      } else {

        if(dblasso_method){
          dblasso(
            target = list(x = DX_t, y = M_t[, i]),
            source = list(x = DX_s, y = M_s[, i]),
            transfer = TRUE
          )$dbcoef.hat[D]
        }
        else{
          fit <- lasso(target = list(x = DX_t, y = M_t[, i]),
                       source = list(x = DX_s, y = M_s[, i]),
                       transfer = TRUE)
          fit[D]
        }

      }
    }

    names(alpha_SIS) <- M
  }


  ## ============================================================
  ## 3. SIS selection
  ## ============================================================
  ab_SIS <- alpha_SIS * beta_SIS
  ID_SIS <- order(abs(ab_SIS), decreasing = TRUE)[seq_len(d_0)]
  M_ID_name_SIS <- M[ID_SIS]

  target_SIS <- target_data[, c(Y, D, M_ID_name_SIS, X), drop = FALSE]
  source_SIS <- if (transfer)
    source_data[, c(Y, D, M_ID_name_SIS, X), drop = FALSE]
  else NULL

  if (verbose)
    message("Top ", length(M_ID_name_SIS),
            " mediators selected: ",
            paste(M_ID_name_SIS, collapse = ", "),
            "  (", format(Sys.time(), "%X"), ")")

  res <- list(
    target_SIS = target_SIS,
    source_SIS = source_SIS,
    M_ID_name_SIS = M_ID_name_SIS
  )
  class(res) <- "SIS"
  return(res)
}

