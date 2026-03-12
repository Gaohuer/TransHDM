#' TransHDM: High-Dimensional Mediation Analysis with Transfer Learning
#'
#' The \code{TransHDM} function performs high-dimensional mediation analysis
#' under a transfer learning framework. It identifies and estimates indirect
#' (mediation) effects of a high-dimensional set of mediators between an
#' exposure and an outcome by integrating a target dataset and a source datasets.
#'
#' @param target_data A data frame containing the target dataset. All variables must be numeric.
#' @param source_data A list of data frames containing source datasets (optional, default: NULL).
#'   All variables must be numeric and have the same column names as target_data.
#' @param Y Character string specifying the outcome variable name.
#' @param D Character string specifying the exposure (treatment) variable name.
#' @param M Character vector specifying mediator variable names.
#' @param X Character vector specifying covariate variable names..
#' @param transfer A logical value (default: \code{FALSE}) indicating whether
#'   to enable transfer learning by incorporating information from
#'   \code{source_data}.
#' @param verbose A logical value (default: \code{TRUE}) controlling whether
#'   progress messages are printed to the console.
#' @param ncore An integer (default: 1) specifying the number of CPU cores to
#'   use for parallel computation when fitting mediator models.
#' @param topN An integer (default: \code{NULL}) specifying the number of
#'   mediators to retain after Sure Independence Screening (SIS).
#'   If \code{NULL}, the number is determined automatically based on
#'   the data dimensions.
#' @param dblasso_SIS A logical value (default: \code{FALSE}) indicating whether
#'   to apply a two-stage procedure combining SIS and debiased Lasso.
#'   When \code{TRUE}, mediators are first screened via SIS and then
#'   debiased Lasso is applied to the reduced set, which is recommended for
#'   ultra-high-dimensional settings.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{contributions}: A data frame of identified mediators containing:
#'     \itemize{
#'       \item \code{mediator}: Mediator name
#'       \item \code{alpha}: Estimated exposure mediator effect
#'       \item \code{alpha_pv}: P-value for the exposure mediator effect
#'       \item \code{beta}: Estimated mediator outcome effect
#'       \item \code{beta_pv}: P-value for the mediator outcome effect
#'       \item \code{alpha_beta}: Estimated indirect (mediation) effect
#'       \item \code{ab_pv}: Joint p-value for the indirect effect
#'       \item \code{pa}: Proportion of the total effect mediated
#'     }
#'   \item \code{effects}: A data frame summarizing the total indirect effect,
#'     direct effect, total effect, and proportion mediated.
#'   \item \code{IDE_est}: A numeric vector of indirect effect estimates for all
#'     specified mediators, with non-selected mediators set to zero.
#'   \item \code{DE_est}: The estimated direct effect of the exposure on the
#'     outcome.
#' }
#'
#' @importFrom stats setNames
#'
#' @references
#' Pan L, Liu Y, Huang C, Lin R, Yu Y, Qin G.
#' Transfer learning reveals the mediating mechanisms of cross-ethnic lipid metabolic
#' pathways in the association between APOE gene and Alzheimer's disease.
#' Brief Bioinform. 2025;26(5):bbaf460.
#' \doi{10.1093/bib/bbaf460}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' # Target data
#' target_data <- gen_simData_homo(n = 50, p_x = 3, p_m = 20, rho = 0.1)$data
#'
#' # Source data
#' source_data <- gen_simData_homo(n = 100, p_x = 3, p_m = 20, rho = 0.1, source = TRUE,
#' transferable = TRUE)$data
#'
#' # Run TransHDM
#' result <- TransHDM(
#'   target_data = target_data,
#'   source_data = source_data,
#'   Y = "Y",
#'   D = "D",
#'   M = paste0("M", 1:20),
#'   X = paste0("X", 1:3),
#'   transfer = TRUE,
#'   ncore = 1,
#'   topN = 10
#' )
#' summary(result)
#' }
#'
#' @export
#'

TransHDM <- function(
    target_data,
    source_data = NULL,
    Y,
    D,
    M,
    X,
    transfer = FALSE,
    verbose = TRUE,
    ncore = 1,
    topN = NULL,
    dblasso_SIS = FALSE
) {

  # ---- resolve variables ----
  Y <- .resolve_cols(target_data, Y, "Y")
  D <- .resolve_cols(target_data, D, "D")
  M <- .resolve_cols(target_data, M, "M")
  X <- .resolve_cols(target_data, X, "X")

  stopifnot(length(Y) == 1, length(D) == 1, length(M) >= 1)

  p_m <- length(M)

  # ---- Step 1: SIS ----
  # if (verbose)
  #   message("Step 1: Sure Independence Screening ...")

  SIS_result <- SIS(
    target_data = target_data,
    source_data = source_data,
    Y = Y,
    D = D,
    M = M,
    X = X,
    topN = topN,
    transfer = transfer,
    verbose = verbose,
    ncore = ncore,
    dblasso_method = dblasso_SIS
  )
  target_SIS <- SIS_result$target_SIS
  source_SIS <- SIS_result$source_SIS

  # ---- Step 2: De-biased Lasso ----
  ## mediators retained after SIS
  M_SIS <- intersect(M, colnames(target_SIS))
  p_m_SIS <- length(M_SIS)

  if (p_m_SIS == 0)
    stop("No mediators retained after SIS. Try increasing topN.")

  if (verbose && p_m_SIS < length(M))
    message("After SIS, ", p_m_SIS, " / ", length(M), " mediators are retained.")

  if (verbose)
    message("Step 2: De-biased Lasso Estimates ... (",
            format(Sys.time(), "%X"), ")")

  # ---------- 2.1 Outcome model (beta) ----------
  outcome_x <- setdiff(colnames(target_SIS), Y)

  target_outcome <- list(
    x = as.matrix(target_SIS[, outcome_x, drop = FALSE]),
    y = target_SIS[[Y]]
  )

  if (transfer) {
    source_outcome <- list(
      x = as.matrix(source_SIS[, outcome_x, drop = FALSE]),
      y = source_SIS[[Y]]
    )
    DBLASSO_fit <- dblasso(
      target = target_outcome,
      source = source_outcome,
      transfer = TRUE
    )
  } else {
    DBLASSO_fit <- dblasso(
      target = target_outcome,
      transfer = FALSE
    )
  }

  DE_est <- DBLASSO_fit$dbcoef.hat[D]

  beta_SIS_est <- setNames(rep(0, length(M)), M)
  P_beta_SIS   <- setNames(rep(1, length(M)), M)

  beta_SIS_est[M_SIS] <- DBLASSO_fit$dbcoef.hat[M_SIS]
  P_beta_SIS[M_SIS]   <- DBLASSO_fit$P.value[M_SIS]

  if (verbose)
    message("Estimation of mediator-outcome effects in the outcome model completed.")


  # ---------- 2.2 Mediator models (alpha) ----------
  alpha_SIS_est <- setNames(rep(0, length(M)), M)
  P_alpha_SIS   <- setNames(rep(1, length(M)), M)

  mediator_x <- c(D, X)
  mediator_x <- intersect(mediator_x, colnames(target_SIS))

  if (ncore == 1) {
    for (m in M_SIS) {
      target_m <- list(
        x = as.matrix(target_SIS[, mediator_x, drop = FALSE]),
        y = target_SIS[[m]]
      )
      if (transfer) {
        source_m <- list(
          x = as.matrix(source_SIS[, mediator_x, drop = FALSE]),
          y = source_SIS[[m]]
        )
        fit <- dblasso(
          target = target_m,
          source = source_m,
          transfer = TRUE
        )
      } else {
        fit <- dblasso(
          target = target_m,
          transfer = FALSE
        )
      }
      alpha_SIS_est[m] <- fit$dbcoef.hat[D]
      P_alpha_SIS[m]   <- fit$P.value[D]
    }

  } else {

    doParallel::registerDoParallel(ncore)

    alpha_res <- foreach::foreach(
      m = M_SIS,
      .combine = rbind,
      .packages = c("glmnet", "MASS", "foreach", "caret", "qvalue", "HDMT")
    ) %dopar% {

      target_m <- list(
        x = as.matrix(target_SIS[, mediator_x, drop = FALSE]),
        y = target_SIS[[m]]
      )

      if (transfer) {
        source_m <- list(
          x = as.matrix(source_SIS[, mediator_x, drop = FALSE]),
          y = source_SIS[[m]]
        )
        fit <- dblasso(
          target = target_m,
          source = source_m,
          transfer = TRUE
        )
      } else {
        fit <- dblasso(
          target = target_m,
          transfer = FALSE
        )
      }

      # c(p = fit$P.value[D], a = fit$dbcoef.hat[D])
      matrix(
        c(fit$P.value[D], fit$dbcoef.hat[D]),
        nrow = 1,
        dimnames = list(NULL, c("p", "a"))
      )
    }

    P_alpha_SIS[M_SIS]   <- alpha_res[, "p"]
    alpha_SIS_est[M_SIS] <- alpha_res[, "a"]

    closeAllConnections()
  }

  if (verbose)
    message("Estimation of exposure-mediator effects in the mediator model completed.")


  # ---- Step 3: Multiple testing ----
  if (verbose)
    message("Step 3: Multiple-testing procedure ... (",
            format(Sys.time(), "%X"), ")")

  PA <- cbind(P_alpha_SIS[M_SIS], P_beta_SIS[M_SIS])
  P_value <- apply(PA, 1, max)

  # input_pvalues <- PA + matrix(runif(length(PA), 0, 1e-10), ncol = 2)
  N0 <- dim(PA)[1] * dim(PA)[2]
  input_pvalues <- PA + matrix(runif(N0, 0, 10^{-10}), dim(PA)[1], 2)

  nullprop <- null_estimation(input_pvalues)

  fdrcut <- HDMT::fdr_est(
    nullprop$alpha00,
    nullprop$alpha01,
    nullprop$alpha10,
    nullprop$alpha1,
    nullprop$alpha2,
    input_pvalues,
    exact = 0
  )

  ID_fdr <- which(fdrcut <= 0.05)
  selected_M <- M_SIS[ID_fdr]

  alpha_hat_est <- alpha_SIS_est[selected_M]
  beta_hat_est  <- beta_SIS_est[selected_M]
  ab_pv <- P_value[ID_fdr]

  # ---- Effects ----
  alpha_beta_hat_est <- alpha_hat_est * beta_hat_est

  IDE_est <- setNames(rep(0, p_m), M)
  IDE_est[selected_M] <- alpha_beta_hat_est

  TE_est <- sum(IDE_est) + DE_est
  PE_est <- sum(IDE_est) / TE_est

  effects <- data.frame(
    effect = c("indirect", "direct", "total", "pe"),
    estimate = c(sum(IDE_est), DE_est, TE_est, PE_est)
  )

  contributions <- data.frame(
    mediator = selected_M,
    alpha = alpha_hat_est,
    alpha_pv = P_alpha_SIS[selected_M],
    beta = beta_hat_est,
    beta_pv = P_beta_SIS[selected_M],
    alpha_beta = alpha_beta_hat_est,
    ab_pv = ab_pv,
    pa = alpha_beta_hat_est / TE_est,
    row.names = NULL
  )

  if (verbose && length(selected_M) > 0)
    message("Identified mediator(s): ",
            paste(selected_M, collapse = ", "))

  fit <- list(
    contributions = contributions,
    effects = effects,
    IDE_est = IDE_est,
    DE_est = DE_est
  )
  class(fit) <- "TransHDM"
  return(fit)
}



.resolve_cols <- function(data, cols, name) {
  if (is.null(cols)) return(NULL)

  if (is.numeric(cols)) {
    if (any(cols < 1 | cols > ncol(data)))
      stop(sprintf("%s index out of range", name))
    return(colnames(data)[cols])
  }

  if (is.character(cols)) {
    if (!all(cols %in% colnames(data))) {
      stop(sprintf(
        "Some %s variables not found in data: %s",
        name,
        paste(setdiff(cols, colnames(data)), collapse = ", ")
      ))
    }
    return(cols)
  }

  stop(sprintf("%s must be column names or column indices", name))
}

