#' Summary of TransHDM Mediation Analysis
#'
#' @param object An object of class \code{"TransHDM"}.
#' @param top Integer, maximum number of mediators to display.
#' @param digits Number of digits for rounding estimates.
#' @param ... Further arguments (unused).
#'
#' @return
#' An object of class \code{"summary.TransHDM"}.
#'
#' @method summary TransHDM
#' @export
summary.TransHDM <- function(object, top = 10, digits = 4, ...) {

  stopifnot(
    is.list(object),
    all(c("effects", "contributions") %in% names(object))
  )

  effects <- object$effects
  contrib <- object$contributions
  n_med <- nrow(contrib)

  ## ---- overall effects ----
  effects_show <- effects
  effects_show$estimate <- round(effects_show$estimate, digits)

  ## ---- mediator table ----
  if (n_med > 0) {
    show_n <- min(top, n_med)

    contrib_show <- contrib[
      order(abs(contrib$alpha_beta), decreasing = TRUE),
      ,
      drop = FALSE
    ][seq_len(show_n), ]

  } else {
    contrib_show <- NULL
  }

  res <- list(
    effects = effects_show,
    contributions = contrib_show,
    n_mediator = n_med,
    top = top,
    digits = digits
  )

  class(res) <- "summary.TransHDM"
  res
}


#' @method print summary.TransHDM
#' @export
print.summary.TransHDM <- function(x, ...) {

  p_cols <- c("alpha_pv", "beta_pv", "ab_pv")
  digits <- x$digits

  cat("====================================================\n")
  cat("        Summary of TransHDM Mediation Analysis        \n")
  cat("====================================================\n\n")

  ## ---- Overall effects ----
  cat("Overall Effects:\n")
  print(x$effects, row.names = FALSE)

  ## ---- Mediator summary ----
  cat("\nIdentified Mediators:\n")
  cat("  Number of selected mediators:", x$n_mediator, "\n")

  if (x$n_mediator == 0) {
    cat("  (No mediators selected at the given FDR level)\n")
    cat("\n====================================================\n")
    return(invisible(x))
  }

  cat("\nTop",
      min(x$top, x$n_mediator),
      "mediators by |alpha * beta|:\n\n")

  contrib_show <- x$contributions

  ## ---- formatting ----
  for (nm in names(contrib_show)) {
    if (is.numeric(contrib_show[[nm]])) {
      if (nm %in% p_cols) {
        contrib_show[[nm]] <- formatC(
          contrib_show[[nm]],
          format = "e",
          digits = digits
        )
      } else {
        contrib_show[[nm]] <- round(contrib_show[[nm]], digits)
      }
    }
  }

  print(contrib_show, row.names = FALSE)

  cat("\nNote:\n")
  cat("  alpha: exposure-mediator effect\n")
  cat("  beta : mediator-outcome effect\n")
  cat("  pa   : proportion of total effect explained\n")
  cat("  P-values are shown in scientific notation\n")

  cat("\n====================================================\n")

  invisible(x)
}


