#' @method summary SIS
#' @export
summary.SIS <- function(object, ...) {

  if (!inherits(object, "SIS")) {
    stop("Object must be of class 'SIS'")
  }

  target_SIS <- object$target_SIS
  source_SIS <- object$source_SIS
  M_sel <- object$M_ID_name_SIS

  res <- list(
    selected_mediators = M_sel,
    n_selected = length(M_sel),
    target_dim = if (is.data.frame(target_SIS))
      c(n = nrow(target_SIS), p = ncol(target_SIS)) else NULL,
    source_dim = if (is.data.frame(source_SIS))
      c(n = nrow(source_SIS), p = ncol(source_SIS)) else NULL,
    transfer = is.data.frame(source_SIS)
  )

  class(res) <- c("summary.SIS", "SIS")
  res
}

#' @method print summary.SIS
#' @export
print.summary.SIS <- function(x, ...) {

  cat("\nSure Independence Screening (SIS) summary\n")
  cat("----------------------------------------\n")

  ## selected mediators
  cat("Selected mediators:\n")
  cat("  Number selected:", x$n_selected %||% 0, "\n")

  if (!is.null(x$selected_mediators) && length(x$selected_mediators) > 0) {
    cat("  Names:", paste(x$selected_mediators, collapse = ", "), "\n")
  } else {
    cat("  Names: (none)\n")
  }

  ## target info
  if (!is.null(x$target_dim)) {
    cat("\nTarget data after SIS:\n")
    cat("  Sample size (n):", x$target_dim["n"], "\n")
    cat("  Number of variables (p):", x$target_dim["p"], "\n")
  }

  ## source info
  if (!is.null(x$source_dim)) {
    cat("\nSource data after SIS:\n")
    cat("  Sample size (n):", x$source_dim["n"], "\n")
    cat("  Number of variables (p):", x$source_dim["p"], "\n")
    cat("\nTransfer learning: YES\n")
  } else {
    cat("\nTransfer learning: NO\n")
  }

  cat("----------------------------------------\n")

  invisible(x)
}


