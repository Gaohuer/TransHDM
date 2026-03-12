#' Summary of Lasso Regression
#'
#' @param object A numeric vector of lasso coefficients with names.
#' @param ... Further arguments (currently not used).
#'
#' @return An object of class \code{"summary.lasso"}.
#'
#' @method summary lasso
#' @export
summary.lasso <- function(object, ...) {

  if (!is.numeric(object) || is.null(names(object))) {
    stop("Invalid lasso object")
  }

  coef <- object
  intercept <- coef["(Intercept)"]
  beta <- coef[names(coef) != "(Intercept)"]

  nonzero_idx <- which(beta != 0)
  nonzero_coef <- beta[nonzero_idx]

  res <- list(
    intercept = intercept,
    p = length(beta),
    n_nonzero = length(nonzero_coef),
    selected = names(nonzero_coef),
    l1 = sum(abs(beta)),
    l2 = sqrt(sum(beta^2)),
    max_abs = if (length(beta) > 0) max(abs(beta)) else 0
  )

  class(res) <- "summary.lasso"
  res
}


#' @method print summary.lasso
#' @export
print.summary.lasso <- function(x, ...) {

  cat("\nLasso regression summary\n")
  cat("-------------------------\n")

  ## intercept
  cat("Intercept:\n")
  cat(" ", formatC(x$intercept, digits = 6), "\n\n")

  ## model size
  cat("Model size:\n")
  cat("  Number of predictors:", x$p, "\n")
  cat("  Nonzero coefficients:", x$n_nonzero, "\n\n")

  ## selected variables
  cat("Selected variables:\n")
  if (x$n_nonzero > 0) {
    cat(" ", paste(x$selected, collapse = ", "), "\n\n")
  } else {
    cat("  (none)\n\n")
  }

  ## norms
  cat("Coefficient norms (excluding intercept):\n")
  cat("  L1 norm:", formatC(x$l1, digits = 6), "\n")
  cat("  L2 norm:", formatC(x$l2, digits = 6), "\n")
  cat("  Max |coef|:", formatC(x$max_abs, digits = 6), "\n")

  cat("-------------------------\n")

  invisible(x)
}

