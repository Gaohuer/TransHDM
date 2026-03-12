#' Summary of Debiased Lasso Inference
#'
#' @param object An object of class \code{"dblasso"}.
#' @param ... Further arguments (currently not used).
#'
#' @return An object of class \code{"summary.dblasso"}.
#'
#' @method summary dblasso
#' @export
summary.dblasso <- function(object, ...) {

  if (!inherits(object, "dblasso")) {
    stop("Object must be of class 'dblasso'")
  }

  est <- object$dbcoef.hat
  se  <- object$se.est
  z   <- est / se
  p   <- object$P.value
  CI  <- object$CI

  tab <- data.frame(
    Estimate  = est,
    Std.Error = se,
    Z.value   = z,
    Pr...z..  = p,
    CI.Lower  = CI$lb,
    CI.Upper  = CI$ub,
    row.names = names(est)
  )

  res <- list(
    coef_table = tab,
    n_coef = length(est)
  )

  class(res) <- "summary.dblasso"
  res
}


#' @method print summary.dblasso
#' @export
print.summary.dblasso <- function(x, ...) {

  cat("Debiased Lasso Inference Summary\n")
  cat(strrep("-", 32), "\n", sep = "")

  stats::printCoefmat(
    x$coef_table,
    digits = 4,
    signif.stars = TRUE
  )

  invisible(x)
}

