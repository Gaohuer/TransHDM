#' Summary of Source Detection Results
#'
#' @param object An object of class \code{"source_detection"}.
#' @param ... Further arguments (unused).
#'
#' @return
#' An object of class \code{"summary.source_detection"}.
#'
#' @method summary source_detection
#' @export
summary.source_detection <- function(object, ...) {

  n_source <- length(object$source.loss)
  n_transfer <- length(object$transfer.source.id)

  if (n_source > 0) {
    source_table <- data.frame(
      Source = seq_len(n_source),
      SourceLoss = round(object$source.loss, 4),
      T_index = round(object$T_index, 4),
      Transferable = ifelse(
        seq_len(n_source) %in% object$transfer.source.id,
        "YES", "NO"
      )
    )
  } else {
    source_table <- NULL
  }

  res <- list(
    n_source = n_source,
    n_transfer = n_transfer,
    target.valid.loss = round(object$target.valid.loss, 4),
    threshold = round(object$threshold, 4),
    source_table = source_table
  )

  class(res) <- "summary.source_detection"
  res
}


#' @method print summary.source_detection
#' @export
print.summary.source_detection <- function(x, ...) {

  cat("Source detection summary\n")
  cat(strrep("-", 24), "\n", sep = "")
  cat("Number of sources:", x$n_source, "\n")
  cat("Number of transferable sources:", x$n_transfer, "\n\n")

  cat("Target validation loss:\n")
  cat("  mean =", x$target.valid.loss, "\n\n")

  cat("Threshold:\n")
  cat(" ", x$threshold, "\n\n")

  if (!is.null(x$source_table)) {
    cat("Source-wise comparison:\n")
    print(x$source_table, row.names = FALSE)
  } else {
    cat("No source data provided.\n")
  }

  invisible(x)
}


