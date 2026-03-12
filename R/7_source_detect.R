#' Detect Transferable Source Data via Cross-Validation
#'
#' Determines whether external source datasets can be effectively transferred to the target data
#' by comparing predictive performance using LASSO regression under a transfer learning framework.
#'
#' @param target_data A data frame containing the target dataset. All variables must be numeric.
#' @param source_data A list of data frames containing source datasets (optional, default: NULL).
#'   All variables must be numeric and have the same column names as target_data.
#' @param Y Character string specifying the outcome variable name.
#' @param D Character string specifying the exposure (treatment) variable name.
#' @param M Character vector specifying mediator variable names.
#' @param X Character vector specifying covariate variable names.
#' @param kfold Integer (default: 5). Number of folds for cross-validation.
#' @param C0 Numeric (default: 0.05). Threshold constant for determining transferability.
#'   Larger values make the criterion more lenient.
#' @param verbose Logical (default: TRUE). Whether to print progress messages.
#'
#' @return A list containing:
#' \itemize{
#'   \item{transfer.source.id: Indices of source datasets deemed transferable}
#'   \item{source.loss: Mean validation loss for each source dataset}
#'   \item{target.valid.loss: Mean validation loss using target-only model}
#'   \item{T_index: Difference between source loss and target-only loss for each source}
#'   \item{threshold: Calculated transferability threshold}
#'   \item{loss.cv: Full k-fold cross-validation loss matrix}
#' }
#'
#' @examples
#' ## Reproducible example
#' set.seed(123)
#'
#' # Generate synthetic target data
#' target_data <- data.frame(
#'   Y  = rnorm(200),
#'   D  = rnorm(200),
#'   M1 = rnorm(200),
#'   M2 = rnorm(200),
#'   X1 = rnorm(200)
#' )
#'
#' # Generate synthetic source data
#' source1 <- data.frame(
#'   Y  = rnorm(300),
#'   D  = rnorm(300),
#'   M1 = rnorm(300),
#'   M2 = rnorm(300),
#'   X1 = rnorm(300)
#' )
#'
#' source2 <- data.frame(
#'   Y  = rnorm(250),
#'   D  = rnorm(250),
#'   M1 = rnorm(250),
#'   M2 = rnorm(250),
#'   X1 = rnorm(250)
#' )
#'
#' # Run source detection
#' result <- source_detection(
#'   target_data = target_data,
#'   source_data = list(source1, source2),
#'   Y = "Y",
#'   D = "D",
#'   M = c("M1", "M2"),
#'   X = "X1",
#'   kfold = 5,
#'   C0 = 0.05,
#'   verbose = FALSE
#' )
#'
#' # Get Summary
#' summary(result)
#'
#' # Transferable source indices
#' result$transfer.source.id
#'
#' # Compare validation losses
#' data.frame(
#'   Source = c(paste0("Source", seq_along(result$source.loss)), "Target"),
#'   Loss   = c(result$source.loss, result$target.valid.loss)
#' )
#'
#' @export
source_detection <- function(target_data,
                             source_data = NULL,
                             Y,
                             D,
                             M,
                             X,
                             kfold = 5,
                             C0 = 0.05,
                             verbose = TRUE) {

  # Helper: K-fold split
  target_data_split <- kfold_split(target_data, kfold = kfold)

  # Ensure column names are characters
  Y <- if(is.numeric(Y)) colnames(target_data)[Y] else Y
  D <- if(is.numeric(D)) colnames(target_data)[D] else D
  M <- if(is.numeric(M)) colnames(target_data)[M] else M
  if(!is.null(X)) X <- if(is.numeric(X)) colnames(target_data)[X] else X

  # Internal loss function
  loss <- function(data, coef_beta) {
    n <- nrow(data)

    # Prepare design matrix dynamically
    design_cols <- c(D, M)
    if(!is.null(X)) design_cols <- c(design_cols, X)
    X_mat <- as.matrix(cbind(1, data[, design_cols, drop = FALSE]))
    mu <- X_mat %*% as.matrix(coef_beta)
    residuals <- data[[Y]] - mu
    Sigma <- sum(residuals^2) / n
    l.Y <- n/2 * log(2*pi) + n/2 * log(Sigma) + n/2
    return(l.Y)
  }

  # Cross-validation loop
  loss.cv <- t(sapply(1:kfold, function(k) {
    train <- target_data_split$train_set[[k]]
    test <- target_data_split$test_set[[k]]

    source.loss <- c()
    # if(transfer && !is.null(source_data)) {
    source.loss <- sapply(1:length(source_data), function(j) {
      train_list <- list(x = train[, c(D, M, X), drop = FALSE], y = train[[Y]])
      source_list <- list(x = source_data[[j]][, c(D, M, X), drop = FALSE], y = source_data[[j]][[Y]])
      coef_beta <- lasso(train_list, source_list, transfer = TRUE)
      loss(test, coef_beta)
    })
    # }

    # Target-only loss
    train_list <- list(x = train[, c(D, M, X), drop = FALSE], y = train[[Y]])
    coef_beta <- lasso(train_list, transfer = FALSE)
    target.loss <- loss(test, coef_beta)

    c(source.loss, target.loss)
  }))


  # Aggregate losses
  source.loss <- if(ncol(loss.cv) > 1) colMeans(loss.cv)[1:(ncol(loss.cv)-1)] else numeric(0)
  target.valid.loss <- colMeans(loss.cv)[ncol(loss.cv)]
  T_index <- source.loss - target.valid.loss
  threshold <- C0 * max(target.valid.loss, 0.01)
  transfer.source.id <- which(T_index <= threshold)

  if(verbose) {
    message("Transfer candidate sources: ", paste(transfer.source.id, collapse = ", "))
  }

  # Output as a structured list
  res <- list(
    transfer.source.id = transfer.source.id,
    source.loss = source.loss,
    target.valid.loss = target.valid.loss,
    T_index = T_index,
    threshold = threshold,
    loss.cv = loss.cv
  )
  class(res) <- "source_detection"
  res
}
