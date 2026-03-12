#' Fit LASSO Regression with Transfer Learning
#'
#' Fits a LASSO (Least Absolute Shrinkage and Selection Operator) regression model
#' under a transfer learning framework. Supports feature selection and coefficient
#' estimation by combining target data and source data.
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
#' @param transfer A logical value (default: FALSE) indicating whether to enable
#'   transfer learning mode (combine source data with target data).
#' @param lambda A string (default: 'lambda.1se') specifying the criterion for
#'   selecting regularization parameter:
#'   \itemize{
#'     \item{'lambda.min': Lambda value that gives minimum cross-validation error}
#'     \item{'lambda.1se': Largest lambda value within 1 standard error of the minimum error}
#'   }
#'
#' @return A numeric vector \code{coef} containing LASSO coefficient estimates (including intercept).
#'
#'
#' @examples
#' # Prepare target and source data
#' target <- list(x = matrix(rnorm(100 * 20), 100, 20), y = rnorm(100))
#' source <- list(x = matrix(rnorm(200 * 20), 200, 20), y = rnorm(200))
#'
#' # Non-transfer mode
#' coef_no_transfer <- lasso(target = target, transfer = FALSE, lambda = 'lambda.min')
#' print(coef_no_transfer)
#' summary(coef_no_transfer)
#'
#' # Transfer learning mode
#' coef_transfer <- lasso(target = target, source = source, transfer = TRUE, lambda = 'lambda.min')
#' print(coef_transfer)
#' summary(coef_transfer)
#'
#' @export
#'
lasso<-function(target, source = NULL, transfer=FALSE,lambda='lambda.1se'){

  orig_xnames <- colnames(target$x)
  if (is.null(orig_xnames)) {
    orig_xnames <- paste0("X", seq_len(ncol(target$x)))
  }
  coef.names <- c("(Intercept)", orig_xnames)

  x_t<-target$x
  y_t<-target$y
  x_s<-source$x
  y_s<-source$y

  orig_xnames <- colnames(target$x)
  if (is.null(orig_xnames)) {
    orig_xnames <- paste0("X", seq_len(ncol(target$x)))
  }

  if(!transfer){
    cv_glm<-glmnet::cv.glmnet(x=as.matrix(cbind(1,x_t)),
                      y=y_t,
                      alpha =1,
                      family="gaussian")

    coef<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm[[lambda]])],
            cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])][-1])
  }else{
    cv_glm<-glmnet::cv.glmnet(x=as.matrix(cbind(1,rbind(x_t,x_s))),
                      y=c(y_t,y_s),
                      alpha =1,
                      family="gaussian")
    w<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm[[lambda]])],
         cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])][-1])

    offset <- as.numeric(as.matrix(cbind(x_t)) %*% w[-1] + w[1])

    cv_glm<-glmnet::cv.glmnet(x=as.matrix(cbind(1,x_t)),
                      y=y_t,
                      alpha =1,
                      offset = offset,
                      family="gaussian")

    delta<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm[[lambda]])],
             cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])][-1])

    coef<-w+delta
  }
  names(coef) <- coef.names
  class(coef) <- "lasso"
  return(coef)
}

