#' K-Fold Cross-Validation Data Splitting
#'
#' Splits input data into k folds for cross-validation, generating training and test sets
#' for each fold. Particularly useful for mediator selection stability assessment
#' in high-dimensional mediation analysis.
#'
#' @param data A data frame or matrix containing the dataset to be split.
#'            Rows represent observations, columns represent variables.
#' @param kfold Integer (default: 3). Number of folds for cross-validation.
#'            Must be >= 2 and <= nrow(data).
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item{train_set: List of length kfold, each element is a training subset}
#'   \item{test_set: List of length kfold, each element is the corresponding test subset}
#' }
#'
#' @keywords internal

kfold_split<-function(data, kfold=3){
  if(kfold<2){
    stop("kfold should be at least 2")
  }
  n=dim(data)[1]
  split_indices <- caret::createFolds(1:n, k = kfold, list = TRUE)
  train_set<-list()
  test_set<-list()
  for(k in 1:kfold){
    train_set[[k]] <- data[-split_indices[[k]], ]
    test_set[[k]] <- data[split_indices[[k]], ]

  }
  return(list(train_set=train_set,test_set=test_set))
}
