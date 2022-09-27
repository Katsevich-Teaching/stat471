#' classification_metrics
#'
#' @param test_responses test responses
#' @param test_predictions test predictions
#' @param positive_class the name of the positive class
#' @param C_FP cost of making a false positive
#' @param C_FN cost of making a false negative
#'
#' @return tibble of classification metrics
#' @export
classification_metrics <- function(test_responses, test_predictions, positive_class = NULL, C_FP = NULL, C_FN = NULL){
  test_responses <- as.factor(test_responses)
  lvls <- levels(test_responses)
  if(length(lvls) != 2){
    stop("This function is designed only for binary classification problems.")
  }
  if(!is.null(positive_class)){
    if(lvls[2] != positive_class){
      lvls = rev(lvls)
    }
  }
  test_responses <- as.numeric(factor(test_responses, levels = lvls))-1
  test_predictions <- as.numeric(factor(test_predictions, levels = lvls))-1

  # compute misclassification error
  misclass_err <- mean(test_responses != test_predictions)

  # compute weighted misclassification error, if applicable
  if(!is.null(C_FP) & !is.null(C_FN)){
    weighted_misclass_err <- mean(C_FP*(test_responses == 0 & test_predictions == 1) +
                                    C_FN*(test_responses == 1 & test_predictions == 0))
  } else{
   weighted_misclass_err <- NA
  }

  # compute confusion matrix
  conf_matrix <- table(test_responses, test_predictions)

  # compute true positive rate
  TP <- conf_matrix["1", "1"]
  P <- sum(conf_matrix["1",])
  TPR <- TP/P

  # compute true negative rate
  TN <- conf_matrix["0", "0"]
  N <- sum(conf_matrix["0",])
  TNR <- TN/N

  # compute F-score
  F_score <- 1/(mean(c(1/TPR, 1/TNR)))

  # return
  tibble::tibble(
    misclass_err = misclass_err,
    w_misclass_err = weighted_misclass_err,
    TPR = TPR,
    TNR = TNR,
    `F` = F_score
  )
}
