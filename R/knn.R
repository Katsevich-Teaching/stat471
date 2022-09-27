#' Weighted K-Nearest Neighbors
#'
#' Carry out weighted KNN, e.g. to handle imbalanced classes.
#'
#' @param formula Formula like one would use for `lm()`
#' @param training_data The training data, in the form of a data frame
#' @param test_data The test data, in the form of a data frame
#' @param k The number of nearest neighbors
#' @param weights A vector of weights
#' @param algorithm Nearest neighbor search algorithm.
#'
#' @return An object with fields `predictions` and `probabilities`. `predictions` is a vector
#' of the same length as the number of test observations indicating the predicted classes, and
#' `probabilities` is a tibble with columns `observation`, `class`, and `probability`.
#' @export
#'
#' @examples
#' formula = y ~ x
#' training_data <- data.frame(x = 1:6, y = factor(c(1, 1, 1, 2, 2, 3)))
#' test_data <- data.frame(x = 1:6)
#' weights <- exp(1:6)
#' k <- 3
#' knn(formula = formula,
#' training_data = training_data,
#' test_data = test_data,
#' k = k,
#' weights = weights)
knn <- function(formula,
                training_data,
                test_data,
                k = 1,
                weights = NULL,
                algorithm = c("kd_tree", "cover_tree", "brute")){

  # if weights is null, set it to all ones
  if(is.null(weights)){
    weights <- rep(1, nrow(training_data))
  }

  # extract training class labels
  cl <- formula |>
    stats::model.frame(data = training_data) |>
    stats::model.response() |>
    as.factor()

  # extract training features
  train <- stats::model.matrix(object = formula, data = training_data)

  # extract test features
  Terms <- formula |>
    stats::terms() |>
    stats::delete.response()
  m <- stats::model.frame(Terms, test_data, xlev = formula$xlevels)
  test <- stats::model.matrix(Terms, m, contrasts.arg = formula$contrasts)

  # apply knn function from FNN package to get the indices of the neighbors
  nn.index <- FNN::knn(train, test, cl, k, algorithm = algorithm) |>
    attr("nn.index")

  # define function to get predicted probabilities for each observation
  pred_prob <- function(nn.indices) {
    sapply(
      levels(cl),
      function(level){
        sum(weights[nn.indices]*(cl[nn.indices] == level))/sum(weights[nn.indices])
      }
    )
  }

  # apply function to each row of nn.index and then tidy
  probabilities <- nn.index |>
    apply(1, pred_prob) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "observation") |>
    tidyr::pivot_longer(-observation,
                        names_to = "class",
                        values_to = "probability") |>
    dplyr::mutate(observation = as.integer(observation),
           class = as.factor(class))

  # get predictions as classes with maximum probability
  predictions <- probabilities |>
    dplyr::group_by(observation) |>
    dplyr::filter(probability == max(probability)) |>
    dplyr::filter(dplyr::row_number() == 1) |> # in case there are ties, choose first one
    dplyr::ungroup() |>
    dplyr::arrange(observation) |>
    dplyr::pull(class)

  # return predictions and probabilities
  list(predictions = predictions, probabilities = probabilities)
}

#' Weighted K-Nearest Neighbors visualization
#'
#' Carry out weighted KNN and plot the results.
#'
#' @param formula Formula like one would use for `lm()`
#' @param training_data The training data, in the form of a data frame
#' @param test_data The test data, in the form of a data frame
#' @param k The number of nearest neighbors
#' @param weights A vector of weights
#' @param algorithm Nearest neighbor search algorithm.
#'
#' @return A `ggplot` object
#' @export
visualize_knn <- function(formula,
                          training_data,
                          test_data,
                          k = 1,
                          weights = NULL,
                          algorithm = c("kd_tree", "cover_tree", "brute")){
  # get values of the feature for the training dataset
  feature_train <- formula |>
    stats::terms() |>
    stats::delete.response() |>
    stats::get_all_vars(data = training_data)

  # get values of the response for the training dataset
  response_train <- formula |>
    stats::model.frame(data = training_data) |>
    stats::model.response()

  # get values of the feature for the test dataset
  feature_test <- formula |>
    stats::terms() |>
    stats::delete.response() |>
    stats::get_all_vars(data = test_data)

  # get values of the response for the test dataset
  response_test <- formula |>
    stats::model.frame(data = test_data) |>
    stats::model.response()

  # make sure there is exactly one feature
  if(!(ncol(feature_train) == 1 & ncol(feature_test) == 1)){
    stop("This function is designed for the case when there is exactly one feature.")
  } else{
    feature_train <- feature_train |> dplyr::pull()
    feature_test <- feature_test |> dplyr::pull()
  }

  # make sure response is binary
  if(length(unique(response_train)) != 2){
    stop("This function is designed for the case when the response is binary.")
  }

  # run KNN to make predictions on training and test datasets
  knn_results_train <- knn(formula, training_data, training_data, k, weights, algorithm)
  knn_results_test <- knn(formula, training_data, test_data, k, weights, algorithm)

  # get names of negative and positive classes
  neg_class_name <- levels(response_train)[1]
  pos_class_name <- levels(response_train)[2]

  # get names of feature and response
  variables <- formula |>
    stats::terms() |>
    attr("variables")
  feature_name <- variables[3] |> as.character()
  response_name <-variables[2] |> as.character()

  # collate information and plot
  rbind(
    tibble::tibble(
      feature = feature_train,
      response = as.numeric(response_train) - 1,
      probability = knn_results_train$probabilities |>
        dplyr::filter(class == pos_class_name) |>
        dplyr::pull(probability),
      prediction = knn_results_train$predictions,
      set = "train"
    ),
    tibble::tibble(
      feature = feature_test,
      response = as.numeric(response_test) - 1,
      probability = knn_results_test$probabilities |>
        dplyr::filter(class == pos_class_name) |>
        dplyr::pull(probability),
      prediction = knn_results_test$predictions,
      set = "test"
    )
  ) |>
    dplyr::mutate(set = factor(set, levels = c("train", "test"))) |>
    ggplot2::ggplot(ggplot2::aes(x = feature)) +
    ggplot2::geom_jitter(ggplot2::aes(y = response, colour = prediction),
                width = 0, height = 0.1) +
    ggplot2::geom_line(ggplot2::aes(y = probability)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") +
    ggplot2::facet_grid(~set) +
    ggplot2::scale_y_continuous(breaks = c(0,1),
                       labels = c(sprintf("%s (0)", neg_class_name),
                                  sprintf("%s (1)", pos_class_name))) +
    ggplot2::scale_colour_discrete(guide = ggplot2::guide_legend(title = "Prediction")) +
    ggplot2::labs(x = feature_name, y = response_name) +
    ggplot2::theme(legend.position = "bottom")
}
