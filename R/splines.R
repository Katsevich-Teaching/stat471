#' generate_spline_data
#'
#' Generate `resamples` datasets from Y = f(X)+epsilon for
#' X sampled on an equi-spaced grid of `n` samples on (0, 2*pi),
#' and noise standard deviation `sigma`.
#'
#'
#' @param f Function taking values on (0, 2*pi)
#' @param sigma Noise standard deviation
#' @param n Number of points to generate
#' @param resamples Number of resamples
#'
#' @return Data frame with all the datasets
#' @export
generate_spline_data <- function(f, sigma, n, resamples){
  tibble::tibble(
    x = rep(seq(0, 2 * pi, length.out = n), resamples),
    y = f(x) + stats::rnorm(n * resamples, sd = sigma),
    resample = rep(1:resamples, each = n)
  )
}

#' fit_spline_model
#'
#' @param data Data with columns x and y
#' @param df Degrees of freedom to use
#'
#' @return Fitted spline model
fit_spline_model = function(data, df) {
  stats::lm(y ~ splines::ns(x, df = df), data = data)
}

#' fit_spline_models
#'
#' Fit natural spline models with various degrees of freedom to one or more datasets.
#'
#' @param train_data_resamples Data frame with columns `x`, `y`, and optionally, `resample` if multiple datasets are given
#' @param df_values Vector of values of degrees of freedom for spline fits.
#'
#' @return Data frame with columns `x`, `y`, `resample`, `df`, `pred`
#' @export
fit_spline_models <- function(train_data_resamples, df_values){
  # add resample column if it is not present
  if(is.null(train_data_resamples$resample)){
    resample_column <- FALSE
    train_data_resamples <- train_data_resamples |>
      dplyr::mutate(resample = 1)
  } else{
    resample_column <- TRUE
  }

  # fit all the spline models
  fitted_spline_models <- train_data_resamples |>
    tidyr::crossing(df = df_values) |>
    dplyr::group_by(resample, df) |>
    tidyr::nest() |>
    dplyr::mutate(model = purrr::map2(data, df, fit_spline_model)) |>
    dplyr::mutate(fitted = purrr::map2(data, model, modelr::add_predictions)) |>
    dplyr::select(resample, df, fitted) |>
    tidyr::unnest(fitted) |>
    dplyr::ungroup() |>
    dplyr::select(x, y, resample, df, pred)

  # remove resample column if it was not originally present
  if(!resample_column){
    fitted_spline_models <- fitted_spline_models |> dplyr::select(-resample)
  }

  # return
  fitted_spline_models
}

#' compute_bias_variance_ETE
#'
#' @param fitted_spline_models Data frame outputed by `fit_spline_models()`
#' @param f The true function
#' @param sigma The noise standard deviation
#'
#' @return A data frame with columns `df`, `mean_sq_bias`, `mean_variance`,
#' and `expected_test_error`
#' @export
compute_bias_variance_ETE <- function(fitted_spline_models, f, sigma){
  fitted_spline_models |>
    dplyr::mutate(true_fit = f(x)) |>
    dplyr::group_by(df, x) |>
    dplyr::summarise(bias = mean(pred - true_fit),
              variance = stats::var(pred)) |>
    dplyr::summarise(mean_sq_bias = mean(bias^2),
              mean_variance = mean(variance)) |>
    dplyr::mutate(expected_test_error = mean_sq_bias + mean_variance + sigma^2)
}

#' bias_variance_ETE_plot
#'
#' Create a bias/variance/ETE plot based on output of `compute_bias_variance_ETE`
#'
#' @param bias_variance_ETE The output of `compute_bias_variance_ETE`
#'
#' @return `ggplot` object containing bias/variance/ETE plot
#' @export
bias_variance_ETE_plot <- function(bias_variance_ETE){
  bias_variance_ETE |>
    tidyr::pivot_longer(-df, names_to = "metric", values_to = "error") |>
    ggplot2::ggplot(ggplot2::aes(x = df, y = error, colour = metric)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}

#' cross_validate_spline
#'
#' Run cross-validation to select the degrees of freedom
# for a natural spline fit.
#'
#' @param x vector of x coordinates in training data
#' @param y vector of y coordinates in training data
#' @param nfolds number of folds for cross-validation
#' @param df_values vector of values of degrees of freedom to try
#'
#' @return An object containing the CV table, CV plot, df.1se and df.min
#' @export
cross_validate_spline <- function(x, y, nfolds, df_values) {
  # a few checks of the inputs
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))

  # divide training data into folds
  n <- length(x)
  train_data <- tibble::tibble(x, y)
  folds <- sample(rep(1:nfolds, length.out = n))
  train_data <- train_data |> dplyr::mutate(fold = folds)

  # create a matrix for out-of-fold predictions
  num_df_values <- length(df_values)
  out_of_fold_predictions <-
    matrix(0, n, num_df_values) |>
    `colnames<-`(paste0("y_hat_", df_values)) |>
    tibble::as_tibble(.name_repair = 'unique')

  # iterate over folds
  for (current_fold in 1:nfolds) {
    # out-of-fold data will be used for training
    out_of_fold_data <- train_data |> dplyr::filter(fold != current_fold)
    # in-fold data will be used for validation
    in_fold_data <- train_data |> dplyr::filter(fold == current_fold)

    # iterate over df
    for (i in 1:num_df_values) {
      df <- df_values[i]

      # train on out-of-fold data
      formula <- sprintf("y ~ splines::ns(x, df = %d)", df)
      spline_fit <- stats::lm(formula = formula, data = out_of_fold_data)

      # predict on in-fold data
      out_of_fold_predictions[folds == current_fold, i] <-
        stats::predict(spline_fit, newdata = in_fold_data)
    }
  }

  # add the out-of-fold predictions to the data frame
  results <- train_data |> dplyr::bind_cols(out_of_fold_predictions)

  # compute the CV estimate and standard error
  cv_table <- results |>
    tidyr::pivot_longer(-c(x, y, fold),
      names_to = "df",
      names_prefix = "y_hat_",
      names_transform = list(df = as.integer),
      values_to = "yhat"
    ) |>
    dplyr::group_by(df, fold) |>
    dplyr::summarise(cv_fold = mean((yhat - y)^2)) |> # CV estimates per fold
    dplyr::summarise(
      cv_mean = mean(cv_fold),
      cv_se = stats::sd(cv_fold) / sqrt(nfolds)
    )

  df.min <- cv_table |>
    dplyr::filter(cv_mean == min(cv_mean)) |>
    dplyr::summarise(min(df)) |>
    dplyr::pull()

  df.1se_cv_threshold <- cv_table |>
    dplyr::filter(df == df.min) |>
    dplyr::summarise(cv_max = cv_mean + cv_se) |>
    dplyr::pull()

  df.1se <- cv_table |>
    dplyr::filter(cv_mean <= df.1se_cv_threshold) |>
    dplyr::summarise(min(df)) |>
    dplyr::pull()

  # plot the results, along with the previously computed validation error
  cv_plot <- cv_table |>
    ggplot2::ggplot(ggplot2::aes(x = df, y = cv_mean, ymin = cv_mean - cv_se, ymax = cv_mean + cv_se)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = min(cv_mean)), linetype = "dashed") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = df.1se_cv_threshold), linetype = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = df.1se), linetype = "dotted") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = df.min), linetype = "dotted") +
    ggplot2::xlab("Degrees of freedom") +
    ggplot2::ylab("CV error") +
    ggplot2::theme_bw()

  # return CV table and plot
  return(list(
    cv_table = cv_table,
    cv_plot = cv_plot,
    df.1se = df.1se,
    df.min = df.min
  ))
}
