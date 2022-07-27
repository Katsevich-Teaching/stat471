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
cross_validate_spline = function(x, y, nfolds, df_values){
  # a few checks of the inputs
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))

  # divide training data into folds
  n = length(x)
  train_data = tibble::tibble(x,y)
  folds = sample(rep(1:nfolds, length.out = n))
  train_data = train_data |> dplyr::mutate(fold = folds)

  # create a matrix for out-of-fold predictions
  num_df_values = length(df_values)
  out_of_fold_predictions =
    matrix(0, n, num_df_values) |>
    tibble::as_tibble() |>
    stats::setNames(paste0('y_hat_', df_values))

  # iterate over folds
  for(current_fold in 1:nfolds){
    # out-of-fold data will be used for training
    out_of_fold_data = train_data |> dplyr::filter(fold != current_fold)
    # in-fold data will be used for validation
    in_fold_data = train_data |> dplyr::filter(fold == current_fold)

    # iterate over df
    for(i in 1:num_df_values){
      df = df_values[i]

      # train on out-of-fold data
      formula = sprintf("y ~ splines::ns(x, df = %d)", df)
      spline_fit = stats::lm(formula = formula, data = out_of_fold_data)

      # predict on in-fold data
      out_of_fold_predictions[folds == current_fold, i] =
        stats::predict(spline_fit, newdata = in_fold_data)
    }
  }

  # add the out-of-fold predictions to the data frame
  results = train_data |> dplyr::bind_cols(out_of_fold_predictions)
  results

  # compute the CV estimate and standard error
  cv_table = results |>
    tidyr::pivot_longer(-c(x,y,fold),
                 names_to = "df",
                 names_prefix = "y_hat_",
                 names_transform = list(df = as.integer),
                 values_to = "yhat") |>
    dplyr::group_by(df, fold) |>
    dplyr::summarise(cv_fold = mean((yhat-y)^2)) |>  # CV estimates per fold
    dplyr::summarise(cv_mean = mean(cv_fold),
              cv_se = stats::sd(cv_fold)/sqrt(nfolds))

  df.1se = cv_table |>
    dplyr::filter(cv_mean-cv_se <= min(cv_mean)) |>
    dplyr::summarise(min(df)) |>
    dplyr::pull()

  df.min = cv_table |>
    dplyr::filter(cv_mean == min(cv_mean)) |>
    dplyr::summarise(min(df)) |>
    dplyr::pull()

  # plot the results, along with the previously computed validation error
  cv_plot = cv_table |>
    ggplot2::ggplot(ggplot2::aes(x = df, y = cv_mean, ymin = cv_mean-cv_se, ymax = cv_mean+cv_se)) +
    ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::geom_errorbar() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = min(cv_mean)), linetype = "dashed") +
    ggplot2::xlab("Degrees of freedom") + ggplot2::ylab("CV error") +
    ggplot2::theme_bw()

  # return CV table and plot
  return(list(cv_table = cv_table,
              cv_plot = cv_plot,
              df.1se = df.1se,
              df.min = df.min))
}
