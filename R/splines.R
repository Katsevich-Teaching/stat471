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
  stats::lm(y ~ splines::ns(x, df = df-1), data = data)
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
    dplyr::mutate(metric = factor(
      metric,
      levels = c("mean_sq_bias", "mean_variance", "expected_test_error"),
      labels = c("Mean squared bias", "Mean variance", "Expected test error"))) |>
    ggplot2::ggplot(ggplot2::aes(x = df, y = error, colour = metric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Degrees of freedom",
                  y = "Error") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}

#' Create bias-variance tradeoff plot
#'
#' @param f The true function
#' @param sigma The noise level
#' @param n The training sample size
#' @param df_values The degrees of freedom values to use
#' @param resamples The number of resamples to use
#'
#' @return The bias-variance tradeoff plot
#' @export
bias_variance_tradeoff <- function(f, sigma, n, df_values = 2:15, resamples = 100){
  # set seed for reproducibility
  set.seed(1)

  # create bias variance trade-off plot
  generate_spline_data(f, sigma, n, resamples) |>  # generate training datasets
    fit_spline_models(df_values = df_values) |> # fit natural splines to training data
    compute_bias_variance_ETE(f,sigma) |> # compute bias, variance, and ETE
    bias_variance_ETE_plot() |>  # plot the bias, variance, and ETE
    plot()
}

#' Create bias-variance tradeoff plot
#'
#' @param f The true function
#' @param sigma The noise level
#' @param n The training sample size
#' @param df The degrees of freedom to use
#' @param resamples The number of resamples to use
#'
#' @return The plot showing the true function and the spline fits
#' @export
spline_simulation_single <- function(f, sigma, n, df, resamples = 100){
  RESAMPLES_TO_PLOT <- 10
  # set seed for reproducibility
  set.seed(1)

  train_data_resamples <- generate_spline_data(f, sigma, n, resamples)
  training_results <- fit_spline_models(train_data_resamples, df_values = df)
  bias_variance_ETE <- training_results |> compute_bias_variance_ETE(f,sigma)
  training_results_summary <- training_results |>
    dplyr::mutate(true_fit = f(x)) |>
    dplyr::group_by(x) |>
    dplyr::summarise(mean_pred = mean(pred),
              bias = mean(pred - true_fit),
              variance = var(pred))
  training_results |>
    dplyr::filter(resample <= RESAMPLES_TO_PLOT) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = pred)) +
    ggplot2::geom_line(ggplot2::aes(colour = "Training fits",
                  alpha = "Training fits",
                  linetype = "Training fits",
                  group = resample),
              linewidth = 1) +
    ggplot2::stat_function(fun = f,
                           ggplot2::aes(colour = "True function",
                      alpha = "True function",
                      linetype = "True function"),
                  linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = mean_pred,
                  colour = "Mean training fit",
                  alpha = "Mean training fit",
                  linetype = "Mean training fit"),
                  linewidth = 1,
              data = training_results_summary) +
    ggplot2::scale_colour_manual(name = "legend", values = c("blue", "red", "red"),
                        breaks = c("True function",
                                   "Training fits",
                                   "Mean training fit")) +
    ggplot2::scale_linetype_manual(name = "legend", values = c("solid", "solid", "dashed"),
                          breaks = c("True function",
                                     "Training fits",
                                     "Mean training fit")) +
    ggplot2::scale_alpha_manual(name = "legend", values = c(1,0.25,1),
                       breaks = c("True function",
                                  "Training fits",
                                  "Mean training fit")) +
    ggplot2::labs(x = "x",
                  y = "y",
                  title = sprintf("Mean sq bias = %0.2f; Mean variance = %0.2f; ETE = %0.2f",
                                  bias_variance_ETE$mean_sq_bias,
                                  bias_variance_ETE$mean_variance,
                                  bias_variance_ETE$expected_test_error)) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
}

#' Create bias-variance tradeoff plot
#'
#' @param f The true function
#' @param sigma The noise level
#' @param n The training sample size
#' @param df The degrees of freedom to use
#' @param resamples The number of resamples to use
#'
#' @return The plot showing the true function and the spline fits
#' @export
spline_simulation <- function(f, sigma, n, df, resamples = 100){
  RESAMPLES_TO_PLOT <- 10
  # set seed for reproducibility
  set.seed(1)

  varying_params <- list(f = f, sigma = sigma, n = n, df = df)
  varying_param_name <- names(which(sapply(varying_params, function(x) length(x) > 1)))
  if(length(varying_param_name) > 1){
    stop("Only one of the four parameters f, sigma, n, df can be varied.")
  }
  if(length(varying_param_name) == 0){
    return(spline_simulation_single(f, sigma, n, df, resamples))
  }

  varying_param_values <- varying_params[[varying_param_name]]

  results_list <- vector("list", length(varying_param_values))

  all_true_function_data <- list()

  for (i in seq_along(varying_param_values)) {
    varying_params_i <- varying_params
    varying_params_i[[varying_param_name]] <- varying_param_values[[i]]
    f <- varying_params_i[["f"]]
    sigma <- varying_params_i[["sigma"]]
    n <- varying_params_i[["n"]]
    df <- varying_params_i[["df"]]

    train_data_resamples <- generate_spline_data(f, sigma, n, resamples)
    training_results <- fit_spline_models(train_data_resamples, df_values = df)
    bias_variance_ETE <- training_results |> compute_bias_variance_ETE(f,sigma)
    training_results_summary <- training_results |>
      dplyr::mutate(true_fit = f(x)) |>
      dplyr::group_by(x) |>
      dplyr::summarise(mean_pred = mean(pred),
                       bias = mean(pred - true_fit),
                       variance = var(pred))

    if (varying_param_name == "f") {
      x_grid <- seq(0, 2 * pi, length.out = 1000)
      true_function_data <- data.frame(x = x_grid, y = f(x_grid)) |>
        dplyr::mutate(varying_parameter = paste0("f = ", deparse(body(varying_param_values[[i]]))))
      all_true_function_data <- rbind(all_true_function_data, true_function_data)
    }

    # store data, summary, and metrics in a list
    results_list[[i]] <- list(
      training_results = training_results,
      training_results_summary = training_results_summary,
      bias_variance_ETE = bias_variance_ETE
    )
  }

  all_training_results <- dplyr::bind_rows(lapply(results_list, `[[`, 'training_results'), .id = "idx")
  all_training_results_summary <- dplyr::bind_rows(lapply(results_list, `[[`, 'training_results_summary'), .id = "idx")
  all_bias_variance_ETE <- dplyr::bind_rows(lapply(results_list, `[[`, 'bias_variance_ETE'), .id = "idx")

  if(varying_param_name == "f"){
    varying_param_values <- sapply(varying_param_values, function(fun)(deparse(body(fun))))
  }
  all_training_results <- all_training_results |>
    dplyr::mutate(varying_parameter = paste(varying_param_name, varying_param_values[as.integer(idx)], sep = " = ")) |>
    dplyr::mutate(varying_parameter = factor(varying_parameter, levels = paste(varying_param_name, varying_param_values, sep = " = ")))
  all_training_results_summary <- all_training_results_summary |>
    dplyr::mutate(varying_parameter = paste(varying_param_name, varying_param_values[as.integer(idx)], sep = " = ")) |>
    dplyr::mutate(varying_parameter = factor(varying_parameter, levels = paste(varying_param_name, varying_param_values, sep = " = ")))
  all_bias_variance_ETE <- all_bias_variance_ETE |>
    dplyr::mutate(varying_parameter = paste(varying_param_name, varying_param_values[as.integer(idx)], sep = " = ")) |>
    dplyr::mutate(varying_parameter = factor(varying_parameter, levels = paste(varying_param_name, varying_param_values, sep = " = ")))

  plot <- all_training_results |>
    dplyr::filter(resample <= RESAMPLES_TO_PLOT) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = pred)) +
    ggplot2::geom_line(ggplot2::aes(colour = "Training fits",
                                    alpha = "Training fits",
                                    linetype = "Training fits",
                                    group = resample),
                       linewidth = 1) +
    # ggplot2::stat_function(fun = f,
    #                        ggplot2::aes(colour = "True function",
    #                                     alpha = "True function",
    #                                     linetype = "True function"),
    #                        linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = mean_pred,
                                    colour = "Mean training fit",
                                    alpha = "Mean training fit",
                                    linetype = "Mean training fit"),
                       linewidth = 1,
                       data = all_training_results_summary) +
    ggplot2::scale_colour_manual(name = "legend", values = c("blue", "red", "red"),
                                 breaks = c("True function",
                                            "Training fits",
                                            "Mean training fit")) +
    ggplot2::scale_linetype_manual(name = "legend", values = c("solid", "solid", "dashed"),
                                   breaks = c("True function",
                                              "Training fits",
                                              "Mean training fit")) +
    ggplot2::scale_alpha_manual(name = "legend", values = c(1,0.25,1),
                                breaks = c("True function",
                                           "Training fits",
                                           "Mean training fit")) +
    # Add statistics as text annotation
    ggplot2::geom_text(data = all_bias_variance_ETE,
                       ggplot2::aes(label = sprintf("Mean sq bias = %0.2f; Mean variance = %0.2f; ETE = %0.2f",
                                  mean_sq_bias, mean_variance, expected_test_error),
                  x = -Inf, y = -Inf),
              hjust = -0.1, vjust = -0.5, size = 2.5) +
    # Add faceting based on the varying parameter
    ggplot2::facet_wrap(~varying_parameter, ncol = 1) +
    ggplot2::labs(x = "x",
                  y = "y") +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())

  if (varying_param_name == "f") {
    plot <- plot +
      ggplot2::geom_line(data = all_true_function_data, ggplot2::aes(x = x, y = y, colour = "True function",
                                                   linetype = "True function",
                                                   alpha = "True function"), linewidth = 1)
  } else {
    plot <- plot +
      ggplot2::stat_function(fun = f, ggplot2::aes(colour = "True function",
                                 linetype = "True function",
                                 alpha = "True function"), linewidth = 1)
  }
  plot
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
  # set seed for reproducibility
  set.seed(6)

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
      formula <- sprintf("y ~ splines::ns(x, df = %d)", df - 1)
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
    dplyr::summarise(cv_fold = mean((yhat - y)^2), .by = c(df, fold)) |> # CV estimates per fold
    dplyr::summarise(
      cv_mean = mean(cv_fold),
      cv_se = stats::sd(cv_fold) / sqrt(nfolds),
      .by = df
    )

  df.min <- cv_table |>
    dplyr::filter(cv_mean == min(cv_mean)) |>
    dplyr::summarise(min(df)) |>
    dplyr::pull()

  df.1se.thresh <- cv_table |>
    dplyr::filter(df == df.min) |>
    dplyr::summarise(cv_mean + cv_se) |>
    dplyr::pull()

  df.1se <- cv_table |>
    dplyr::filter(cv_mean <= df.1se.thresh) |>
    dplyr::summarise(min(df)) |>
    dplyr::pull()

  # plot the results, along with the previously computed validation error
  cv_plot <- cv_table |>
    ggplot2::ggplot(ggplot2::aes(x = df, y = cv_mean, ymin = cv_mean - cv_se, ymax = cv_mean + cv_se)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = min(cv_mean)), linetype = "dashed") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = df.1se.thresh), linetype = "dashed") +
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
