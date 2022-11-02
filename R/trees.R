#' Extract CV information from decision tree fit
#'
#' @param tree_fit Tree fit object returned by rpart
#'
#' @return An object containing the CV table, CV plot, nleaves.1se and nleaves.min
#' @export
cv_tree <- function(tree_fit){
  # get the CV table
  cv_table <- tree_fit$cptable |>
    as.data.frame() |>
    dplyr::rename(cv_mean = xerror,
           cv_se = xstd) |>
    dplyr::select(CP, nsplit, cv_mean, cv_se) |>
    dplyr::mutate(nleaves = nsplit + 1, .before = cv_mean)

  nleaves.min <- cv_table |>
    dplyr::filter(cv_mean == min(cv_mean)) |>
    dplyr::summarise(min(nleaves)) |>
    dplyr::pull()

  nleaves.1se_cv_threshold <- cv_table |>
    dplyr::filter(nleaves == nleaves.min) |>
    dplyr::summarise(cv_max = cv_mean + cv_se) |>
    dplyr::pull()

  nleaves.1se <- cv_table |>
    dplyr::filter(cv_mean <= nleaves.1se_cv_threshold) |>
    dplyr::summarise(min(nleaves)) |>
    dplyr::pull()

  CP.1se <- cv_table |>
    dplyr::filter(nleaves == nleaves.1se) |>
    dplyr::pull(CP)

  CP.min <- cv_table |>
    dplyr::filter(nleaves == nleaves.min) |>
    dplyr::pull(CP)

  # plot the results, along with the previously computed validation error
  cv_plot <- cv_table |>
    ggplot2::ggplot(ggplot2::aes(x = nleaves, y = cv_mean, ymin = cv_mean - cv_se, ymax = cv_mean + cv_se)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = min(cv_mean)), linetype = "dashed") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = nleaves.1se_cv_threshold), linetype = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = nleaves.1se), linetype = "dotted") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = nleaves.min), linetype = "dotted") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::xlab("Number of terminal nodes") +
    ggplot2::ylab("CV error") +
    ggplot2::theme_bw()

  # return CV table and plot
  return(list(
    cv_plot = cv_plot,
    nleaves.1se = nleaves.1se,
    nleaves.min = nleaves.min,
    CP.1se = CP.1se,
    CP.min = CP.min
  ))
}
