#' Extract a tidy tibble of coefficients from model fit
#'
#' @param fit_object A fit object coming from `lm()`, `glm()`, `glmnet()`, or `cv.glmnet()`
#' @param s An optional parameter `s` for `glmnet()` and `cv.glmnet()`
#'
#' @return A tibble with columns `feature` and `coefficient`
#' @export
coef_tidy <- function(fit_object, s = NULL){
  # determine the type of fit that was used based on class attribute
  fit_class <- attr(fit_object, "class")
  if("cv.glmnet" %in% fit_class){
    fit_type <- "cv.glmnet"
    if(is.null(s)){
      s = "lambda.1se"
    }
  } else if("glmnet" %in% fit_class){
    fit_type <- "glmnet"
    if(is.null(s)){
      stop("The s argument must be specified for glmnet fits")
    } else{
      if(!is.numeric(s) | length(s) != 1){
        stop("The s argument must be numeric and length 1 for glmnet fits")
      }
    }
  } else if("lm" %in% fit_class | "glm" %in% fit_class){
    fit_type <- "unpenalized"
  } else{
    stop("Fit object must be obtained from either
           lm(), glm(), glmnet(), or cv.glmnet()")
  }

  if(fit_type %in% c("lm", "glm")){
    fit_object |>
      stats::coef() |>
      as.list() |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(tidyr::everything(),
                          names_to = "feature",
                          values_to = "coefficient")
  } else if(fit_type %in% c("glmnet", "cv.glmnet")){
    fit_object |>
      stats::coef(s = s) |>
      as.matrix() |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "feature") |>
      tibble::as_tibble() |>
      dplyr::rename(coefficient = "s1")
  }
}
