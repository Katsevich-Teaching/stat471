#' Computed mixed correlation matrix
#'
#' Like a regular correlation matrix, but allowing arbitrary column types.
#'
#' @param data A data frame or tibble.
#'
#' @return Correlation matrix
#' @export
#'
#' @examples
#' mixed_cor(ggplot2::diamonds)
mixed_cor <- function(data){
  # check if input is a data frame
  if(!is.data.frame(data)){
    stop("The input must be a data frame!")
  }

  # convert to tibble and remove NA rows
  data <- data |> tibble::as_tibble() |> stats::na.omit()

  # check if there are enough rows to compute a correlation
  if(nrow(data) <= 1){
    stop("There are not enough rows without NAs to compute correlation!")
  }

  # initialize correlation matrix
  p <- ncol(data)
  corrmat <- matrix(NA, p, p, dimnames = list(colnames(data), colnames(data)))

  # iterate over all pairs of variables
  for(i in 1:p){
    for(j in 1:p){
      x <- data[,i,drop=TRUE]
      y <- data[,j,drop=TRUE]
      # compute correlation between x and y depending on their types
      corrmat[i,j] <-
        if(is.numeric(x) & is.numeric(y)) {
          stats::cor(x, y, method = "pearson")
        } else if (is.numeric(x) & is.ordered(y)) {
          stats::cor(x, as.numeric(y), method = "spearman")
        } else if(is.numeric(x) & !is.ordered(y)){
          stats::summary.lm(stats::lm(x ~ y))$r.squared
        } else if(is.ordered(x) & is.numeric(y)){
          stats::cor(as.numeric(x), y, method = "spearman")
        } else if(is.ordered(x) & is.ordered(y)){
          stats::cor(as.numeric(x), as.numeric(y), method = "spearman")
        } else if(!is.ordered(x) & is.numeric(y)){
          stats::summary.lm(stats::lm(y ~ x))$r.squared
        } else{
          cramerV(x,y)
        }
    }
  }

  # return
  corrmat
}

# compute Cramer's V (code from https://rdrr.io/cran/rcompanion/src/R/cramerV.r)
cramerV <- function(x, y){
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){y=as.vector(y)}
  N      <- length(x)
  Chi.sq <- suppressWarnings(stats::chisq.test(x, y, correct=FALSE)$statistic)
  Phi    <- Chi.sq / N
  Row    <- length(unique(x))
  C      <- length(unique(y))
  CV     <- sqrt(Phi / min(Row-1, C-1))
  CV |> unname()
}
