#'#########################################################################################
#' Efficient scale function.
#'
#' https://www.r-bloggers.com/a-faster-scale-function/
#'########################################################################################
#'
#' @source https://www.r-bloggers.com/a-faster-scale-function/
#'
#' @importFrom Matrix colMeans
#' @importFrom matrixStats colSds
#'
#' @param x a table/dataframe
#' @param center a logical specifying if the data should be centered
#' @param scale a logical specifying if the data should be scaled
#' @param add_attr  a logical
#' @param rows an array of numeric values
#' @param cols an array of numeric values
#'
#' @return x a table/dataframe

colScale <- function(x,
                     center = TRUE,
                     scale = TRUE,
                     add_attr = TRUE,
                     rows = NULL,
                     cols = NULL) {

  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  ################
  # Get the column means
  ################
  cm = Matrix::colMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = matrixStats::colSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t( (t(x) - cm) / csd )
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
} # end scale
