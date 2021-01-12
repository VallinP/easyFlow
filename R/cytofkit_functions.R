#'#########################################################################################
#' Cytofkit's functions
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'
#' Adapted from :
#' Chen et al., 2016, PLOS computational biology
#' https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005112
#'########################################################################################
#'
#' @source https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005112
#'
#' @importFrom methods is
#' @importFrom stats quantile
#' @importFrom flowCore read.flowSet flowSet read.FCS fsApply exprs flowFrame transform
#' @importClassesFrom flowCore flowFrame
#'
#' @param x a dataframe/table of numeric value
#' @param channels an array of character strings specifying the channels to use
#' @param m a numeric value
#' @param q a numeric value
#'
#' @param value a dataframe/table of numeric value
#' @param cofactor a numeric value
#'
#' @return value a dataframe/table of numeric value
# flowCore logicleTransform transformList

autoLgcl <- function(x, channels, m = 4.5, q = 0.05) {
  if (!is(x, "flowFrame"))
    stop("x has to be an object of class \"flowFrame\"")
  if (missing(channels))
    stop("Please specify the channels to be logicle transformed")
  indx <- channels %in% colnames(x@exprs)
  if (!all(indx))
    stop(paste("Channels", channels[!indx], "were not found in the FCS file.\n ",
               sep = " "))

  trans <- lapply(channels, function(p) {
    data <- x@exprs[, p]
    w <- 0
    t <- max(data)
    ndata <- data[data < 0]
    ## use 1.5 * IQR to filter outliers in negative values
    nThres <- quantile(ndata, 0.25) - 1.5 * IQR(ndata)
    ndata <- ndata[ndata >= nThres]
    transId <- paste(p, "autolgclTransform", sep = "_")
    if (length(ndata)) {
      r <- .Machine$double.eps + quantile(ndata, q)
      ## Check to avoid failure of negative w
      if (10^m * abs(r) <= t) {
        w <- 0
      } else {
        w <- (m - log10(t/abs(r)))/2
        if(is.nan(w) || w>2) {
          warning(paste0("autoLgcl failed for channel: ", p, "; using default logicle transformation!"))
          w <- 0.1
          t <- 4000
          m <- 4.5
        }
      }
    }

    flowCore::logicleTransform(transformationId = transId,
                     w = w, t = t, m = m, a = 0)
  })

  flowCore::transformList(channels, trans)

} # end autological

cytofAsinh <- function(value, cofactor = cofactor) {
  value <- value-1
  loID <- which(value < 0)
  if(length(loID) > 0)
    value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  value <- value / cofactor
  value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
  return(value)
} # end cytofAsinh


