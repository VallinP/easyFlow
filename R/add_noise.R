#'#########################################################################################
#' Function to  add noise to raw data.
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @importFrom stats rnorm
#' @import magrittr
#'
#' @param df a dataframe/table of numeric values
#' @param clustering.markers an array of character strings defining the colname to use
#'
#' @return df a dataframe/table of numeric values

add_noise <- function(df, clustering.markers){

  dups <- duplicated(df)

  i <- 1
  if(sum(dups) >= 1){
    message("Some dups events were detected. Adding noise...")
    while (sum(dups) >= 1){

      df2 <- df[dups, clustering.markers]
      df3 <- df[ , clustering.markers]

      rand <- stats::rnorm(n=ncol(df2)*nrow(df2), mean = 0, sd = 1) / 1000000
      rand <- matrix(rand,
                     nrow = nrow(df2),
                     ncol = ncol(df2))

      df2 <- df2 + rand
      df3[dups, ] <- df2

      dups <- sum(duplicated(df3))
      #print(paste0("Try ", i, " dups: ", dups))
      i <- i+1
    }

    colnames(df2)
    df <- cbind(df2, df[, !colnames(df) %in% clustering.markers])
  }

  return(df)
}
