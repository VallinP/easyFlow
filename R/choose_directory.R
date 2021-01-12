#'#########################################################################################
#' Function to choose a project directory.
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @import magrittr
#' @import utils
#' @importFrom tcltk tk_choose.dir
#'
#' @param caption a character string
#'
#' @return dir

choose_directory <- function(caption = "Select a directory") {
  setwd("~/")

  if (exists('utils::choose.dir')) {
    dir <- choose.dir(caption)
  } else {
    dir <- tcltk::tk_choose.dir(caption)
  }
  return(dir)
}
