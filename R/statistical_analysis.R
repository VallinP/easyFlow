#'#########################################################################################
#' Statistical tests functions and plots
#'
#' Adapted from :
#' Guillaume Gautreau et al., Bioinformatics, March 2017
#' https://academic.oup.com/bioinformatics/article/33/5/779/2662307
#' https://github.com/tchitchek-lab/SPADEVizR
#'
#' Adapted from :
#' Nowicka M, Crowell H, Robinson M (2019).
#' cytofWorkflow: CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets.
#' https://f1000research.com/articles/6-748
#' https://github.com/markrobinsonuzh/cytofWorkflow.
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @source https://academic.oup.com/bioinformatics/article/33/5/779/2662307
#' @source https://github.com/tchitchek-lab/SPADEVizR
#' @source https://f1000research.com/articles/6-748
#' @source https://github.com/markrobinsonuzh/cytofWorkflow
#'
#' @param dir a character string specifying the project directory
#'
#' @export

#SPADEVizR importResultsFromFCS assignContext qcSmallClusters qcUniformClusters identifyAC identifyDAC identifyCC classifyAbundanceProfiles countViewer heatmapViewer phenoViewer boxplotViewer kineticsViewer streamgraphViewer MDSViewer distogramViewer

statistical_analysis <- function(dir = NULL){

  {
    if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

    message("\n")
    message(" ########## Running statistical analysis ########## ")

  }

  #############
  # SPADEVizR #
  #############
  easyFlow:::SPADEvizR(dir)

  #################
  # CytofWorkFlow #
  #################
  easyFlow:::CytofWorkFlow(dir)

}
