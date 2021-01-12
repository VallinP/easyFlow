##########################################################################################
#' Integrated and easy to run workFlow
#' for unsupervised data analysis.
#' raw data > preprocessing > clustering > dimensionality reduction > statistical analysis
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#########################################################################################
#'
#' @export

easyFlow.analysis <- function(){

  ############# Preprocessing data #############
  {
    dir <- easyFlow::preprocessing()
  } # End Preprocessing data


  ############# SOM #############
  {
    easyFlow::som_clustering(dir = dir)
  } # End SOM


  ############# Hierarchical clustering #############
  {
    easyFlow::metaClustering_HC(dir = dir)
  } # End hierarchical clustering


  ############# Dimensionality reduction #############
  if(do.dimred == TRUE){
    easyFlow::dimensional_reduction(dir = dir)
  } # end FItSNE


  ############# Set metaclusters #############
  {
    easyFlow::set_metaclusters(dir = dir)
  } # end save results


  ############# Marker Enrichment Modeling #############
  if(do_mem == TRUE){
    easyFlow::MEM_labels(dir = dir)
  } # End MEM


  ############# MetaCyto #############
  if(do_metacyto == TRUE){
    metacyto(dir = dir)
  } # End MetaCyto


  ############# Save results #############
  {
    easyFlow::save_results(dir = dir)
  } # end save results


  ############# Reports & statistics #############
  {
    easyFlow::report(dir = dir)
    if(do_mem == TRUE){easyFlow::MEM_report(dir = dir)}
    if(compute_stat == TRUE){easyFlow:::statistical_analysis(dir = dir)}
  }

} # end easyFlow analysis
