#'#########################################################################################
#' Save results functions
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @importFrom flowCore write.FCS flowFrame
#'
#' @param dir a character string specifying the project directory
#'
#' @export

save_results <- function(dir = NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  # Unsampling results
  #tsne.coord <- tsne.coord[order(sampling.id),]
  #clus[[1]] <- clus[[1]][order(sampling.id)]
  #metacluster <- metacl[, n.metaclusters][som.cluster]#[order(sampling.id)]

  # Load data
  {
    message("\n")
    message(" ########## Saving results ########## ")
    message("Saving parameters... ")
    suppressWarnings(dir.create(paste0(dir,"/Rdata/")))

    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/som_map.RData"))
    load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))
  } # End load data


  if(do.dimred == TRUE){
    load(paste0(dir,"/Rdata/tsne.RData"))

    # Add new pars to data_table
    new.parameters <- matrix(cbind(as.numeric(tsne.coord[,1]),
                                   as.numeric(tsne.coord[,2]),
                                   as.numeric(som.cluster),
                                   as.numeric(metacluster)),
                             ncol = 4)
    colnames(new.parameters) <- c("tsne_1", "tsne_2", "som_cluster", "cluster")

    exclude.markers <- c(exclude.markers, "tsne_1", "tsne_2", "som_cluster", "cluster")
    exclude.markers.names <- c(exclude.markers.names, "tsne_1", "tsne_2", "som_cluster", "cluster")

  } else {
    # Add new pars to data_table
    new.parameters <- matrix(cbind(as.numeric(som.cluster),
                                   as.numeric(metacluster)),
                             ncol = 2)
    colnames(new.parameters) <- c("som_cluster", "cluster")

    exclude.markers <- c(exclude.markers, "som_cluster", "cluster")
    exclude.markers.names <- c(exclude.markers.names, "som_cluster", "cluster")
  } # End if do_dimred
  #head(new.parameters)


  # Save all objects in a Rdata file
  save(list = c("dir", "FItSNE.dir", "files", "parameters", "metadata",
                "data_type","cofactor",
                "compute.quantiles","quantile_thresh", "compute.scaling",
                "data_table", "new.parameters",
                "transform.markers.names",
                "clustering.markers", "clustering.markers.names",
                "activation.markers", "activation.markers.names",
                "exclude.markers","exclude.markers.names",
                "markers", "markers.names",
                "assignments", "condition1.files", "condition2.files",
                "variable", "status", "lineage_markers", "functional_markers",
                "files_extract", "samples", "ix.samples", "assignments_extract",
                "color_conditions", "ix.bc", "sample_ids", "sample_test", "ID_event"
  ),
  file = paste0(dir,"/Rdata/easyFlow.RData"),
  compress = "gzip"
  )


  # Save fcs files
  {
    message("Saving fcs files with new parameters...")
    data_table <- cbind(data_table, new.parameters)
    output.dir <- paste0(dir,"/results/FCS/")
    suppressWarnings(dir.create(output.dir))

    for (i in 1:length(files)){
      message(paste0("Processing ",files[i]))

      # Read data
      data_table.extract <- data_table[data_table[,"ID_file"]==i,]

      # Write fcs file
      tube.name <- gsub(".fcs","_easyFlow.fcs",files[i])
      flowCore::write.FCS(flowCore::flowFrame(data_table.extract), filename = paste0(output.dir , tube.name))

    } # End For Files

  } # End save fcs files

}
