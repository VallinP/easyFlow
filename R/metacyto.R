#'#########################################################################################
#' Metacyto
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'
#' Adapted from
#' Hu et al., Cytometry Part A, 2018
#' https://linkinghub.elsevier.com/retrieve/pii/S2211-1247(18)31080-5
#' https://bioconductor.org/packages/release/bioc/html/MetaCyto.html
#'########################################################################################
#'
#' @source https://linkinghub.elsevier.com/retrieve/pii/S2211-1247(18)31080-5
#' @source https://bioconductor.org/packages/release/bioc/html/MetaCyto.html
#'
#' @import magrittr
#' @import MetaCyto
#'
#' @param dir a character string specifying the project directory
#'
#' @export


metacyto <- function(dir = NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  message("\n")
  message(" ########## Running MetaCyto ########## ")
  message("Loading data...")
  load(paste0(dir,"/Rdata/easyFlow.RData"))
  load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))

  # Subsampling clusters
  {

    # Select event_per_clus IDs in each clusters (1000events)
    message("Extracting reference clusters (deep clusters profile)...")
    ID <- NULL
    event_per_clus <- 100
    message(paste0("clusters : ", length(unique(som.cluster)) ))
    message(paste0("max events per metacluster : ", event_per_clus ))
    ID_cluster <- sort(unique(som.cluster))

    for(i in ID_cluster){
      if(sum(as.numeric(som.cluster) == i) > 0){
        clustering.map.id <- which(as.numeric(som.cluster) == i)
      } else {
        clustering.map.id <- NULL
      }

      if(is.null(ID)){

        if(length(clustering.map.id) > event_per_clus){
          ID <- sample(clustering.map.id, event_per_clus, replace = FALSE)
        }
        if(length(clustering.map.id) <= event_per_clus && length(clustering.map.id) > 0){
          ID <- clustering.map.id
        }

      } else {

        if(length(clustering.map.id) > event_per_clus){
          ID <- c(ID, sample(clustering.map.id, event_per_clus, replace = TRUE))
        }
        if(length(clustering.map.id) <= event_per_clus && length(clustering.map.id) > 0){
          ID <- c(ID, clustering.map.id)
        }

      }
    } # End for cluster

    # Cluster Map
    subsample_df <- data_table[ID,]
    subsample_som.cluster <- som.cluster[ID]

    # Check nrow for each ref cluster
    #for(lev in levels(factor(som.cluster))){print(paste0("som.cluster ", lev, " nrow ", sum(subsample_som.cluster==lev)))}

    ###### TO DO ############
    # Validate the subsampling
    #subsample_group <- c(rep(1, nrow(subsample_df)), rep(2, nrow(data_table)))
    #subsample_df <- rbind(subsample_df[, clustering.markers], data_table[, clustering.markers])

    #for(i in ID_cluster){
    #  subsample_df2 <- subsample_df[c(subsample_som.cluster == i, som.cluster == i), ]
    #  subsample_group2 <- subsample_group[c(subsample_som.cluster == i, som.cluster == i)]
    #} # End for cluster

  } # End reference clusters

  # MetaCyto Labels
  message("Resolving metaclusters with MetaCyto...")

  # label each clusters
  #metacyto.res <- MetaCyto::labelCluster(fcsFrame = flowFrame(data_table[, clustering.markers]), clusterList = som.cluster)
  metacyto.res <- MetaCyto::labelCluster(fcsFrame = flowFrame(subsample_df[, clustering.markers]), clusterList = subsample_som.cluster)

  levels(factor(metacyto.res$clusterLabel))
  metacyto.res$cutoff

  # Add a user-defined label (CCR7+ CD8+ T cells) to the cluster_label. This step allows users
  # to include their favorate cell types into the analysis.
  #cluster_label$clusterLabel=c(cluster_label$clusterLabel,"CD3+|CD4-|CD8B+|CCR7+")
  #MetaCyto::labelSummary(metacyto.res, minStudy = 2)

  # Find the hyper-ractangle clusters corresponding to the labels
  cluster_list <- MetaCyto::searchCluster(fcsFrame=flowFrame(data_table[, clustering.markers]),
                                          clusterLabel=metacyto.res$clusterLabel,
                                          cutoff=metacyto.res$cutoff)
  #cluster_list$clusterList

  # Find ungated events
  ID <- NULL
  for(i in 1:length(metacyto.res$clusterLabel)){
    ID <- c(ID, cluster_list$clusterList[[i]])
  }
  ID <- sort(unique(ID))
  length(ID)

  ID_ungated <- c(1:nrow(data_table))
  '%!in%' <- function(x,y)!('%in%'(x,y))
  ID_ungated <- ID_ungated[ID_ungated %!in% ID]

  if(sum(ID_ungated) != 0){
    cluster_list$clusterList[["ungated"]] <- ID_ungated
  }

  # Format metacyto metaclusters
  metacyto_names <- rev(sort(levels(factor(names(cluster_list$clusterList)))))
  names(metacyto_names) <- c(1: length(metacyto_names))
  metacyto_label <- rep(NA, nrow(data_table))
  metacyto_mc <- rep(NA, nrow(data_table))
  ID <- 0

  i <- names(metacyto_names)[2]
  for(i in names(metacyto_names)){
    ix <- cluster_list$clusterList[[metacyto_names[i]]]
    ix <- ix[ix %!in% ID]
    metacyto_label[ix] <- as.character(metacyto_names[i])
    metacyto_mc[ix] <- as.numeric(i)
    ID <- unique(sort(c(ID, ix)))
  }

  #summary(metacyto_mc[metacyto_mc != 1])
  #hist(metacyto_mc[metacyto_mc != 1])

  # Save Hierarchical Consensus Clustering
  message("Saving metacyto results...")
  suppressWarnings(dir.create(paste0(dir,"/Rdata/")))
  save(list = c("metacyto.res",
                "cluster_list",
                "metacyto_label",
                "metacyto_mc",
                "metacyto_names"
  ),
  file = paste0(dir,"/Rdata/metacyto.RData"),
  compress = "gzip"
  )

  # Format metacyto results
  metacluster <- metacyto_mc
  n.metaclusters <- length(unique(metacyto_label))


}
