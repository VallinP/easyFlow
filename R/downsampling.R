DCP <- function(dir =NULL, downsample_events = NULL){

  # Load data
  {
    if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

    message("\n")
    message(" ########## Downsampling data ########## ")

    message("Loading data...")
    if(is.null(downsample_events)){source(paste0(dir,"/attachments/easyFlow.R"))}
    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/som_map.RData"))
    source(FItSNE.dir, chdir=T)

    subsample_data_table <- NULL
    set.seed(seed)

  } # end load data


  # Downsampling data
  {
    # Select event_per_clus IDs in each clusters (100 clusters x 500 events = 50.000events)
    message("Extracting reference clusters (deep clusters profile)...")
    ID <- NULL
    i <- sort(unique(som.cluster))[1]
    event_per_clus <- floor(downsample_events / length(as.numeric(levels(som.cluster))))
    message(paste0("clusters : ", length(as.numeric(levels(som.cluster))) ))
    message(paste0("max events per cluster : ", event_per_clus ))
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
    subsample_data_table <- data_table[ID,]
    subsample_som.cluster <- som.cluster[ID]

  } # End downsampling data

} # end DCP function
