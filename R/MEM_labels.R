#'#########################################################################################
#' Function to calculate MEM (Marker enrichment modeling) score.
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @import magrittr
#' @importFrom MEM MEM
#'
#' @param dir a character string specifying the project directory
#'
#' @export

MEM_labels <- function(dir = NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  {
    message("\n")
    message(" ########## Running Marker Enrichment Modeling ########## ")

    message("Loading data...")
    source(paste0(dir,"/attachments/easyFlow.R"))
    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))

    MEM_outdir <- paste0(dir, "/results/MEM/")
    #current_dir <- getwd()

    suppressWarnings(dir.create(paste0(MEM_outdir,"/output files"), recursive = TRUE))
    #setwd(MEM_outdir)
    #getwd()
  }

  # Subsampling clusters
  {
    # Select event_per_clus IDs in each clusters (1000events)
    message("Extracting reference clusters (deep clusters profile)...")
    ID <- NULL
    event_per_clus <- 1000
    message(paste0("metaclusters : ", length(unique(metacluster)) ))
    message(paste0("max events per metacluster : ", event_per_clus ))
    ID_cluster <- sort(unique(metacluster))

    for(i in ID_cluster){
      if(sum(as.numeric(metacluster) == i) > 0){
        clustering.map.id <- which(as.numeric(metacluster) == i)
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
    data <- data_table[ID,]
    cluster_ids <- metacluster[ID]

    # Check nrow for each ref cluster
    #for(lev in levels(factor(metacluster))){print(paste0("metacluster ", lev, " nrow ", sum(subsample_metacluster==lev)))}

    ###### TO DO ############
    # Validate the subsampling
    #subsample_group <- c(rep(1, nrow(subsample_df)), rep(2, nrow(data_table)))
    #subsample_df <- rbind(subsample_df[, clustering.markers], data_table[, clustering.markers])

    #for(i in ID_cluster){
    #  subsample_df2 <- subsample_df[c(subsample_metacluster == i, metacluster == i), ]
    #  subsample_group2 <- subsample_group[c(subsample_metacluster == i, metacluster == i)]
    #} # End for cluster


  } # End reference clusters


  {
    # Format data
    expr <- data
    data <- cbind(data[, clustering.markers], cluster_ids)
    colnames(data) <- c(clustering.markers, "cluster")


    ## Run MEM. For detailed explanation of the arguments, enter ?MEM in console.
    # To test MEM on example data (normal human PBMC subsets; see ?PBMC), modify the following line of code to read MEM(PBMC,transform=FALSE...)
    # If any of the arguments choose.markers, choose.ref, and/or rename.markers = TRUE you will be prompted to enter the column indices of markers, reference cluster number, or new marker names, respectively
    message("Computing MEM scores...")
    #?MEM
    MEM_values <- MEM::MEM(data,
                           #transform = FALSE,
                           #cofactor = 1,
                           #choose.markers = FALSE,
                           #choose.ref = FALSE,
                           #rename.markers = FALSE,
                           #file.is.clust = FALSE,
                           #add.fileID = FALSE,
                           IQR_thresh = "auto",
                           output_prescaled_MEM = F)
  }

  {
    #MEM values will be written to text file in "output files" folder (created by build.heatmaps(), can be found in your working directory)
    heatmap_data <- (MEM_values[[5]])[[1]]

    MEM_vals_scale <- as.matrix(round(heatmap_data, 0))

    create.labels <- function (MEM_vals_scale, display.thresh, heatmap_data) {
      posRownamesList = list()
      posRownamesMatrix = matrix()
      negRownamesMatrix = matrix()
      for (i in 1:nrow(heatmap_data)) {
        posMarkers = MEM_vals_scale[i, ][MEM_vals_scale[i, ] >=
                                           display.thresh]
        negMarkers = MEM_vals_scale[i, ][MEM_vals_scale[i, ] <
                                           0 & abs(MEM_vals_scale[i, ]) >= display.thresh]
        posMarkersOrdered = posMarkers[order(posMarkers, decreasing = TRUE)]
        negMarkersOrdered = negMarkers[order(abs(negMarkers),
                                             decreasing = TRUE)]
        posVector = rep("+", length(posMarkers))
        negVector = rep("-", length(negMarkers))
        posSignedVec = as.matrix(paste(posVector, abs(round(posMarkersOrdered,
                                                            2)), sep = ""))
        posMarkersRanked = as.matrix(names(posMarkersOrdered))
        posRownames = t(paste(posMarkersRanked, posSignedVec,
                              sep = ""))
        posRownamesMatrix[i] = apply(posRownames, 1, paste,
                                     collapse = "_")
        negSignedVec = as.matrix(paste(negVector, abs(round(negMarkersOrdered,
                                                            2)), sep = ""))
        negMarkersRanked = as.matrix(names(negMarkersOrdered))
        negRownames = t(paste(negMarkersRanked, negSignedVec,
                              sep = ""))
        negRownamesMatrix[i] = apply(negRownames, 1, paste,
                                     collapse = "_")
      }
      new_rownames = paste(row.names(MEM_vals_scale), ":", posRownamesMatrix,
                           negRownamesMatrix, sep = "_")
      return(new_rownames)
    }

    new_rownames <- create.labels(MEM_vals_scale, display.thresh=0,
                                  heatmap_data)

    write.table(new_rownames, paste0(MEM_outdir, "output files/enrichment score-rownames.txt"),
                sep = ",", row.names = FALSE)
  }


  {
    # Save MEM results
    message("Saving MEM results...")
    #setwd(current_dir)
    #getwd()

    save(
      MEM_values,
      MEM_outdir,

      # RData compression parameters
      compress = TRUE, compression_level=1,

      file=paste0(dir, sep="/Rdata/MEM.RData"))
  }

  {
    MEM_metacluster <- c(1:nrow(MEM_vals_scale))
    MEM.df <- cbind(MEM_metacluster, MEM_vals_scale)
    colnames(MEM.df) <- c("Metacluster", colnames(MEM_vals_scale))

    # Sort MEM.df for each markers per their MEM score
    for (i in 3: length(colnames(MEM.df))){
      MEM.df <- MEM.df[order(as.numeric(MEM.df[,i])),]
    } # End for length(clustering.markers)

    # Set an ID to each row to later sort
    i <- data.frame(matrix(c(1:nrow(MEM.df)), ncol=1))
    colnames(i) <- "ID"
    MEM.df <- cbind(i, MEM.df)
    MEM.df

    #MEM.df <- MEM.df[order(as.numeric(MEM.df[,2])),]

    # Extract events per MEM.score as a df
    i <- data.frame(rep(x = NA, times = nrow(MEM.df)))
    MEM.df <- cbind(MEM.df, i)
    colnames(MEM.df) <- c(colnames(MEM.df[,c(1:ncol(MEM.df)-1)]), "MEM")
    for (i in 1:nrow(MEM.df)){
      for (j in 3:(ncol(MEM.df)-1)){
        if(j==3){
          if(MEM.df[i,j] < 0){
            MEM.df[i,"MEM"] <- paste0(colnames(MEM.df)[j], MEM.df[i,j])
          } else {
            MEM.df[i,"MEM"] <- paste0(colnames(MEM.df)[j], "+", MEM.df[i,j])
          }

        } else {
          if(MEM.df[i,j] < 0){
            MEM.df[i,"MEM"] <- paste0(MEM.df[i,"MEM"], "_", colnames(MEM.df)[j], MEM.df[i,j])
          } else {
            MEM.df[i,"MEM"] <- paste0(MEM.df[i,"MEM"], "_", colnames(MEM.df)[j], "+", MEM.df[i,j])
          }

        }
      }
    }

    # Write a reorganised table for MEM scores
    write.table(MEM.df, file=paste0(dir, "/results/MEM/output files/enrichment score_reorganised.csv"), row.names = F, sep=",")
  }

}
