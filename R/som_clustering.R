##########################################################################################
#' Clustering using self organizing map
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#########################################################################################
#'
#' @importFrom stats prcomp
#' @importFrom FlowSOM SOM ReadInput BuildSOM BuildMST
#' @importFrom kohonen som supersom
#' @importFrom Rsomoclu Rsomoclu.train Rsomoclu.kohonen
#'
#' @param dir a character string specifing the project directory
#'
#' @export

som_clustering <- function(dir = NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  message("\n")
  message(" ########## Running self organizing map ########## ")

  source(paste0(dir,"/attachments/easyFlow.R"))
  load(paste0(dir,"/Rdata/easyFlow.RData"))
  clustering_method <- "FlowSOM"

  df <- data_table[, clustering.markers]

  # set seed
  set.seed(seed)

  # Compute PCA
  PCA.res <- stats::prcomp(df)


  # Training
  if(use_map == "train"){

    # Build the self-organizing map

    if(clustering_method == "FlowSOM"){
      message("SOM method : FlowSOM")
      flowsom.res <- FlowSOM::SOM(data = as.matrix(df),
                                  xdim = flowSOM_dim, ydim = flowSOM_dim,
                                  rlen = 10, mst = 1, alpha = c(0.05, 0.01),
                                  distf = 2, # Distance function (1=manhattan, 2=euclidean, 3=chebyshev, 4=cosine)
                                  #radius = stats::quantile(nhbrdist, 0.67) * c(1, 0),
                                  init = FALSE, silent = FALSE, codes = NULL, importance = NULL
      )

      if(F){
        # df <- data_Rtsne[, clustering.markers]
        flowsom.res <- FlowSOM::ReadInput(flowFrame(df), compensate = FALSE, spillover = NULL,
                                          transform = FALSE, toTransform = NULL,
                                          scale = FALSE, scaled.center = FALSE, scaled.scale = FALSE, silent = FALSE)
        flowsom.res <- FlowSOM::BuildSOM(fsom = flowsom.res, colsToUse=clustering.markers,
                                         xdim=flowSOM_dim, ydim=flowSOM_dim, rlen=10)

        flowsom.res <- FlowSOM::BuildMST(flowsom.res)

      } #end false

      # Save SOM result
      message("Saving SOM results...")
      suppressWarnings(dir.create(paste0(dir,"/Rdata/")))
      save("flowsom.res",
           file = paste0(dir,"/Rdata/flowsom_map.RData"),
           compress = "gzip"
      )

      rownames(flowsom.res$grid) <- NULL
      rownames(flowsom.res$codes) <- paste0("V", c(1:nrow(flowsom.res$codes)))
      names(flowsom.res$radius[[1]]) <- c("66.66667%", "")
      colnames(flowsom.res$grid) <- c("x", "y")

      flowsom.grid <- list(pts = as.matrix(flowsom.res$grid),
                           xdim = flowsom.res$xdim,
                           ydim = flowsom.res$ydim,
                           topo = "rectangular",
                           neighbourhood.fct = "gaussian",
                           toroidal = FALSE)
      class(flowsom.grid) <- "somgrid"

      som.res <- list(data = list(df),
                      unit.classif = flowsom.res$mapping[,1],
                      distances = flowsom.res$mapping[,2],
                      grid = flowsom.grid,
                      codes = list(flowsom.res$codes),
                      changes = rep(1, flowsom.res$rlen),
                      alpha = flowsom.res$alpha[[1]],
                      radius = flowsom.res$radius[[1]],
                      user.weights = 1,
                      distance.weights = 1,
                      whatmap = 1,
                      maxNA.fraction = 0,
                      dist.fcts = "sumofsquares"
      )
      class(som.res) <- "kohonen"

    }


    if(clustering_method == "som"){
      message("SOM method : Kohonen som")
      som.res <- kohonen::som(
        as.matrix(df),
        rlen = 10, cores = -1,
        alpha = c(0.05, 0.01),
        #radius = quantile(nhbrdist, 2/3),
        grid=somgrid(flowSOM_dim,flowSOM_dim,grid.type)
      )
    }

    if(clustering_method == "supersom"){
      message("SOM method : supersom")
      som.res <- kohonen::supersom(
        as.matrix(df),
        rlen = 10, cores = -1,
        alpha = c(0.05, 0.01),
        #radius = quantile(nhbrdist, 2/3),
        grid=somgrid(flowSOM_dim,flowSOM_dim,grid.type)
      )
    }

    if(clustering_method == "Rsomoclu"){
      message("SOM method : Rsomoclu")
      som.res <- Rsomoclu::Rsomoclu.train(
        as.matrix(df),
        nSomX	= flowSOM_dim,
        nSomY	= flowSOM_dim,
        nEpoch = 10,
        radius0 = 10,
        radiusN = 1,
        radiusCooling = "linear",
        scale0 = 0.05,
        scaleN = 0.01,
        scaleCooling = "linear",
        kernelType = 0, #Kernel type 0: Dense CPU 1: Dense GPU 2: Sparse CPU (default: 0)
        mapType = "planar",
        gridType = grid.type,
        compactSupport = FALSE,
        codebook = NULL,
        neighborhood = "gaussian",
        stdCoeff = 0.5
      )
      som.res <- Rsomoclu::Rsomoclu.kohonen(as.matrix(df),som.res)
    }

    if(clustering_method == "somoclu"){
      message("SOM method : somoclu")
      system(
        "mpirun -np 4 somoclu -k 0 -e 10 -l 0.05 -L 0.01 -x 20 -y 20 /home/vallin/Documents/Mosmann_rare/results/table_data.csv /home/vallin/Documents/Mosmann_rare/results/Mosmann"
      )
    }

    # save som map
    som.cluster <- factor(x=som.res$unit.classif, levels=as.numeric(levels(factor(som.res$unit.classif))))
    clustering.markers.som <- clustering.markers
    clustering.markers.names.som <- clustering.markers.names

    suppressWarnings(dir.create(paste0(dir,"/Rdata/")))
    save(list = c("clustering.markers.som",
                  "clustering.markers.names.som",
                  "som.res", "som.cluster",
                  "flowSOM_dim", "grid.type", "seed",
                  "clustering_method", "use_map"
    ),
    file = paste0(dir,"/Rdata/som_map.RData"),
    compress = "gzip"
    )

  } # end if use_map = train


  if(use_map == "recall"){

    caption.value <- "Select a compatible som map"
    #som.file <- paste0(choose_directory(caption.value),"/som_map.RData")
    message("Loading som map...")
    load(som.file)
    message("Checking compatibility...")
    if(!all(clustering.markers %in% clustering.markers.som ||
            clustering.markers.som %in% clustering.markers ||
            clustering.markers.names %in% clustering.markers.names.som ||
            clustering.markers.names.som %in% clustering.markers.names)){
      stop("This som is not compatible with the current panel")
    } else {
      stop("To Do")
    }

  } # end if use_map = recall

}
