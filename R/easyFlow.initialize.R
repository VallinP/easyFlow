#'#########################################################################################
#' Initialize easyFlow
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @import magrittr
#' @importFrom flowCore read.flowSet read.FCS
#'
#' @export

easyFlow.initialize <- function(){

  caption.value <- "Select your project directory"
  dir <- easyFlow:::choose_directory(caption.value)
  fcs.file <- list.files(dir, pattern = "*.fcs")[1]

  if(is.na(fcs.file)){stop("No fcs files found.")}

  # Files to use
  files <- list.files(dir, pattern = "*.fcs$", recursive = F)
  files

  # Attachment directory
  suppressWarnings(dir.create(paste0(dir, "/attachments/")))

  # FItSNE directory
  FItSNE.dir <- easyFlow:::choose_FItSNE.dir(dir)

  # Panel file
  message("Reading fcs files...")
  fs <- NULL
  if(!suppressMessages(class(try(flowCore::read.flowSet(paste0(dir,"/",files)), silent=TRUE)) == "try-error")){
    fs <- suppressMessages(flowCore::read.flowSet(files = paste0(dir,"/",files),
                                            transformation=FALSE, truncate_max_range = FALSE))
    ff <- flowCore::read.FCS(paste0(dir,"/",files)[1])
    ff
    df <- data.frame(ff@parameters@data)
    df[,c("range", "minRange", "maxRange")] <- NULL
    rownames(df) <- NULL

    df <- cbind(c(1:nrow(df)),df)

    df[as.logical(is.na(df[,"desc"])),"desc"] <- df[as.logical(is.na(df[,"desc"])),"name"]
    colnames(df) <- c("Name", "Parameter", "Antigen")

    tmp <- data.frame(matrix(data = NA, nrow = nrow(df), ncol = 4))
    colnames(tmp) <- c("colsToTransform", "clustering.cols", "activation.cols", "exclude.cols")
    df <- cbind(df, tmp)

  } else {
    stop("The fcs files provided aren't from the same panel")
  }

  if(!file.exists(paste0(dir, "/attachments/parameters.csv"))){
    message("No previous parameters.csv file detected. Writing a new one...")
    write.table(df, file = paste0(dir,"/attachments/parameters.csv"), sep = ",", row.names = F, col.names = T)
    message("Please, complete parameters.csv file before next step")
  } else {
    if(all(as.character(df[,"Parameter"]) %in% read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"Parameter"]) &&
       all(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"Parameter"] %in% df[,"Parameter"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"Antigen"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"colsToTransform"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"clustering.cols"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"activation.cols"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"exclude.cols"]) ){
      message("Existing parameter.csv file detected.")
    } else {
      message("Previous parameters.csv file detected, but invalid. Writing a new one...")
      write.table(df, file = paste0(dir,"/attachments/parameters.csv"), sep = ",", row.names = F, col.names = T)
      message("Please, complete parameters.csv file before next step")
    }
  }

  # Metadata file
  df2 <- matrix(data = NA, ncol = 8, nrow = length(files))
  colnames(df2) <- c("Name", "bc", "ind", "tp", "variable", "status", "FlowJo_workspace",	"gate")
  df2[, "Name"] <- files

  if(!file.exists(paste0(dir, "/attachments/metadata.csv"))){
    message("No previous metadata.csv file detected. Writing a new one...")
    write.table(df2, file = paste0(dir,"/attachments/metadata.csv"), sep = ",", row.names = F, col.names = T)
    message("Please, complete metadata.csv file before next step")
  } else {
    if(all(files %in% read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"Name"]) &&
       all(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"Name"] %in% files) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"bc"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"ind"]) ){
      message("Existing metadata.csv file detected.")
    } else {
      message("Previous metadata.csv file detected, but invalid. Writing a new one...")
      write.table(df2, file = paste0(dir,"/attachments/metadata.csv"), sep = ",", row.names = F, col.names = T)
      message("Please, complete metadata.csv file before next step")
    }
  }


  # easyFlow.parameters
  if(!file.exists(paste0(dir, "/attachments/easyFlow.R"))){
    message("easyFlow.R initialized with default parameters")
    sink(file=paste0(dir, "/attachments/easyFlow.R"))


    cat("# Preprocessing")
    cat("\n")
    cat("\n")
    cat("## data type")
    cat("\n")
    cat("data_type <- \"FCM\" ")
    cat("\n")
    cat("seed <- 123456")
    cat("\n")
    cat("\n")

    cat("## FlowAI")
    cat("\n")
    cat("do_flowai <- TRUE")
    cat("\n")
    cat("do_pregating <- TRUE")
    cat("\n")
    cat("\n")

    cat("## data transformation # cofactor")
    cat("\n")
    cat("cofactor <- 5")
    cat("\n")
    cat("\n")

    cat("compute.quantiles <- F")
    cat("\n")
    cat("quantile_thresh <- 1/100")
    cat("\n")
    cat("compute.scaling <- F")
    cat("\n")
    cat("\n")

    cat("# FlowSOM")
    cat("\n")
    cat("flowSOM_dim <- 20")
    cat("\n")
    cat("grid.type <- \"rectangular\" ")
    cat("\n")
    cat("use_map <- \"train\" ")
    cat("\n")
    cat("n.metaclusters <- 30")
    cat("\n")
    cat("HCC_method <- \"average\" #ward.D, simple, average, complete")
    cat("\n")
    cat("use_elbow_point <- FALSE")
    cat("\n")
    cat("HC_unimodal <- FALSE")
    cat("\n")
    cat("do_mem <- FALSE")
    cat("\n")
    cat("do_metacyto <- FALSE")
    cat("\n")

    cat("\n")
    cat("# FItSNE")
    cat("\n")
    cat("do.dimred <- TRUE")
    cat("\n")
    cat("perplexity.value <- 50")
    cat("\n")
    cat("theta.value <- 0.5")
    cat("\n")
    cat("freedom.value <- 1")
    cat("\n")
    cat("n.iter <- 3000")
    cat("\n")
    cat("init.matrix <- FALSE")
    cat("\n")
    cat("tsne_prediction <- \"none\" # none, complete, k1init, k1final")
    cat("\n")
    cat("tsne_pred_sampling <- \"CS\" # cluster_stratified or \"DCP\" deep_cluster_profile ")
    cat("\n")
    cat("tsne_pred_n_events <- 100000")
    cat("\n")
    cat("tsne_algo <- TRUE")  # FItSNE = TRUE, BH-tSNE = FALSE
    cat("\n")

    cat("\n")
    cat("# Statistics")
    cat("\n")
    cat("compute_stat <- FALSE")
    cat("\n")
    cat("\n")
    cat("condition1 <- \"NA\"")
    cat("\n")
    cat("condition2 <- \"NA\"")
    cat("\n")
    cat("\n")
    sink()



  } else {
    message("easyFlow.R already exists")
  }

  message("easyFlow was successfully initialized.")

} # End easyFlow.initialize()
