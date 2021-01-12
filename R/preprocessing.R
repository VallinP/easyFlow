#'#########################################################################################
#' Preprocessing data
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @import magrittr
#' @importFrom flowCore read.flowSet flowSet read.FCS fsApply exprs flowFrame transform
#' @importFrom stringr str_replace_all
#' @importFrom matrixStats colQuantiles
#'
#' @return  dir a character string specifying the project directory
#'
#' @export

preprocessing <- function(){

  # Check input
  {
    message("\n")
    caption.value <- "Select your project directory"
    message(caption.value)
    dir <- easyFlow:::choose_directory(caption.value)

    # Check fcs files
    fcs.file <- list.files(dir, pattern = "*.fcs")[1]
    if(is.na(fcs.file)){stop("No fcs files found. NB : Only minus .fcs pattern is accepted")}

    # Files to use
    files <- list.files(dir,
                        pattern = "*.fcs$",
                        recursive = F,
                        full.names = FALSE)
    files

    # Check file names
    invalid_char <- FALSE
    # castNum function from https://stat.ethz.ch/pipermail/r-help/2008-September/172617.html
    castNum <- function(n) {
      suppressWarnings(x<-as.numeric(as.character(n)))
      if (is.na(x)){
        return(n)
      }else{
        return(x)
      }
    }
    for(i in 1:length(files)){
      if(is.numeric(castNum(substr(files[i], 1, 1)))){invalid_char <- TRUE}
    }
    if(invalid_char){stop("fcs filenames cannot begin by a number")}

    # Punctuation characters "[[:punct:]]"; !"#$%&’()*+,-./:;<=>?@[]^_`{|}~
    #  Patrice Vallin List (special chars & metachars)
    pattern <- "²|&|é|-|è|ç|à|=|~|#|@|¤¨£|ù|%|µ|<|>|,|;|:|/|§|!|\\.|\\$|\\*|\\+|\\?|\\^|\\||\\(|\\)|\\{|\\}|\\[|\\]"
    if(TRUE %in% grepl(pattern, gsub(".fcs$", "", files))){stop("fcs filenames cannot contain special charaters and/or metacharacters. Please rename them and re-initialize before to continue.")}

    # Checking FItSNE dir
    message("\n")
    FItSNE.dir <- easyFlow:::choose_FItSNE.dir(dir)
    source(FItSNE.dir, chdir=T)

    # Check parameters.csv
    if(!file.exists(paste0(dir, "/attachments/parameters.csv"))){
      stop("parameters.csv file not found!")
    } else {
      parameters <- read.csv(file=paste0(dir,"/attachments/parameters.csv"))
    }

    # Check metadata.csv
    if(!file.exists(paste0(dir, "/attachments/metadata.csv"))){
      stop("metadata.csv file not found!")
    } else {
      metadata <- read.csv(file=paste0(dir,"/attachments/metadata.csv"))
    }

    # Load easyFlow.R
    if(!file.exists(paste0(dir, "/attachments/easyFlow.R"))){
      stop("easyFlow.R file not found!")
    } else {
      source(paste0(dir,"/attachments/easyFlow.R"))
    }

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

    # Check parameter.csv file
    if(file.exists(paste0(dir, "/attachments/parameters.csv")) &&
       all(as.character(df[,"Parameter"]) %in% read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"Parameter"]) &&
       all(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"Parameter"] %in% df[,"Parameter"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"Antigen"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"colsToTransform"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"clustering.cols"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"activation.cols"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/parameters.csv"))[,"exclude.cols"]) ){
      message("Valid parameter.csv file detected.")
    } else {
      stop("Invalid parameter.csv file.")
    }

    # Check Metadata file
    df2 <- matrix(data = NA, ncol = 6, nrow = length(files))
    colnames(df2) <- c("Name", "bc", "ind", "tp", "variable", "status")
    df2[, "Name"] <- files

    if(file.exists(paste0(dir, "/attachments/metadata.csv")) &&
       all(tolower(files) %in% tolower(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"Name"])) &&
       all(tolower(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"Name"]) %in% tolower(files)) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"bc"]) &&
       !anyNA(read.csv(file=paste0(dir,"/attachments/metadata.csv"))[,"ind"]) ){
      message("Valid metadata.csv file detected.")
    } else {
      stop("Invalid metadata.csv file.")
    }

    # Check easyFlow.R file
    # easyFlow.parameters
    if(!file.exists(paste0(dir, "/attachments/easyFlow.R")) ||
       # FlowAI & pregating
       !class(do_flowai) == "logical" || !class(do_pregating) == "logical" ||
       # Preprocessing
       !(data_type == "FCM" || data_type == "CYTOF") || !class(cofactor) == "numeric" ||
       !class(seed) == "numeric" || !class(compute.quantiles) == "logical" ||
       !class(quantile_thresh) == "numeric" || !class(compute.scaling) == "logical" ||
       # FlowSOM
       !class(flowSOM_dim) == "numeric" || !grid.type == "rectangular" ||
       !use_map == "train" || !class(n.metaclusters) == "numeric" ||
       # FItSNE
       !class(do.dimred) == "logical" ||!class(perplexity.value) == "numeric" || !class(theta.value) == "numeric" ||
       !class(freedom.value) == "numeric" || !class(n.iter) == "numeric" ||
       !class(init.matrix) == "logical" ||
       !(tsne_prediction == "none" || tsne_prediction == "k1init" ||
         tsne_prediction == "k1final" ||tsne_prediction == "complete") ||
       !(tsne_pred_sampling == "CS" || tsne_pred_sampling == "DCP") ||
       !(class(tsne_pred_n_events)  == "numeric" && tsne_pred_n_events >= 100) ||
       !class(tsne_algo) == "logical" ||
       # Statistics
       !class(compute_stat) == "logical"){
      stop("Invalid easyFlow.R file")
    } else {
      message("Valid easyFlow.R file detected.")
    }

    # Retrieve parameters
    transform.markers.names <- as.character(parameters[which(parameters[,"colsToTransform"]==1),"Parameter"])
    clustering.markers.names <- as.character(parameters[which(parameters[,"clustering.cols"]==1),"Parameter"])
    activation.markers.names <- as.character(parameters[which(parameters[,"activation.cols"]==1),"Parameter"])
    markers.names <- unique(clustering.markers.names, activation.markers.names)
    exclude.markers.names <- as.character(parameters[which(parameters[,"exclude.cols"] == 1),"Parameter"])

    # Format parameters
    parameters[, "Antigen"] <- stringr::str_replace_all(as.character(parameters[, "Antigen"]), "[²&é~{}()'[']|è`ç^à=°+$£¤%µù*><,;:!?./§-]", ".")
    clustering.markers <- as.character(parameters[which(parameters[,"clustering.cols"]==1),"Antigen"])
    activation.markers <- as.character(parameters[which(parameters[,"activation.cols"]==1),"Antigen"])
    exclude.markers <- as.character(parameters[which(parameters[,"exclude.cols"] == 1),"Antigen"])
    markers <- unique(clustering.markers, activation.markers)

  } # End Check input


  # FlowAI QC
  {
    easyFlow:::flowai(dir = dir,
                      files = files,
                      exclude.markers.names = exclude.markers.names,
                      metadata = metadata,
                      do_flowai = do_flowai,
                      do_pregating = do_pregating)
  }


  # Read FCS
  {
    if(do_flowai){
      fs <- suppressMessages(flowCore::read.flowSet(files = paste0(dir, "/results/FCS/QC/HighQuality/", gsub(".fcs","_hiQ.fcs", files)),
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

    } # end if do_flowai

    # Create results folder
    suppressWarnings(
      dir.create(paste0(dir,"/Rdata/")) +
        dir.create(paste0(dir,"/results/"))
    )

    message("\n")
    message(" ########## Preprocessing data ########## ")
    ff_data <- data.frame(flowCore::fsApply(fs, flowCore::exprs))

  } # End read FCS


  # Add ID tag
  {
    # Add an unique ID
    message("Adding ID tag...")
    ff_data <- cbind(ff_data, c(1:nrow(ff_data)))
    names(ff_data) <- c(as.character(parameters[,"Parameter"]), "ID_event") #colnames(ff)
    exclude.markers <- c(exclude.markers, "ID_event")
    exclude.markers.names <- c(exclude.markers.names, "ID_event")

    # Add ID_File tag
    ID_file <- NULL
    for(i in 1:length(files)){
      if(is.null(ID_file)){
        ID_file <- rep(i, nrow(fs[[i]]))
      } else {
        ID_file <- c(ID_file, rep(i, nrow(fs[[i]])) )
      }
    }
    ff_data <- cbind(ff_data, ID_file)
    exclude.markers <- c(exclude.markers, "ID_file")
    exclude.markers.names <- c(exclude.markers.names, "ID_file")

    # Randomizing rows
    #sampling.id <- sample(c(1:nrow(ff_data)), nrow(ff_data), replace = F)
    #sum(duplicated(sampling.id))

    # prepare data for Rtsne (matrix format required)
    data_table <- as.matrix(ff_data) #[sampling.id,])
  } # End Add ID tag


  if(toupper(data_type)=="FCM" && !sum(parameters[,"colsToTransform"])==0){
    message("Transforming data (autoLogical)...")
    #estim.lgcl <- estim.lgcl(flowFrame(data_table), channels = clustering.markers.names, m=cofactor)
    #data_table <- transform(data_table, estim.lgcl)
    auto.lgcl <- suppressWarnings(easyFlow:::autoLgcl(flowCore::flowFrame(data_table), channels = transform.markers.names, m=cofactor))
    ff_Rtsne <- flowCore::transform(flowCore::flowFrame(data_table), auto.lgcl)
    data_table <- ff_Rtsne@exprs
  } # end if FCM


  if(toupper(data_type)=="CYTOF" && !sum(parameters[,"colsToTransform"])==0){
    message("Transforming data (cytofAsinh)...")
    #data_table[,transform.markers.names] <- easyFlow:::cytofAsinh(as.matrix(data_table[,transform.markers.names]), cofactor = cofactor)
    data_table[,transform.markers.names] <- asinh(as.matrix(data_table[,transform.markers.names])/cofactor)
  } # End CyTOF

  # Rename clustering markers
  colnames(data_table)[colnames(data_table) %in% clustering.markers.names] <- clustering.markers
  colnames(data_table)[colnames(data_table) %in% activation.markers.names] <- activation.markers
  colnames(data_table)[colnames(data_table) %in% exclude.markers.names] <- exclude.markers

  # Add noise
  #data_table <- add_noise(data_table, clustering.markers)

  # Compute quantiles
  if(compute.quantiles==TRUE){
    message("Quantiles scaling...")
    # Calculate Quantiles for each marker
    rng <- matrixStats::colQuantiles(data_table[,clustering.markers], probs = c(0+quantile_thresh, 1-quantile_thresh))
    data_table[,clustering.markers] <- t((t(data_table[,clustering.markers]) - rng[, 1]) / (rng[, 2] - rng[, 1]))

    data_table[data_table[,clustering.markers] < 0] <- 0
    data_table[data_table[,clustering.markers] > 1] <- 1

  } # end compute.quantiles


  # Compute scaling
  if(compute.scaling==TRUE){
    message("Scaling data...")
    data_table[,clustering.markers] <- easyFlow:::colScale(data_table[,clustering.markers],
                                                           center = TRUE,
                                                           scale = TRUE,
                                                           add_attr = TRUE,
                                                           rows = NULL,
                                                           cols = NULL)
    summary(data_table)
  } # end compute.scaling


  {
    # Format assignments
    assignments <- data.frame(metadata)[,c("Name","bc","tp","ind")]
    row.names(assignments) <- as.character(assignments[,"Name"])
    assignments[,"Name"] <- NULL

    # Format conditions
    suppressWarnings({
      condition1.files <- gsub(".fcs", "_easyFlow", rownames(assignments[assignments[,"bc"]==condition1,]))
      condition2.files <- gsub(".fcs", "_easyFlow", rownames(assignments[assignments[,"bc"]==condition2,]))
    })


    # Format variable
    variable <- as.numeric(metadata[,"variable"])
    names(variable) <- row.names(assignments)
    variable

    # Format status
    status <- as.logical(metadata[,"status"])
    names(status) <- row.names(assignments)
    status

    lineage_markers <- clustering.markers
    functional_markers <- activation.markers

    # Format assignements
    assignments$bc <- factor(assignments$bc)

    row.names(assignments) %>%
      gsub(".fcs", "_easyFlow", . ) ->
      row.names(assignments)

    # Format variable
    names(variable) %>%
      gsub(".fcs", "_easyFlow", . ) ->
      names(variable)

    # Format status
    names(status) %>%
      gsub(".fcs", "_easyFlow", . ) ->
      names(status)

    # Set a file list
    files_extract <- paste0(dir, "/results/FCS/", gsub(".fcs$", "_easyFlow.fcs", basename(files)))
    files_extract

    basename(files_extract) %>%
      gsub(".fcs$","",.) ->
      samples
    samples

    # Select assignments
    ix.samples <- row.names(assignments) %in% samples  # samples %in% row.names(assignments)
    ix.samples

    assignments_extract <- assignments[ix.samples,]
    assignments_extract$bc <- factor(assignments$bc[ix.samples])
    assignments_extract$bc
  }
  ############# End : Section General CytofWorkFlow ##############

  {
    ###########################################################################
    # Section : Select events associated with specific biological conditions
    ###########################################################################

    ## Define colors for conditions
    color_conditions <- c("#6A3D9A", "#FF7F00", "grey70")
    names(color_conditions) <- factor(c(condition1, condition2, "total")) #levels(assignments_extract$bc)

    # Select a subset of biological conditons
    ix.bc <- ((assignments_extract$bc==condition1) | (assignments_extract$bc==condition2))
    assignments_extract <- assignments_extract[ix.bc,]

    assignments_extract$bc <- factor(assignments_extract$bc, levels = c(condition1, condition2)[!c(condition1, condition2)=="NA"])
    assignments_extract$bc

  }
  # End : section - Select events associated with biological conditions
  #####################################################################

  {
    ###########################################
    ####### Section General CytofWorkFlow #####
    ###########################################
    ## Generate sample IDs corresponding to each cell in the `expr` matrix
    sample_ids <- rep(gsub(".fcs$","_easyFlow",fs@phenoData@data$name), flowCore::fsApply(fs, nrow))
    sample_test <- rep(ix.samples, flowCore::fsApply(fs, nrow))
    sample_ids <- sample_ids[sample_test]

    ## Extract expression
    ID_event <- data_table[sample_test,"ID_event"]

  }
  ############################################################################
  # End section General
  #############################################################################



  save(list = c("dir", "FItSNE.dir", "files", "parameters", "metadata",
                "data_type","cofactor",
                "compute.quantiles","quantile_thresh", "compute.scaling",
                "data_table",
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

  message("Preprocessing done.")

  return(dir)

}
