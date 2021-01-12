#'#########################################################################################
#' Compute a quality control of each samples and extract events in the specified gate of the metadata file
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @source http://bioconductor.org/packages/release/bioc/html/flowAI.html
#' @source https://academic.oup.com/bioinformatics/article/32/16/2473/2240408
#'
#' @source https://www.bioconductor.org/packages/release/bioc/html/CytoML.html
#' @source https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.23663
#'
#' @param dir a character string specifying the project directory
#' @param files a array of character string specifying the fcs filenames
#' @param metadata a dataframe containing the metadata
#' @param exclude.markers.names a array of character string specifying the markers excluded of the analysis
#' @param do_flowai a logical specifying if the flowAI should be computed
#' @param do_pregating a logical specifying if the pregating should be computed
#'
#' @import magrittr
#' @importFrom flowAI flow_auto_qc
#' @importFrom flowCore read.FCS compensate write.FCS
#' @importFrom flowWorkspace GatingSet getData getIndices getCompensationMatrices
#' @importFrom CytoML getSamples parseWorkspace openWorkspace


flowai <- function(dir, files, metadata, exclude.markers.names, do_flowai, do_pregating){

  if(do_flowai == TRUE){

    {
      message("\n")
      message(" ########## Running FlowAI quality control ########## ")

      suppressWarnings(
        dir.create(paste0(dir, "/results/FCS/QC/report"), recursive = TRUE) +
          dir.create(paste0(dir, "/results/FCS/QC/HighQuality"), recursive = TRUE)
      )

    }


    ##############
    # QC control #
    ###############
    {

      i <- 1
      for (i in 1:length(files)){

        file <- files[i]

        message(paste0("FCS file ", i, " / ", length(files)))

        flowAI::flow_auto_qc(paste0(dir, "/", file),

                             remove_from = "all",

                             #             output = 1,

                             timeCh = NULL, #Character string corresponding to the name of the Time Channel in the set of FCS files. By default is NULL and the name is retrieved automatically.

                             second_fractionFR = .1, # The fraction of a second that is used to split the time channel in order to recreate the flow rate. Set it to "timestep" if you wish to recreate the flow rate at the maximum resolution allowed by the flow cytometry instrument.

                             alphaFR = 0.01,# The level of statistical significance used to accept anomalies detected by the ESD method. The default value is 0.01.

                             decompFR = TRUE, # Logical indicating whether the flow rate should be decomposed in the trend and cyclical components. Default is TRUE and the ESD outlier detection will be executed on the trend component penalized by the magnitude of the cyclical component. If it is FALSE the ESD outlier detection will be executed on the original flow rate.

                             ChExcludeFS = exclude.markers.names, # Character vector with the names or name patterns of the channels that you want to exclude from the signal acquisition check. The default option, c("FSC", "SSC"), excludes the scatter parameters. If you want to include all the parameters in the analysis use NULL.

                             outlier_binsFS = FALSE, #logical indicating whether outlier bins (not events) have to be removed before the changepoint detection of the signal acquisition check. The default is FALSE.

                             pen_valueFS = 200, #The value of the penalty for the changepoint detection algorithm. This can be a numeric value or text giving the formula to use; for instance, you can use the character string "1.5*log(n)", where n indicates the number of cells in the FCS file. The higher the penalty value the less strict is the detection of the anomalies. The default is 200.

                             max_cptFS = 3, # The maximum number of changepoints that can be detected for each channel. The default is 3.

                             ChExcludeFM = exclude.markers.names, #Character vector with the names or name patterns of the channels that you want to exclude from the signal acquisition check. The default option, c("FSC", "SSC"), excludes the scatter parameters. If you want to include all the parameters in the analysis use NULL.

                             sideFM = "both", # dynamic range check has to be executed on both limits, the upper limit or the lower limit

                             neg_valuesFM = 1, #Scalar indicating the method to use for the removal of the anomalies from the lower limit of the dynamic range. Use 1 to remove negative outliers or use 2 to truncate the negative values to the cut-off indicated in the FCS file.

                             html_report = "_QCreport", # paste0("report/", gsub(pattern=".fcs", replacement="",basename(file)),"_QC/"),

                             #             mini_report = paste0("report/", gsub(pattern=".fcs", replacement="",basename(file)),"_QCmini"),

                             fcs_QC = FALSE,

                             fcs_highQ =  "_hiQ", # paste0("HighQuality/", gsub(pattern=".fcs", replacement="",basename(file)),"_hiQ"),

                             fcs_lowQ = FALSE,

                             folder_results = dir
        )


        # Move html report to a new folder
        from_html.file <- paste0(dir, "/", gsub(pattern=".fcs", "_QCreport.html", file))
        to_html.file <- paste0(dir, "/results/FCS/QC/report/", gsub(pattern=".fcs", "_QCreport.html", basename(file)))

        file.rename(from = from_html.file,
                    to = to_html.file)


        # Move "hiQ" fcs file to a new folder
        from_hiQ.file <- paste0(dir, "/", gsub(pattern=".fcs", "_hiQ.fcs", file))
        to_hiQ.file <- paste0(dir, "/results/FCS/QC/HighQuality/", gsub(pattern=".fcs", "_hiQ.fcs", basename(file)))

        file.rename(from = from_hiQ.file,
                    to = to_hiQ.file)

      } # End for files

      # Set a list of files
      hiQ_files <- paste0(dir, "/results/FCS/QC/HighQuality/", gsub(pattern=".fcs", "_hiQ.fcs", basename(files)))
      hiQ_dir <- paste0(dir, "/results/FCS/QC/HighQuality/")

      # Move QCmini.txt
      file.rename(from = paste0(dir, "/QCmini.txt"),
                  to = paste0(dir, "/results/FCS/QC/QCmini.txt")
      )

    } # End FlowAI QC

  } # End do_flowAI


  #################################################################################################################
  #
  # Extract selected gate events
  #
  ###############################################################################
  # Gatinset
  # https://github.com/gfinak/OpenCytoTutorial/blob/master/OpenCytoTutorial.Rmd
  ################################################################################

  if(do_pregating == TRUE){

    message("\n")
    message(" ########## Pregating events ########## ")
    message(" ")

    files %>%
      basename(.) ->
      #gsub(".fcs$","",.)
      tube.names

    # Check all FCS files are defined in metadata table
    if( !all(metadata[,"Name"] %in% tube.names)){stop("Not all fcs files are defined in the metadata")}

    # Select corresponding gatingFiles
    metadata[,"FlowJo_workspace"] %>%
      paste0(dir,"/attachments/gating/", . ) %>%
      unique(.) ->
      gatingFiles

    selected.gates <- metadata[,"gate"]

    i <- 1
    offset <- 0
    # Read all GatingFiles
    for (i in 1:length(unique(gatingFiles))){

      ###########################
      # Import FlowJo workspace
      ###########################

      # Select a gatingFile
      gatingFile <- gatingFiles[i]
      gatingFile

      # Select files for a level of gatingFiles
      files.selected <- files[metadata[,"FlowJo_workspace"] %in% basename(gatingFile)]
      files.selected

      # Open FlowJo workspace
      ws <- CytoML::openWorkspace(gatingFile)
      ws
      summary(ws)

      # Get samples name
      samples_list <- CytoML::getSamples(ws)$name

      basename(files.selected) %>%
        gsub("_hiQ.fcs", ".fcs", .) ->
        files.selected_list

      if(!all(files.selected_list %in% samples_list) ){ #|| !all(samples_list %in% files.selected_list) ){
        message("### ERROR ###")
        message("Tube Names :")
        message(gsub(".fcs$", ",  ", files.selected_list[!files.selected_list %in% samples_list]))
        message("are absent of the FlowJo workspace:")
        message(basename(gatingFile))
        message(" ")
        message("This FlowJo workspace contains the following samples :")
        message(paste(samples_list, collapse = ",  "))
        message(" ")
        message("#########################################################################################")
        message("IMPORTANT : Tubenames and filenames must be identical for FlowJo workspace parsing step.")
        message("#########################################################################################")
        stop("Stop gating extraction")
      }


      ####################################
      # Parse workspace on hiQ fcs files (auto comp, transform and apply gates)
      ####################################
      # Parse workspace (require original fcs files)
      message(paste0("Parsing workspace ", basename(gatingFile), " (",i,"/",length(levels(factor(gatingFiles))),")"))
      suppressMessages(
        gs_tmp <- CytoML::parseWorkspace(ws,
                                         name = 1,
                                         path = dir,
                                         isNcdf = FALSE,
                                         includeGates = TRUE
                                         )
        )

      #name = 1 : The name or index of the group of samples (1=All samples)
      #isNcdf = FALSE,  logical specifying if you would like to use netcdf to store the data, or if you would like to keep all the flowFrames in memory.
      #subset=1 #numeric vector specifying the subset of samples in a group to import.
      #requiregates=TRUE
      #execute=TRUE #a logical specifying if the gates, transformations, and compensation should be immediately calculated
      #compensation=NULL
      #,keywords=

      message("Workspace parsed")

      # Explore all fcs files in the current gatingFile
      j <- 1
      for (j in 1:length(files.selected_list)){
        message(paste0("Processing ",files.selected_list[[j]]))

        # Create a gatingSet for the current fcs file
        #print("Applying GatingSet on hiQ fcs...")
        suppressMessages(
          gs <- flowWorkspace::GatingSet(x = gs_tmp[[j]], #character or flowSet or GatingHierarchy
                                         y = gsub(".fcs","_hiQ.fcs", files.selected_list[[j]]), #character or missing
                                         path = paste0(dir, "/results/FCS/QC/HighQuality/"),
                                         isNcdf = FALSE,
                                         includeGates = FALSE
          )
        )
                                         #includeGates logical whether to parse the gates or just simply extract the flowJo stats
                                         #guids, #character vectors to uniquely identify each sample (Sometime FCS file names alone may not be unique)
                                         #sampNloc = "keyword", #character scalar indicating where to get sampleName(or FCS filename) within xml workspace. It is either from "keyword" or "sampleNode".
                                         #xmlParserOption, #integer option passed to xmlTreeParse
                                         #wsType="win" #character workspace type, can be value of "win", "macII", "vX", "macIII"
                                         #                execute = FALSE, #a logical specifying if the gates, transformations, and compensation should be immediately calculated
                                         #compensation=NULL


        #############################
        # Plot the GatingHierarchy
        #plot(gs[[1]])
        # Population  statistics  for  all  populations  in  the sample
        #getPopStats(gs)
        # Get compensation matrice
        #getCompensationMatrices(gs[[1]])
        # Retrieve the samples name from gs
        #samples <- sampleNames(gs)
        # Get the nodes list
        #nodelist<-getNodes(gs, prefix="all", path = "full") # chemin complet
        #nodelist <- getNodes(gs, prefix="all", path = 1) # chemin abbr?g?
        ################################

        selected.gate <- as.character(selected.gates[which(metadata[, "Name"] == files.selected_list[[j]])])

        # Get events gated on selected gate
        fs <- flowWorkspace::getData(gs, selected.gate)
        ff <- fs[[1]] # ff contenant les events tri?s

        # Extract event index for selected gate
        selected.ids <- flowWorkspace::getIndices(gs[[1]], selected.gate)


        ###################################
        # Compensate hiQ fcs files
        ###################################
        # Extract compensation matrix
        comp.matrix <- flowWorkspace::getCompensationMatrices(gs_tmp[[1]])

        # Compensate hiQ fcs
        ff <- flowCore::read.FCS(paste0(dir, "/results/FCS/QC/HighQuality/", gsub(".fcs","_hiQ.fcs",files.selected_list[[j]])))

        if(!is.null(comp.matrix)){
          ff <- flowCore::compensate(ff, comp.matrix)
        }
        ff

        # Remove internal compensation matrix
        comp <- ff@description$SPILL
        if (!is.null(comp)){
          ff@description$SPILL <- NULL
        }
        comp <- ff@description$`$SPILLOVER`
        if (!is.null(comp)){
          ff@description$`$SPILLOVER` <- NULL
        }
        ff

        ###################################
        # Isolate selected events
        ###################################
        #print("Saving fcs with selected events only...")

        # Filter comp fcs exprs data using selected.events.index
        ff@exprs <- ff@exprs[selected.ids,]
        ff

        # Write selected fcs
        flowCore::write.FCS(ff, filename = paste0(dir, "/results/FCS/QC/HighQuality/", gsub(".fcs","_hiQ.fcs",files.selected_list[[j]])) )

      } # End do_pregating

    } # End For File in selectd.Files

  } # End For GatingFile in GatingFiles

}
