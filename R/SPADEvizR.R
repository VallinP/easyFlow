#'#########################################################################################
#' SPADEVizR QC & statistical functions
#'
#' Adapted from :
#' Guillaume Gautreau et al., Bioinformatics, March 2017
#' https://academic.oup.com/bioinformatics/article/33/5/779/2662307
#' https://github.com/tchitchek-lab/SPADEVizR
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @source https://academic.oup.com/bioinformatics/article/33/5/779/2662307
#' @source https://github.com/tchitchek-lab/SPADEVizR
#'
#' @import magrittr
#' @importFrom SPADEVizR importResultsFromFCS assignContext qcSmallClusters qcUniformClusters identifyAC identifyDAC identifyCC classifyAbundanceProfiles countViewer heatmapViewer phenoViewer boxplotViewer kineticsViewer streamgraphViewer MDSViewer distogramViewer
#'
#' @param dir

SPADEvizR <- function(dir = NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  message("\n")
  message(" ########## SPADEVizR ########## ")

  message("Loading data...")
  source(paste0(dir,"/attachments/easyFlow.R"))
  load(paste0(dir,"/Rdata/easyFlow.RData"))
  load(paste0(dir,"/Rdata/som_map.RData"))
  load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))

  ################
  # SPADEVizR QC #
  ################
  {
    results <- SPADEVizR::importResultsFromFCS(path=paste0(dir, "/results/FCS/"),
                                               exclude.markers = exclude.markers[!exclude.markers == "cluster"],
                                               clustering.markers = clustering.markers, #.names,
                                               probs = c(0.05, 0.95),
                                               trans = "none", # if FALSE, apply arcsinh transf.
                                               quantile.approximation = FALSE,
                                               th.min_cells = 0,
                                               assignments = NULL)

    results <- SPADEVizR::assignContext(results, assignments = assignments)

    nClus <- n.metaclusters
    clusters <- levels(factor(metacluster)) # as.character(1:nClus)
    pdf.width <- length(samples)/2
    if (pdf.width <= 6){pdf.width <- 6}

    suppressWarnings(
      dir.create(paste0(dir, "/results/SPADEVizR_QC/"), recursive = TRUE) +
        dir.create(paste0(dir, "/results/SPADEVizR/plots/"), recursive = TRUE)
    )

    message("\n")
    message("Generating qc small clusters")
    suppressMessages(
      small <- SPADEVizR::qcSmallClusters(results, width=pdf.width, height = (5+nClus/4), PDFfile = paste0(dir, "/results/SPADEVizR_QC/SPADEVizR-qcreport-SmallClusters_heatmap.pdf"), th = 50)
    )
    message("\n")
    message("Generating Uniform Phenotypes QC")
    suppressMessages(
      uniform <- SPADEVizR::qcUniformClusters(results,
                                              density.PDFfile = paste0(dir, "/results/SPADEVizR_QC/SPADEVizR-qcreport-UniformClusters_density.pdf"),
                                              heatmap.PDFfile = paste0(dir, "/results/SPADEVizR_QC/SPADEVizR-qcreport-UniformClusters_heatmap.pdf"),
                                              heatmap.PDFfile.dim = c(length(results@marker.names),(length(results@cluster.names)+5)),
                                              uniform.test = "unimodality")
    )

    ####################### End : QC ####################

  } #end SPADEVizR QC



  ###################
  # SPADEVizR stats #
  ###################
  {
    pdf.width <- (5 + nClus)

    resultsAC <- NULL
    resultsDAC <- NULL
    resultsCC <- NULL
    results_AP <- NULL
  }

  {
    # 4. SPADEVizR analysis
    #######################

    # 4.1 Identification of cell clusters having an abundance greater than a specific value (Abundant Clusters)
    if (inherits(try(SPADEVizR::identifyAC(results, samples = samples, mu = 1, th.pvalue = 0.01),
                     silent=TRUE), "try-error") == FALSE){
      message("\n")
      resultsAC <- SPADEVizR::identifyAC(results, samples = samples, mu = 1, th.pvalue = 0.01)
      # visualizes the identified Abundant Clusters stored in the AC object
      pdf(file = paste0(dir, "/results/SPADEVizR/plots/results_AbondantClusters.pdf"))
      SPADEVizR::plot(resultsAC)
      dev.off()
    } else {
      message("\n")
      message("SPADEVizR cannot perform report results_AbondantClusters")
    }


    # 4.2 Identification of cell clusters having an abundance different between two biological conditions (Differentially Abundant Clusters)
    if(inherits(try(SPADEVizR::identifyDAC(results, #a 'Results' object
                                           condition1 = condition1.files, #a character vector providing the sample names defined as the first condition
                                           condition2 = condition2.files, #a character vector providing the sample names defined as the second condition
                                           use.percentages = TRUE, #a logical specifying if the computations should be performed on percentage
                                           method = "t.test", #a character specifying the name of the statistical test to use "t.test" or "wilcox.test"
                                           method.adjust = NULL, #a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method)
                                           method.paired= TRUE, #use.method.paired, #a logical indicating if the statistical test must be performed in a paired manner
                                           th.pvalue = 0.05, #a numeric specifying the p-value threshold
                                           th.fc = 2),
                    silent=TRUE), "try-error") == FALSE){
      message("\n")
      resultsDAC <- SPADEVizR::identifyDAC(results, #a 'Results' object
                                           condition1 = condition1.files, #a character vector providing the sample names defined as the first condition
                                           condition2 = condition2.files, #a character vector providing the sample names defined as the second condition
                                           use.percentages = TRUE, #a logical specifying if the computations should be performed on percentage
                                           method = "t.test", #a character specifying the name of the statistical test to use "t.test" or "wilcox.test"
                                           method.adjust = NULL, #a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method)
                                           method.paired= TRUE, #use.method.paired, #a logical indicating if the statistical test must be performed in a paired manner
                                           th.pvalue = 0.05, #a numeric specifying the p-value threshold
                                           th.fc = 2) #a numeric specifying the fold-change threshold
      resultsDAC@results$significant
      # visualizes the identified Differentially Abundant Clusters stored in the DAC object
      pdf(file = paste0(dir, "/results/SPADEVizR/plots/results_DiffAbondantClusters.pdf"))
      SPADEVizR::plot(resultsDAC)
      dev.off()
    } else {
      message("SPADEVizR cannot perform report results_DiffAbondantClusters")
    }


    # 4.3 Identification of cell clusters having an abundance correlated with a biological variable (Correlated Clusters)
    # To compute the correlation, "variable" must be numerical values
    if((inherits(try(SPADEVizR::identifyCC(results,
                                           variable = variable, #a numerical named vector providing the correspondence between a sample name (in rownames) and the specific numerical phenotype
                                           use.percentages = TRUE, #a logical specifying if the computations should be performed on percentage
                                           method = "pearson", #a character indicating the correlation method to use: "pearson", "spearman"
                                           method.adjust = NULL, #a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method)
                                           th.correlation = 0.75, #a numeric specifying the absolute value of the correlation coefficient threshold
                                           th.pvalue = 0.05),
                     silent=TRUE), "try-error") == FALSE) && !anyNA(variable)){
      message("\n")
      resultsCC <- SPADEVizR::identifyCC(results,
                                         variable = variable, #a numerical named vector providing the correspondence between a sample name (in rownames) and the specific numerical phenotype
                                         use.percentages = TRUE, #a logical specifying if the computations should be performed on percentage
                                         method = "pearson", #a character indicating the correlation method to use: "pearson", "spearman"
                                         method.adjust = NULL, #a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method)
                                         th.correlation = 0.75, #a numeric specifying the absolute value of the correlation coefficient threshold
                                         th.pvalue = 0.05) #a numeric specifying the p-value threshold
      # visualizes the identified Correlated Clusters stored in the CC object
      pdf(file = paste0(dir, "/results/SPADEVizR/plots/results_CorrelatedClusters.pdf"))
      SPADEVizR::plot(resultsCC)
      dev.off()
    } else {
      resultsCC <- NULL
      message("\n")
      message("SPADEVizR cannot perform report results_CorrelatedClusters")
    }


    # 4.4 Classification of cell clusters based on theirs abundance profiles
    # performs a k-means clustering to identify 6 classes of abundance profiles
    if(inherits(try(SPADEVizR::classifyAbundanceProfiles(results,
                                                         method = "k-means", #a character specifying the clustering method among: "hierarchical_h", "hierarchical_k","k-means","eigencell","clique"
                                                         method.parameter = 6),
                    silent=TRUE), "try-error") == FALSE){
      set.seed(seed)
      message("\n")
      results_AP <- SPADEVizR::classifyAbundanceProfiles(results,
                                                         method = "k-means", #a character specifying the clustering method among: "hierarchical_h", "hierarchical_k","k-means","eigencell","clique"
                                                         method.parameter = 6) #a numeric specifying the numeric value required by the selected method
      # visualizes the classified Abundance Profiles stored in the AP object
      pdf(file = paste0(dir, "/results/SPADEVizR/plots/results_AbundanceProfile.pdf"))
      SPADEVizR::plot(results_AP)
      dev.off()
    } else {
      results_AP <- NULL
      message("SPADEVizR cannot perform report results_AbundanceProfile")
    }

  } # End SPADEVizR part1 : description


  # 5. Visualization methods
  {
    # 5.1. Visualization of the number of cells associated to each cluster (Count Viewer)
    if(inherits(try(SPADEVizR::countViewer(results, samples = samples),
                    silent=TRUE), "try-error") == FALSE){
      pdf(file = paste0(dir, "/results/SPADEVizR/plots/countview.pdf"), width = pdf.width)
      SPADEVizR::countViewer(results, samples = samples)
      dev.off()
    }

    # 5.3. Visualization of cell cluster phenotypes using a categorical heatmap (Heatmap Viewer)
    # displays an heatmap representation summarizing phenotypes for the overall dataset
    pdf(file = paste0(dir, "/results/SPADEVizR/plots/heatmapview.pdf"), width = pdf.width)
    SPADEVizR::heatmapViewer(results,
                             clusters = NULL, #a character vector providing the clusters to be used (all clusters by default=NULL)
                             markers = NULL,#a character vector providing the markers to be used (all markers by default=NULL)
                             num = 5, #a numeric value specifying the number of markers expression categories to be used
                             clustering = "both", #a character specifying which clustering must be build ("markers", "clusters", "both" or "none")
                             tile.color = "black", #a character specifying the border color of the tiles (NA to remove tile borders)
                             show.on_device = TRUE)
    dev.off()

    if(FALSE){
      # 5.4. Visualization of cell cluster phenotypes using parallels coordinates (PhenoViewer)
      # displays the detailed cluster phenotype for cluster xx, ie 10
      #phenoViewer(results, clusters = c("10"))

      pdf(file = paste0(dir, "/results/SPADEVizR/plots/PhenoView.pdf"))
      for(i in 1:results@cluster.number){
        message(paste0("Processing cluster number ", i, "/", results@cluster.number))
        SPADEVizR::phenoViewer(results,
                               samples = NULL, #a character vector providing the sample names to used (all samples by default)
                               clusters = c(as.character(i)), #a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
                               markers = NULL, #a character vector specifying the markers to be displayed assignments
                               show.mean = "both", #a character specifying if marker means expression should be displayed, possible value are among: "none", "only" or "both"
                               show.on_device = TRUE,
                               sort.markers = TRUE) #a logical specifying if the markers must be sorted by names in the representation
      }
      dev.off()
    } # end FALSE


    # 5.5. Visualization of cell cluster abundances in different biological conditions (Boxplot Viewer)
    # displays the abundance of the cluster 33 for all samples of each biological condition
    message("\n")
    message(paste0("Printing BoxPlotView..."))
    pdf(file = paste0(dir, "/results/SPADEVizR/plots/BoxPlotView.pdf"))
    for(i in 1:results@cluster.number){
      message(paste0("Processing cluster ", i, "/", results@cluster.number))
      SPADEVizR::boxplotViewer(results,
                               samples = NULL,
                               clusters = c(as.character(i)),
                               use.percentages = TRUE, #a logical specifying if the visualization must be performed on percentage
                               show.legend = TRUE,
                               show.violin = TRUE, #a logical specifying if the count distribution must be displayed
                               show.on_device = TRUE)
    } # end For
    dev.off()


    # 5.6. Visualization of cell cluster abundance kinetics in different biological conditions (Kinetics Viewer)
    # displays the abundance of the combined clusters 51 and 57 in a kinetic manner
    #kineticsViewer(results, clusters = c("2", "12"))
    if(length(unique(assignments[,"tp"]))>=2){
      pdf(file = paste0(dir, "/results/SPADEVizR/plots/KineticView.pdf"))
      SPADEVizR::kineticsViewer(results,
                                samples = NULL,
                                clusters = as.character(1:5),
                                use.percentages = TRUE, #a logical specifying if the visualization should be performed on percentage
                                show.on_device = TRUE,
                                scale_y = NULL) #a numeric value specifying the maximal value in y axsis
      dev.off()
    }


    # 5.7. Visualization of cell cluster abundance dynamics in different samples (Streamgraph Viewer)
    # displays the relative abundances of the selected set of clusters following the order previously defined
    pdf(file = paste0(dir, "/results/SPADEVizR/plots/streamgraphView.pdf"), height = pdf.width, width = length(samples)*2)
    SPADEVizR::streamgraphViewer(results, samples = samples, clusters = clusters)
    # The same could be done in a relative manner using the `use.relative = TRUE` parameter
    SPADEVizR::streamgraphViewer(results, samples = samples, clusters = clusters, use.relative = TRUE)
    dev.off()


    # 5.8. Visualization of cell cluster or sample abundance similarities (MDS Viewer)
    if(inherits(try(SPADEVizR::MDSViewer(results,
                                         space = "samples",
                                         clusters = clusters,
                                         use.percentages = TRUE, #a logical specifying if the visualization should be performed on percentage
                                         dist.method = "euclidean", #a character string containing the name of the distance measure to use
                                         show.on_device = FALSE),
                    silent=TRUE), "try-error") == FALSE){

      pdf(file = paste0(dir, "/results/SPADEVizR/plots/MDSView.pdf"))
      # displays the similarities between all samples based on the selected cluster abundances
      SPADEVizR::MDSViewer(results,
                           space = "samples",
                           clusters = clusters,
                           use.percentages = TRUE, #a logical specifying if the visualization should be performed on percentage
                           dist.method = "euclidean", #a character string containing the name of the distance measure to use
                           show.on_device = TRUE)
      # displays the similarities between all samples based on the selected cluster abundances
      SPADEVizR::MDSViewer(results,
                           space = "clusters",
                           clusters = clusters,
                           use.percentages = TRUE, #a logical specifying if the visualization should be performed on percentage
                           dist.method = "euclidean", #a character string containing the name of the distance measure to use
                           show.on_device = TRUE)
      # Note: In MDS representations, the Kruskal Stress (KS) indicates the percentage of information lost during the dimensionality reduction process
      dev.off()
    }


    # 5.9. Visualization of pairwise marker co-expressions (Distogram Viewer)
    # specifies a set of samples to use
    #samples <- c("PBD08_BB078", "PBD08_BB231", "PBD08_BC641", "PBD08_BD619", "PBD08_BD620")

    # specifies a set of clusters to use

    pdf(file = paste0(dir, "/results/SPADEVizR/plots/distogramView.pdf"))
    # displays a distogram representation showing all marker co-expressions filtered by the selected samples and clusters
    SPADEVizR::distogramViewer(results,
                               samples = samples,
                               clusters = clusters,
                               markers = NULL,
                               show.on_device = TRUE)
    dev.off()


    # Save SPADEVizR results
    suppressWarnings(dir.create(paste0(dir, "/Rdata/")))
    save(file = paste0(dir, "/Rdata/SPADEVizR_results.RData"),
         list= c("results", "resultsAC", "resultsDAC", "resultsCC", "results_AP")
    )

  } # End SPADEVizR part2 : Visualization methods


  # 6. Modeling methods

  if(!anyNA(variable) && !anyNA(status) && FALSE){

    {
      #########################################
      # All Clusters
      ##########################################

      resultsGLM <- NULL
      resultsCPHM <- NULL
      resultsRFM <- NULL


      # 6. Modeling methods on all clusters
      DAC.sig <- as.data.frame(resultsDAC@results[, c("cluster","significant")])

      if(nrow(DAC.sig) < floor(length(samples)/2 - 1) ){

        suppressWarnings(dir.create(paste0(dir, "/results/SPADEVizR/GLM/")))
        pdf(file = paste0(dir, "/results/SPADEVizR/GLM/All_Clusters.pdf"))

        # 6.1. Prediction of biological outcomes using generalized linear models

        # generates the GLM baed on a set of differentially abundant clusters
        resultsGLM <- SPADEVizR::generateGLM(results, variable = variable, clusters = as.character(DAC.sig[, "cluster"]))

        # displays a summary of the generalized linear model
        summary(resultsGLM$model)

        # displays the coefficients associated to each cluster
        SPADEVizR::plot(resultsGLM$plot.clusters)

        # displays the variable predictions for each sample based on a model constructed with the clusters
        SPADEVizR::plot(resultsGLM$plot.samples)



        # 6.2. Prediction of biological outcomes using Cox proportional hazards regression models

        # generates the Cox model based on a set of differentially abundant clusters
        resultsCPHM <- SPADEVizR::generateCPHM(results, variable = variable, status = status, clusters = as.character(DAC.sig[, "cluster"]))

        # displays a summary of the Cox model
        summary(resultsCPHM$model)

        # displays the provided survival curve
        SPADEVizR::plot(resultsCPHM$plot.provided_survival_curve)

        # displays the abundance coefficients associated to each cluster
        SPADEVizR::plot(resultsCPHM$plot.clusters)

        # displays the variable predictions for each sample based on a model constructed with the clusters
        SPADEVizR::plot(resultsCPHM$plot.samples)



        # 6.3. Prediction of biological outcomes using random forest models

        # generates the random forest model based on a set of differentially abundant clusters
        resultsRFM <- SPADEVizR::generateRFM(results, variable = variable, status = status, clusters = as.character(DAC.sig[, "cluster"]))

        # displays variable importance
        SPADEVizR::plot(resultsRFM$plot.vimp)

        # displays the minimum depth of each variable
        SPADEVizR::plot(resultsRFM$plot.minimal_depth)

        # displays the variable predictions for each sample based on a model constructed with the clusters
        SPADEVizR::plot(resultsRFM$plot.samples)


        dev.off()

        # 6.4 Save SPADEVizR results
        save(file = paste0(dir, "/Rdata/SPADEVizR_All_Clusters_results.RData"),
             list= c("resultsGLM", "resultsCPHM", "resultsRFM")
        )

      } else {

        message("\n")
        message("the number of DAC must be inferior to the number of observation (samples) minus 1")

      }# End if DAC.sig < n samples/2 -1


    } # End SPADEVizR models All clusters

    {
      #########################################
      # per cluster analysis
      #########################################

      resultsGLM <- NULL
      resultsCPHM <- NULL
      resultsRFM <- NULL

      suppressWarnings(dir.create(paste0(dir, "/results/SPADEVizR/GLM/")))

      DAC.sig <- as.data.frame(resultsDAC@results[, c("cluster","significant")])

      i <- 11
      for(i in 1:nrow(DAC.sig)){

        # select a DAC significant
        if (DAC.sig[i,"significant"]==TRUE){

          pdf(file = paste0(dir, "/results/SPADEVizR/GLM/Cluster_", i, ".pdf"))

          # 6.1. Prediction of biological outcomes using generalized linear models

          # generates the GLM baed on a set of differentially abundant clusters
          resultsGLM <- SPADEVizR::generateGLM(results, variable = variable, clusters = as.character(DAC.sig[i, "cluster"]))

          # displays a summary of the generalized linear model
          summary(resultsGLM$model)

          # displays the coefficients associated to each cluster
          SPADEVizR::plot(resultsGLM$plot.clusters)

          # displays the variable predictions for each sample based on a model constructed with the clusters
          SPADEVizR::plot(resultsGLM$plot.samples)


          # 6.2. Prediction of biological outcomes using Cox proportional hazards regression models

          # generates the Cox model based on a set of differentially abundant clusters
          resultsCPHM <- SPADEVizR::generateCPHM(results, variable = variable, status = status, clusters = as.character(DAC.sig[i, "cluster"]))

          # displays a summary of the Cox model
          summary(resultsCPHM$model)

          # displays the provided survival curve
          SPADEVizR::plot(resultsCPHM$plot.provided_survival_curve)

          # displays the abundance coefficients associated to each cluster
          SPADEVizR::plot(resultsCPHM$plot.clusters)

          # displays the variable predictions for each sample based on a model constructed with the clusters
          SPADEVizR::plot(resultsCPHM$plot.samples)


          # 6.3. Prediction of biological outcomes using random forest models

          # generates the random forest model based on a set of differentially abundant clusters
          resultsRFM <- SPADEVizR::generateRFM(results, variable = variable, status=status, clusters = as.character(DAC.sig[i, "cluster"]))

          # displays variable importance
          SPADEVizR::plot(resultsRFM$plot.vimp)

          # displays the minimum depth of each variable
          SPADEVizR::plot(resultsRFM$plot.minimal_depth)

          # displays the variable predictions for each sample based on a model constructed with the clusters
          SPADEVizR::plot(resultsRFM$plot.samples)

          dev.off()
        }

        # 6.4 Save SPADEVizR results
        save(file = paste0(dir, "/results/FlowSOM/SPADEVizR/RData/Cluster_", i, ".results.RData"),
             list= c("resultsGLM", "resultsCPHM", "resultsRFM")
        )

      }

    } # End SPADEVizR models per cluster

  } else {

    message("\n")
    message("SPADEVizR prediction models : status and/or variable contains NAs in metadata.csv")

  } # End if Status and/or variable contains NAs

}
