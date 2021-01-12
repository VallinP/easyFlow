##########################################################################################
#' Set metaclusters
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#########################################################################################
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom NbClust NbClust
#'
#' @param dir a character string specifing the project directory
#'
#' @export

set_metaclusters <- function(dir = NULL){

  {
    if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}
    elbow_point_res <- NULL

    message("\n")
    message(" ########## Running set final metaclusters ########## ")

    message("Loading data...")
    source(paste0(dir,"/attachments/easyFlow.R"))
    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/som_map.RData"))
    load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))
    if(do.dimred == TRUE && file.exists(paste0(dir,"/Rdata/tsne.RData"))){load(paste0(dir,"/Rdata/tsne.RData"))}

    suppressWarnings(dir.create(paste0(dir, "/results/set_metaclusters/density/"), recursive = TRUE) +
                       dir.create(paste0(dir, "/results/set_metaclusters/tSNE report/"), recursive = TRUE))

  }


  # From FlowSOM R package
  {
    SSE <- function (data, clustering)
    {
      if (class(clustering) != "numeric")
        clustering <- as.numeric(as.factor(clustering))
      c_wss <- 0
      for (j in seq_along(clustering)) {
        if (sum(clustering == j) > 1) {
          c_wss <- c_wss + (nrow(data[clustering == j, , drop = FALSE]) -
                              1) * sum(apply(data[clustering == j, , drop = FALSE],
                                             2, stats::var))
        }
      }
      c_wss
    }
    # From FlowSOM R package
    findElbow <- function (data)
    {
      n <- length(data)
      data <- as.data.frame(cbind(1:n, data))
      colnames(data) <- c("X", "Y")
      min_r <- Inf
      optimal <- 1
      for (i in 2:(n - 1)) {
        f1 <- stats::lm(Y ~ X, data[1:(i - 1), ])
        f2 <- stats::lm(Y ~ X, data[i:n, ])
        r <- sum(abs(c(f1$residuals, f2$residuals)))
        if (r < min_r) {
          min_r <- r
          optimal <- i
        }
      }
      optimal
    }
  } # End flowsom functions


  # Test unimodality of metaclusters
  {
    message("Checking unimodality of metaclusters...")

    unimodal.res <- NULL
    unimodal.res <- as.list(rep(TRUE, n.metaclusters))

    for(i in rev(2:n.metaclusters)){

      classif <- metacl[,i][som.res$unit.classif]

      unimodal.res[[i]] <- easyFlow:::SPADEVizR_unimodality(df = som.res$data[[1]],
                                                            cluster = classif,
                                                            clusters = NULL,
                                                            clustering.markers = clustering.markers,
                                                            uniform.test = "unimodality", #"unimodality", "spread", "both"
                                                            th.pvalue = 0.05,
                                                            th.IQR = 2,
                                                            density.PDFfile = NULL, #paste0(dir, "/results/set_metaclusters/", i, "_metaclusters_UniformClusters_density.pdf"), #"qcUniformClusters_density.pdf",
                                                            #density.PDFfile.dim = NULL, # c(17, 10),
                                                            heatmap.PDFfile = NULL, # paste0(dir, "/results/set_metaclusters/", i, "_metaclusters_UniformClusters_heatmap.pdf"), #"qcUniformClusters_heatmap.pdf"
                                                            #tile.color = "black",
                                                            verbose = TRUE)
    }

    unimodal.perc <- 0
    for(i in 2:n.metaclusters){
      unimodal.perc[i] <- unimodal.res[[i]]$perc
    }

    # Code from FlowSOM package, Metaclustering function
    res <- unimodal.perc
    #res[1] <- 0
    res <- rev(res)
    smooth <- 0.2
    for (i in 2:(n.metaclusters-1)) { res[i] <- (1 - smooth) * res[i] + (smooth/2) * res[i - 1] + (smooth/2) * res[i + 1] }

    #elbow_point <- findElbow(res)
    #elbow_point2 <- which(rev(res) == res[elbow_point]) #(n.metaclusters-1) - elbow_point

    #message("Based on elbow point of unimodal percentages")
    #message("Found ", elbow_point2, " metaclusters")
    #message("Unimodal percentage : ", round(res[elbow_point], 2), "%")

    nbclus <- NbClust::NbClust(res, method = "ward.D", min.nc = 2, max.nc = 2, index = c("silhouette"))
    silh.score <- min(which(nbclus$Best.partition == 2))
    silh.score2 <- which(rev(res) == res[max(which(nbclus$Best.partition == 1))]) #(n.metaclusters-1) - silh.score
    message("Based on sihouette score of unimodal percentages")
    message("Found ", silh.score2, " metaclusters")
    message("Unimodal percentage : ", round(res[silh.score2], 2), "%")

    png(filename = paste0(dir, "/results/set_metaclusters/Unimodality_percentage.png"), width = 480, height = 480)
    plot(1:n.metaclusters, unimodal.perc, type = "b", xlab = "Number of Clusters",
         ylab = "Percentage of unimodal metaclusters")
    abline(v = (silh.score2 + .01), col = "blue")
    dev.off()

    #n.metaclusters <- elbow_point2

  } # end HC_unimodal



  # Elbow point
  elbow_point_res <- NULL
  {

    #if(n.metaclusters > 50){max <- 50} else {max <- n.metaclusters}
    max <- n.metaclusters
    elbow_point_res <- rep(0, max)
    elbow_point_res[1] <- SSE(som.res$codes[[1]], rep(1, nrow(som.res$codes[[1]])))
    for (i in 2:max) {
      c <- SOM.metaclustering[[i]]$consensusClass
      elbow_point_res[i] <- SSE(som.res$codes[[1]], c)
    }

    smooth <- 0.2
    for (i in 2:(max - 1)) {
      elbow_point_res[i] <- (1 - smooth) * elbow_point_res[i] + (smooth/2) * elbow_point_res[i - 1] + (smooth/2) * elbow_point_res[i + 1]
    }
    elbow_point_ix <- findElbow(elbow_point_res)
    #n.metaclusters <- findElbow(elbow_point_res)

    png(filename = paste0(dir, "/results/set_metaclusters/Elbow_point.png"), width = 480, height = 480)
    plot(1:max, elbow_point_res, type = "b", xlab = "Number of Clusters",
         ylab = "Within groups sum of squares")
    abline(v = (elbow_point_ix + .01), col = "blue")
    dev.off()

  } # Elbow point


  # UniformClusters_heatmap.pdf
  {
    message("Uniform Metaclusters Heatmap report...")
    pdf(file = paste0(dir, "/results/set_metaclusters/UniformClusters_heatmap.pdf"),
        width = 17, height = 10)

    for(i in rev(2:n.metaclusters)){

      classif <- metacl[,i][som.res$unit.classif]
      #message(paste0("For ", i, " metaclusters"))

      unimodal_tmp <- easyFlow:::SPADEVizR_unimodality(df = som.res$data[[1]], #codes[[1]],
                                                       cluster = classif,
                                                       clusters = NULL,
                                                       clustering.markers = clustering.markers,
                                                       uniform.test = "unimodality", #"unimodality", "spread", "both"
                                                       th.pvalue = 0.05,
                                                       th.IQR = 2,
                                                       density.PDFfile = NULL, #paste0(dir, "/results/set_metaclusters/density_plot/", i, "_metaclusters_UniformClusters_density.pdf"), #"qcUniformClusters_density.pdf",
                                                       #density.PDFfile.dim =  c(17, 10),
                                                       heatmap.PDFfile = TRUE, #paste0(dir, "/results/set_metaclusters/density_plot/", i, "_metaclusters_UniformClusters_heatmap.pdf"), #"qcUniformClusters_heatmap.pdf"
                                                       tile.color = "black",
                                                       verbose = FALSE
      )

    } # end For 2:n.metaclusters

    dev.off()

  } # End Heatmap unimodality


  # Heatmap
  {
    color_clusters <- rep(easyFlow:::color_clusters, ceiling(length(unique(new.parameters[,"som_cluster"]))/length(easyFlow:::color_clusters) ) )

    message(" ")
    message("Metaclusters Expression Heatmap...")
    pdf(file = paste0(dir, "/results/set_metaclusters/Heatmap of metalusters.pdf"),
        width = 17, height = 10)

    for(i in rev(2:n.metaclusters)){

      classif <- metacl[,i][som.res$unit.classif]
      #message(paste0("For ", i, " metaclusters"))

      # Calculate Quantiles for each marker
      expr01 <- som.res$data[[1]]
      rng <- matrixStats::colQuantiles(som.res$data[[1]], probs = c(0.1, 0.9))
      expr01 <- t((t(som.res$data[[1]]) - rng[, 1]) / (rng[, 2] - rng[, 1]))

      expr01[expr01 < 0] <- 0
      expr01[expr01 > 1] <- 1

      # Plot clustering heatmap of manual merged metaclusters
      easyFlow:::plot_clustering_heatmap_wrapper(expr = som.res$data[[1]],
                                                 expr01 = expr01,
                                                 cell_clustering = classif,
                                                 color_clusters = color_clusters)


    } # end print unimodality report

    dev.off()

  } # End Heatmap of metaclusters


  # Metaclusters projection on tSNE
  if(do.dimred == TRUE && file.exists(paste0(dir,"/Rdata/tsne.RData"))){

    color_clusters <- rep(easyFlow:::color_clusters, ceiling(length(unique(new.parameters[,"som_cluster"]))/length(easyFlow:::color_clusters) ) )

    message(" ")
    message("tSNE projection of metaclusters...")
    for(i in rev(2:n.metaclusters)){

      png(file = paste0(dir, "/results/set_metaclusters/tSNE report/", i, " metalusters.png"),
          width = 640, height = 640)

      classif <- metacl[,i][som.res$unit.classif]
      #message(paste0("For ", i, " metaclusters"))

      #  tSNE
      colnames(tsne.coord) <- c("tsne1", "tsne2")
      classif_factor <- as.factor(classif)
      gg_tsne_cl2m <- ggplot2::ggplot(tsne.coord,  ggplot2::aes(x = tsne1, y = tsne2, color = classif_factor)) +
        ggplot2::geom_point(size = 0.8) +
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual(values = color_clusters) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))
      print(gg_tsne_cl2m)

    dev.off()

    } # End for n.metaclusters


  } # End tSNE report


  # UniformClusters_density
  {
    message(" ")
    message("Uniform Metaclusters density report...")
    for(i in rev(2:n.metaclusters)){

      classif <- metacl[,i][som.res$unit.classif]
      #message(paste0("For ", i, " metaclusters"))

      pdf(file = paste0(dir, "/results/set_metaclusters/density/", i, "_metaclusters_UniformClusters_density.pdf"),
          width = 17, height = 10)

      unimodal_tmp <- easyFlow:::SPADEVizR_unimodality(df = som.res$data[[1]], #codes[[1]],
                                                       cluster = classif,
                                                       clusters = NULL,
                                                       clustering.markers = clustering.markers,
                                                       uniform.test = "unimodality", #"unimodality", "spread", "both"
                                                       th.pvalue = 0.05,
                                                       th.IQR = 2,
                                                       density.PDFfile = TRUE, #paste0(dir, "/results/set_metaclusters/density_plot/", i, "_metaclusters_UniformClusters_density.pdf"), #"qcUniformClusters_density.pdf",
                                                       density.PDFfile.dim =  c(17, 10),
                                                       heatmap.PDFfile = NULL, #paste0(dir, "/results/set_metaclusters/heatmap/", i, "_metaclusters_UniformClusters_heatmap.pdf"), #"qcUniformClusters_heatmap.pdf"
                                                       tile.color = "black",
                                                       verbose = FALSE
      )
      dev.off()
    } # end for

  } # End Density unimodality

  # Set final n.metaclusters value
  {
    new_n.metaclusters <- 0
    while(new_n.metaclusters < 2 || new_n.metaclusters > n.metaclusters){
      suppressWarnings(
        new_n.metaclusters <- as.numeric(readline(prompt=paste0("Enter the final n.metaclusters value (between 2 & ",  n.metaclusters, ") : ")))
      )

      if(is.na(new_n.metaclusters) || new_n.metaclusters < 2 || new_n.metaclusters > n.metaclusters){
        print(paste0(new_n.metaclusters, " is not a valid value. Final n.metaclusters must be a numeric betwwen 2 and ", n.metaclusters, "."))
        new_n.metaclusters <- 0
      }

    } # End While

    # Define new metaclustering
    n.metaclusters <- new_n.metaclusters
    metacluster <- metacl[,n.metaclusters][som.res$unit.classif]

    # Save Hierarchical Consensus Clustering
    message("Saving final n.metaclusters value...")
    suppressWarnings(dir.create(paste0(dir,"/Rdata/")))
    save(list = c("n.metaclusters",
                  "SOM.metaclustering",
                  "use_elbow_point",
                  "elbow_point_res",
                  "HC_unimodal",
                  "unimodal.res",
                  "metacl",
                  "metacluster"
    ),
    file = paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"),
    compress = "gzip"
    )

  } # End Set final n.metaclusters value


}
