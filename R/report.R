#'#########################################################################################
#' Report function
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_density facet_wrap facet_grid ggtitle theme theme_bw scale_color_manual scale_colour_gradientn element_rect element_line coord_fixed guides guide_legend
#' @importFrom colorRamps matlab.like2
#' @importFrom matrixStats colQuantiles
#' @importFrom kohonen getCodes add.cluster.boundaries
#'
#' @param dir a character string specifying the project directory
#'
#' @export

report <- function(dir=NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  {
    message("\n")
    message(" ########## Generating pdf reports ########## ")

    # Load data
    message("Loading data...")
    source(paste0(dir,"/attachments/easyFlow.R"), chdir=T)
    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/som_map.RData"))
    load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))

    data_table <- cbind(data_table, new.parameters)
    data_table[,"som_cluster"] <- as.factor(data_table[,"som_cluster"])
    data_table[,"cluster"] <- as.factor(data_table[,"cluster"])

    # order data_table
    ix <- data_table[,"cluster"]
    data_table <- data_table[order(as.numeric(ix)),]

    # Calculate Quantiles for each marker
    expr <- as.matrix(data_table[ , markers])
    rng <- matrixStats::colQuantiles(expr, probs = c(0.01, 0.99))
    expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    expr01[expr01 < 0] <- 0
    expr01[expr01 > 1] <- 1

    color_clusters <- rep(easyFlow:::color_clusters, ceiling(length(unique(new.parameters[,"som_cluster"]))/length(easyFlow:::color_clusters) ) )

  }


  ##########################
  # Metaclustering report
  ##########################
  {
    message("clustering report...")
    expr01 <- cbind(expr01, data_table[, c("som_cluster", "cluster", "ID_file")])

    suppressWarnings(dir.create(paste0(dir, "/results/SOM/"), recursive = TRUE))
    pdf(file = paste0(dir, "/results/SOM/",
                      n.metaclusters," Metaclusters - ", flowSOM_dim*flowSOM_dim, " clusters.pdf"),
        bg = "transparent",
        width = 7, height = 7)

    # Distance entre noeuds utilsant les codebooks
    #d <- dist(som.res$codes[[1]])

    # cah
    # Le poid des noeuds est ignoré
    #cah <- hclust(d, method = "ward.D")
    #plot(cah, hang=-1)

    # Cut the tree
    #tree.groups <- cutree(cah, n.metaclusters)
    #tree.groups <- metacluster

    # Plot the tree w/ groups
    #plot(cah, hang=-1)
    #rect.hclust(cah, k=n.metaclusters, border = color_clusters)

    # Plot the som w/ HC groups
    #plot(som.res, type="mapping", bgcol = color_clusters[tree.groups])

    kohonen:::plot.kohonen(som.res, type="mapping", col = as.integer(),
                           pchs = as.integer(), bgcol = color_clusters[as.integer(metacluster)],
                           main = "mapping plot", shape = "straight", border = NA)
    kohonen::add.cluster.boundaries(som.res, clustering = metacluster)

    #plot(x, type = c("codes", "changes", "counts", "dist.neighbours", "mapping", "property", "quality")

    kohonen:::plot.kohonen(som.res, type = "codes")

    #plot(som.res, type = "changes")

    kohonen:::plot.kohonen(som.res, type="counts", shape = "straight")

    kohonen:::plot.kohonen(som.res, type = "quality", palette.name = terrain.colors, shape = "straight")

    kohonen:::plot.kohonen(som.res, type = "dist.neighbours", palette.name = terrain.colors, shape = "straight")
    #plot(som.res, type = "mapping")
    #plot(som.res, type = "property")

    i <- 1
    for(i in 1:length(clustering.markers)){
      kohonen:::plot.kohonen(som.res, type = "property", property = kohonen::getCodes(som.res, 1)[,i],
                             main = colnames(kohonen::getCodes(som.res, 1))[i], shape = "straight")
      kohonen::add.cluster.boundaries(som.res, clustering = metacluster)
    }

    dev.off()
  } # end metaclustering report


  ###################
  ### CREATE PLOT ###
  ###################

  if(do.dimred == TRUE){

    message("tSNE report...")

    #load(paste0(dir,"/Rdata/tsne.RData"))
    expr01 <- cbind(expr01, data_table[, c("som_cluster", "cluster", "ID_file", "tsne_1", "tsne_2")])

    suppressWarnings(rm(list= c("data_table.extract", "res", "dd", "plot1")))
    data_table.extract <- expr01[, c("som_cluster", "cluster", "ID_file", "tsne_1", "tsne_2")]
    colnames(data_table.extract) <- c("som_cluster", "cluster", "ID_file", "tsne_1", "tsne_2")
    #data_table.extract[,"ID_file"] <- as.factor(data_table.extract[,"ID_file"])
    #data_table.extract[,"som_cluster"] <- as.factor(data_table.extract[,"som_cluster"])
    #data_table.extract[,"cluster"] <- as.factor(data_table.extract[,"cluster"])
    head(data_table.extract)
    tail(data_table.extract)

    # plot 2-dimensional t-SNE projection with FlowSOM clustering
    plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2, color = as.factor(cluster))) +
      ggplot2::geom_point(size = 1) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::ggtitle(paste0(n.metaclusters," Metaclusters - Init ", init.matrix,
                              " - tSNE prediction ", tsne_prediction,
                              " - Perplexity ", paste(perplexity.value, collapse = "-"),
                              " - Freedom ", freedom.value)) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = color_clusters) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4), ncol = 2))

    suppressWarnings(dir.create(paste0(dir, "/results/FItSNE/"), recursive = TRUE))
    png(filename = paste0(dir, "/results/FItSNE/",
                          n.metaclusters," Metaclusters - Init ", init.matrix,
                          " - tSNE prediction ", tsne_prediction,
                          " - Perplexity ", paste(perplexity.value, collapse = "-"),
                          " - Freedom ", freedom.value, ".png"),
        bg = "transparent",
        width = 2048, height = 2048)
    print(plot1)
    dev.off()


    png(file=paste0(dir,"/results/FItSNE/tSNE.png"),
        bg = "transparent",
        width = 2048, height = 2048)

    plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2)) +
      ggplot2::geom_point(size = 1) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::ggtitle(paste0("FitSNE - Init ", init.matrix,
                              " - tSNE prediction ", tsne_prediction,
                              " - Perplexity ", paste(perplexity.value, collapse = "-"),
                              " - Freedom ", freedom.value)) +
      ggplot2::theme_bw()

    print(plot1)

    dev.off()


    png(file=paste0(dir,"/results/FItSNE/tSNE_per_File.png"),
        bg = "transparent",
        height = 2048,
        width = (1600*length(files)+512))

    plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2)) +
      ggplot2::geom_point(size = 0.3) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::ggtitle(paste0("FitSNE - Init ", init.matrix,
                              " - tSNE prediction ", tsne_prediction,
                              " - Perplexity ", paste(perplexity.value, collapse = "-"),
                              " - Freedom ", freedom.value)) +
      ggplot2::theme_bw() +
      ggplot2::facet_grid( . ~ ID_file)

    print(plot1)

    dev.off()


    i <- 1
    #colours <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100)
    colours <- colorRamps::matlab.like2(11)
    suppressWarnings(dir.create(paste0(dir,"/results/FItSNE/Markers report/")))
    for (i in 1:length(markers)){

      png(file=paste0(dir,"/results/FItSNE/Markers report/Marker_", markers[i], " .png"),
          bg = "transparent",
          width = 2048, height = 2048)

      suppressWarnings(rm(list= c("data_table.extract", "res", "dd", "plot1")))
      data_table.extract <- data_table.extract <- expr01[, c(markers[i], "ID_file", "tsne_1", "tsne_2")] #cbind.data.frame(data_table[,markers[i]], data_table[,"ID_file"], data_table[,"tsne_1"], data_table[,"tsne_2"])
      colnames(data_table.extract) <- c("marker.level", "ID_file", "tsne_1", "tsne_2")

      # order data_table.extract
      ix <- data_table.extract[,"marker.level"]
      data_table.extract <- data_table.extract[order(as.numeric(ix)),]

      plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2, color = marker.level)) +
        ggplot2::geom_point(size = 1) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::ggtitle(paste0(markers.names[i])) +
        #  geom_hex(bins = 128)  +
        ggplot2::scale_colour_gradientn(colours = colours) +
        #theme_bw() +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "grey95",
                                                   colour = "grey95",
                                                   size = 0.5, linetype = "solid"),
          panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'solid',
                                                   colour = "grey95"),
          panel.grid.minor = ggplot2::element_line(size = 0.25, linetype = 'solid',
                                                   colour = "grey95")
        )

      print(plot1)

      dev.off()
    } # end Markers


    suppressWarnings(dir.create(paste0(dir,"/results/FItSNE/Markers per file/")))
    i <- 1
    for (i in 1:length(markers)){

      png(file=paste0(dir,"/results/FItSNE/Markers per file/Markers_", markers[i], ".png"),
          height = 2048,
          width = (1600*length(files)+512))

      suppressWarnings(rm(list= c("data_table.extract", "res", "dd", "plot1")))
      data_table.extract <- expr01[, c(markers[i], "ID_file", "tsne_1", "tsne_2")] #cbind.data.frame(data_table[,markers[i]], data_table[,"ID_file"], data_table[,"tsne_1"], data_table[,"tsne_2"])
      colnames(data_table.extract) <- c("marker.level", "ID_file", "tsne_1", "tsne_2")

      # order data_table.extract
      ix <- data_table.extract[,"marker.level"]
      data_table.extract <- data_table.extract[order(as.numeric(ix)),]

      plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2, color = marker.level)) +
        ggplot2::geom_point(size = 0.1) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::ggtitle(paste0(markers.names[i])) +
        #  geom_hex(bins = 128)  +
        ggplot2::scale_colour_gradientn(colours = colours) +
        #theme_bw() +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "grey95",
                                                   colour = "grey95",
                                                   size = 0.5, linetype = "solid"),
          panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'solid',
                                                   colour = "grey95"),
          panel.grid.minor = ggplot2::element_line(size = 0.25, linetype = 'solid',
                                                   colour = "grey95")
        ) +
        ggplot2::facet_grid( . ~ ID_file)


      print(plot1)
      dev.off()
    } # end Markers per file

    # Per metaclusters
    suppressWarnings(dir.create(paste0(dir,"/results/FItSNE/metaclusters report/"), recursive = TRUE))
    png(file=paste0(dir,"/results/FItSNE/metaclusters report/per_metaclusters.png"),
        height = 512,
        width = 512 * length(unique(metacluster)))

    suppressWarnings(rm(list= c("data_table.extract", "res", "dd", "plot1")))
    data_table.extract <- data_table[, c("cluster", "ID_file", "tsne_1", "tsne_2")] #cbind.data.frame(data_table[,"cluster"], data_table[,"ID_file"], data_table[,"tsne_1"], data_table[,"tsne_2"])
    #colnames(data_table.extract) <- c("cluster", "ID_file", "tsne_1", "tsne_2")

    # order data_table.extract
    ix <- data_table.extract[,"cluster"]
    data_table.extract <- data_table.extract[order(as.numeric(ix)),]

    plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2)) +
      ggplot2::geom_point(size = .1) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::ggtitle(paste0("metaclusters report")) +
      #  geom_hex(bins = 128)  +
      # scale_colour_gradientn(colours = colours) +
      ggplot2::theme_bw() +
      ggplot2::facet_grid( . ~ cluster)

    print(plot1)
    dev.off()

    if(FALSE){
      # Per clusters
      png(file=paste0(dir,"/results/FItSNE/metaclusters report/per_clusters.png"),
          height = 128,
          width = 128 * length(unique(som.cluster)))

      suppressWarnings(rm(list= c("data_table.extract", "res", "dd", "plot1")))
      data_table.extract <- data_table[, c("som_cluster", "ID_file", "tsne_1", "tsne_2")] #cbind.data.frame(som.cluster, data_table[,"ID_file"], data_table[,"tsne_1"], data_table[,"tsne_2"])
      colnames(data_table.extract) <- c("cluster", "ID_file", "tsne_1", "tsne_2")

      # order data_table.extract
      ix <- data_table.extract[,"cluster"]
      data_table.extract <- data_table.extract[order(as.numeric(ix)),]

      plot1 <- ggplot2::ggplot(data_table.extract, ggplot2::aes(x = tsne_1, y = tsne_2)) +
        ggplot2::geom_point(size = .3) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::ggtitle(paste0("clusters report")) +
        #  geom_hex(bins = 128)  +
        # scale_colour_gradientn(colours = colours) +
        ggplot2::theme_bw() +
        ggplot2::facet_grid( . ~ cluster)

      print(plot1)
      dev.off()
    }

  } # end tsne report

} # end report
