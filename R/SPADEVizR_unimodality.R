#'#########################################################################################
#' Test unimodality of metaclusters using SPADEVizR function
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'
#' Adapted from :
#' Guillaume Gautreau et al., Bioinformatics, March 2017
#' https://academic.oup.com/bioinformatics/article/33/5/779/2662307
#' https://github.com/tchitchek-lab/SPADEVizR
#'########################################################################################
#'
#' @source https://academic.oup.com/bioinformatics/article/33/5/779/2662307
#' @source https://github.com/tchitchek-lab/SPADEVizR
#'
#' @importFrom flowCore flowSet flowFrame fsApply each_col
#' @importFrom stats quantile
#' @importFrom diptest dip.test
#' @importFrom ggplot2 ggplot aes aes_string arrow geom_point geom_segment geom_vline geom_hline geom_density geom_line geom_ribbon facet_wrap facet_grid ggtitle theme theme_bw scale_fill_manual scale_color_manual scale_colour_gradientn scale_y_continuous scale_x_continuous element_blank element_text element_rect element_line coord_fixed guides position_jitter geom_boxplot stat_summary position_jitterdodge position_dodge unit as_labeller xlab
#'
#' @param df a matrix of numerical values
#' @param clustering.markers an array of character string
#' @param cluster a null value
#' @param clusters a null value
#' @param uniform.test a character string
#' @param th.pvalue a numerical value
#' @param density.PDFfile a character string
#' @param density.PDFfile.dim an array of 2 numerical values
#' @param heatmap.PDFfile a character string
#' @param tile.color
#' @param verbose a logical
#'
#' @return uniform.res

SPADEVizR_unimodality <- function (df = NULL,
                                   cluster = NULL,
                                   clusters = NULL,
                                   clustering.markers = NULL,
                                   uniform.test = "unimodality", #"unimodality", "spread", "both"
                                   th.pvalue = 0.05,
                                   th.IQR = 2,
                                   density.PDFfile = NULL, #"qcUniformClusters_density.pdf",
                                   density.PDFfile.dim = c(8, 5),
                                   heatmap.PDFfile = NULL, #"qcUniformClusters_heatmap.pdf"
                                   tile.color = "black",
                                   verbose = TRUE)
{

  df <- cbind(df[, clustering.markers], cluster)
  clusters <- sort(unique(cluster))
  flowset <- flowCore::flowSet(flowCore::flowFrame(df))

  markersSPADEVizR <- clustering.markers

  accuracy.matrix <- matrix(nrow = length(clusters), ncol = length(markersSPADEVizR),
                            dimnames = list(clusters, markersSPADEVizR))
  min <- floor(min(flowCore::fsApply(flowset[, flowset@colnames !=
                                               "cluster"], flowCore::each_col, base::min, na.rm = TRUE),
                   na.rm = TRUE))
  max <- ceiling(max(flowCore::fsApply(flowset[, flowset@colnames != "cluster"],
                                       flowCore::each_col, base::max, na.rm = TRUE),
                     na.rm = TRUE))
  ordered.markers <- c(gtools::mixedsort(clustering.markers),
                       gtools::mixedsort(setdiff(markersSPADEVizR, clustering.markers)))
  bold.markers <- ifelse(is.element(ordered.markers, clustering.markers),
                         "bold", "plain")
  colored.markers <- ifelse(is.element(ordered.markers, clustering.markers),
                            "blue", "black")
  count <- 0
  for (cluster in clusters) {
    data.facet <- data.frame(ind = c(), subtitle = c(),
                             upper = c(), lower = c(), median = c(), pinnacle = c())
    if (verbose) {
      count <- count + 1
      #message(paste0("Cluster: ", count, " on ", length(clusters)))
    }
    #expressions.list <- vector("list", length(flowset))
    #for (sample in 1:length(flowset)) {
    #  frame <- flowset[[sample]]@exprs
    #  expressions.list[[sample]] <- frame[frame[, "cluster"] ==
    #                                        cluster, colnames(frame) %in% ordered.markers]
    #}
    expressions <- df#do.call(rbind, expressions.list)


    for (marker in ordered.markers) {
      if (nrow(expressions) > 1) {
        marker.expression <- expressions[, marker]
        subtitle <- ""

        if (uniform.test == "unimodality" || uniform.test == "both") {
          suppressMessages(p.value <- suppressMessages(diptest::dip.test(marker.expression)$p.value))
          if (p.value < th.pvalue) {
            uniform <- FALSE
            subtitle <- paste0(subtitle, " Non-unimodale distribution detected (p.value = ",
                               round(p.value, 2), ")")
          }
          else {
            uniform <- TRUE
            subtitle <- paste0(subtitle, " Unimodale distribution detected (p.value = ",
                               round(p.value, 2), ")")
          }
          facet <- data.frame(ind = marker, subtitle = subtitle,
                              upper = NA, lower = NA, median = NA, pinnacle = NA)
        }
        else {
          uniform <- TRUE
        }

        if (uniform.test == "spread" || uniform.test ==
            "both") {
          quantile <- stats::quantile(marker.expression)
          IQR <- quantile[4] - quantile[2]
          pinnacle <- SPADEVizR:::computemode(marker.expression)$y
          if (IQR < th.IQR) {
            uniform <- ifelse(uniform, TRUE, FALSE)
            subtitle <- paste0(subtitle, "\n", " Low spread detected (IQR < ",
                               th.IQR, ")")
          }
          else {
            uniform <- FALSE
            subtitle <- paste0(subtitle, "\n", " High spread detected (IQR > ",
                               th.IQR, ")")
          }
          facet <- data.frame(ind = marker, subtitle = subtitle,
                              upper = quantile[4], lower = quantile[2],
                              median = quantile[3], pinnacle = pinnacle)
        }
        data.facet <- rbind(data.facet, facet)
        accuracy.matrix[cluster, marker] <- uniform
      }
      else {
        accuracy.matrix[cluster, marker] <- NA
      }
    } # marker in ordered.markers


    if (!is.null(density.PDFfile)) {
      data.expressions <- as.data.frame(expressions)
      data.expressions <- utils::stack(data.expressions)
      data.expressions$ind <- factor(data.expressions$ind,
                                     levels = ordered.markers)
      subtitle <- paste(levels(data.expressions$ind),
                        "\n", data.facet$subtitle)
      names(subtitle) <- levels(data.expressions$ind)
      title <- paste0("Densities of marker expressions for cluster : ",
                      cluster, " (", format(nrow(expressions), big.mark = " "),
                      " cells)")
      fill <- ifelse(accuracy.matrix[cluster, ], "green",
                     "red")
      names(fill) <- colnames(accuracy.matrix)
      data.facet$cm <- data.facet$ind %in% clustering.markers
      data.facet$ind <- factor(data.facet$ind, levels = ordered.markers)
      plot <- ggplot2::ggplot() + ggplot2::ggtitle(title) +
        ggplot2::geom_density(data = data.expressions,
                              ggplot2::aes_string(x = "values", fill = "ind"),
                              alpha = 0.1, show.legend = FALSE) +
        ggplot2::scale_fill_manual(name = "uniform", values = fill) +
        ggplot2::scale_x_continuous(limits = c(min, max), expand = c(0.01, 0)) +
        ggplot2::scale_y_continuous(expand = c(0.02, 0))

      if (uniform.test == "spread" || uniform.test == "both") {
        data.facet$IQR <- round(data.facet$upper - data.facet$lower,
                                2)
        plot <- plot + ggplot2::geom_vline(data = data.facet,
                                           ggplot2::aes_string(xintercept = "median"),
                                           color = "blue", size = 0.2) + ggplot2::geom_vline(data = data.facet,
                                                                                             ggplot2::aes_string(xintercept = "upper"),
                                                                                             color = "blue", linetype = "dashed", size = 0.1) +
          ggplot2::geom_vline(data = data.facet, ggplot2::aes_string(xintercept = "lower"),
                              color = "blue", linetype = "dashed", size = 0.1) +
          ggplot2::geom_segment(data = data.facet, ggplot2::aes_string(x = "median",
                                                                       y = "pinnacle*0.75", xend = "upper", yend = "pinnacle*0.75"),
                                color = "blue", size = 0.1, arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"))) +
          ggplot2::geom_segment(data = data.facet,
                                ggplot2::aes_string(x = "median", y = "pinnacle*0.75",
                                                    xend = "lower", yend = "pinnacle*0.75"), color = "blue", size = 0.1,
                                arrow = ggplot2::arrow(length = ggplot2::unit(0.1,"cm"))) +
          ggplot2::geom_label(data = data.facet,
                              ggplot2::aes_string(x = "median", y = "pinnacle*0.75", label = "IQR"), color = "blue", size = 1,
                              label.padding = ggplot2::unit(0.1, "lines"))
      }
      plot <- plot + ggplot2::geom_rect(data = data.facet,
                                        ggplot2::aes_string(color = "cm"), fill = NA,
                                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                                        size = 1.5, show.legend = FALSE) + ggplot2::scale_color_manual(values = c("grey80", "blue")) +
        ggplot2::facet_wrap(~ind, scales = "free", labeller = ggplot2::as_labeller(subtitle)) +
        ggplot2::xlab("marker expressions") + ggplot2::theme_bw() +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 5))
      plot(plot)
    } # End density plot

  } # End cluster in clusters

  quality.matrix.cm <- accuracy.matrix[, colnames(accuracy.matrix) %in%
                                         clustering.markers]
  uniform.cluster <- apply(quality.matrix.cm, 1, all)
  perc <- length(uniform.cluster[uniform.cluster == TRUE])/length(uniform.cluster) * 100
  accuracy.matrix <- data.frame(accuracy.matrix, uniform.cluster = uniform.cluster,
                                check.names = FALSE)
  accuracy.matrix$cluster <- rownames(accuracy.matrix)


  if (!is.null(heatmap.PDFfile)) {
    accuracy.matrix.melted <- reshape2::melt(accuracy.matrix,
                                             id = "cluster")
    colnames(accuracy.matrix.melted) <- c("cluster", "marker",
                                          "uniform")
    ordered.markers <- c(ordered.markers, "uniform.cluster")
    bold.markers <- ifelse(is.element(ordered.markers, clustering.markers),
                           "bold", "plain")
    colored.markers <- ifelse(is.element(ordered.markers,
                                         clustering.markers), "blue", "black")
    accuracy.matrix.melted$cluster <- factor(accuracy.matrix.melted$cluster,
                                             levels = rev(clusters))
    accuracy.matrix.melted$marker <- factor(accuracy.matrix.melted$marker,
                                            levels = ordered.markers, ordered = TRUE)
    details <- ""
    if (uniform.test == "unimodality" || uniform.test ==
        "both") {
      details <- paste0(details, "having unimodal distribution (p.value < ",
                        th.pvalue, ")")
    }
    if (uniform.test == "both") {
      details <- paste0(details, " and ")
    }
    if (uniform.test == "spread" || uniform.test == "both") {
      details <- paste0(details, "having low spread (th.IQR < ",
                        th.IQR, ")")
    }
    title <- paste0("QC: Uniform Phenotypes")
    subtitle <- paste0("percentage of clusters having only uniform clustering marker distributions = ",
                       format(round(perc, 2), nsmall = 2), "% (marker ",
                       details, ")")
    plot <- ggplot2::ggplot(data = accuracy.matrix.melted) +
      ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)),
                                                  "")))) + ggplot2::geom_tile(ggplot2::aes_string(x = "marker",
                                                                                                  y = "cluster", fill = "uniform"), colour = tile.color) +
      ggplot2::scale_fill_manual(values = c(`FALSE` = "red",
                                            `TRUE` = "green"), na.value = "grey50") + ggplot2::scale_x_discrete(expand = c(0,
                                                                                                                           0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::xlab("markers") + ggplot2::ylab("clusters") +
      ggplot2::geom_vline(xintercept = (ncol(accuracy.matrix) -
                                          1.5), colour = "black", size = 2) + ggplot2::geom_vline(xintercept = (length(clustering.markers) +
                                                                                                                  0.5), colour = "grey40", size = 1.5) +
      ggplot2::geom_vline(xintercept = (ncol(accuracy.matrix) - 0.5), colour = "black", size = 1) + ggplot2::geom_vline(xintercept = 0.5,
                                                                                                                        colour = "black", size = 1) +
      ggplot2::geom_hline(yintercept = 0.5, colour = "black", size = 1) +
      ggplot2::geom_hline(yintercept = length(clusters) + 0.5, colour = "black", size = 1) + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1, vjust = 0.5, face = bold.markers,
                                                         color = colored.markers), legend.key = ggplot2::element_blank(),
                     axis.line = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank())
    plot(plot)
  } # End heatmap

  uniform.res <- list(perc = perc, accuracy.matrix = accuracy.matrix)
  return(uniform.res)

}
