#'#########################################################################################
#' Function to plot MEM report (Marker modeling enrichment).
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @import magrittr
#' @importFrom matrixStats colQuantiles
#' @importFrom stats density quantile
#' @importFrom reshape2 melt
#' @importFrom plyr dlply ldply
#' @importFrom ggplot2 ggplot aes geom_point geom_density geom_line geom_ribbon facet_wrap facet_grid ggtitle theme theme_bw scale_color_manual scale_colour_gradientn scale_y_continuous scale_x_continuous element_text element_rect element_line coord_fixed guides
#' @importFrom cowplot ggdraw draw_plot get_legend
#' @importFrom MEM build.heatmaps MEM_RMSD
#'
#' @param dir a character string specifying the project directory
#'
#' @export

MEM_report <- function(dir = NULL){

  {
    if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

    message("\n")
    message(" ########## Ploting MEM reports ########## ")

    message("Loading data...")
    source(paste0(dir,"/attachments/easyFlow.R"))
    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))
    load(paste0(dir,"/Rdata/MEM.RData"))

    set.seed(seed)

    color_clusters <- rep(easyFlow:::color_clusters, ceiling(length(unique(new.parameters[,"som_cluster"]))/length(easyFlow:::color_clusters) ) )

  }

  {
    ###########################################
    ####### Section General CytofWorkFlow #####
    ###########################################

    ## Extract expression
    expr <- data_table[ID_event,] #sample_test
    expr <- as.matrix(expr[,markers])
    #colnames(expr) <- markers

    # Extract cell.clustering
    labels <- metacluster[ID_event]

    expr <- cbind(expr, labels)
    colnames(expr) <- c(markers,"cluster")
  }
  ############################################################################
  # End section General
  #############################################################################


  {
    # Plot Median expr per MEM IDs
    {
      MEM.df <- read.csv(paste0(dir, "/results/MEM/output files/enrichment score_reorganised.csv"), sep=",")

      ## Generate sample IDs corresponding to each cell in the `expr` matrix
      MEM.ids <- MEM.df[order(as.numeric(MEM.df[,"Metacluster"])),"MEM"][labels]

      ## ----plot-merker-expression-distribution ---
      ggdf <- data.frame(sample_id = MEM.ids, expr)
      ggdf <- reshape2::melt(ggdf, id.var = "sample_id",
                   value.name = "expression", variable.name = "antigen")
      mm <- match(sample_ids, row.names(assignments_extract))
      ggdf$condition <- assignments_extract$bc[mm]
      ggdf2 <- ggdf
      ggdf2$condition <- rep("total", times = nrow(ggdf2))

      # Expression plot
      message("Printing MEM report...")
      pdf.width <- 20
      pdf.height <- 20 #length(samples)
      pdf(file=paste0(dir,"/results/MEM/MEM_metacl_exprs_plot.pdf"),
          width = pdf.width, height = pdf.height)

      for (k in 1:length(unique(MEM.df[,"MEM"]))){

        i <- unique(MEM.df[,"MEM"])[k]
        mm <- match(ggdf$sample_id, i)
        ggdf3 <- rbind(ggdf[which(!is.na(ggdf$sample_id[mm])),], ggdf2)

        plot <-  ggplot2::ggplot(ggdf3, ggplot2::aes(x = expression, color = as.factor(condition))) +
          ggplot2::geom_density() +
          ggplot2::facet_wrap(~ antigen, nrow = 4, scales = "free") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                strip.text = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 5)) +
          ggplot2::scale_color_manual(values = color_conditions) +
          ggplot2::ggtitle(paste0("Metacluster ", k," _ MEM label ",i))
        print(plot)

      } # End for levels MEM.df$MEM

      dev.off()

    } # End Plot Median expr per MEM IDs


    # Plot Median expr per MEM scores
    {
      pdf.width <- 15
      pdf.height <- 5
      pdf(file=paste0(dir,"/results/MEM/MEM_score_exprs_plot.pdf"),
          width = pdf.width, height = pdf.height)

      i <- colnames(MEM.df[,c(3:(ncol(MEM.df)-1))])[1]
      for (i in colnames(MEM.df[,c(3:(ncol(MEM.df)-1))])){

        ## Generate sample IDs corresponding to each cell in the `expr` matrix
        #sample_ids <- rep(gsub(".fcs$","",fs@phenoData@data$name), fsApply(fs, nrow))
        MEM.ids <- MEM.df[order(as.numeric(MEM.df[,"Metacluster"])),i][labels]

        ## ----plot-merker-expression-distribution ---
        expr.extract <- data.frame(matrix(expr[,i], ncol=1))
        colnames(expr.extract) <- i
        ggdf <- data.frame(sample_id = MEM.ids, expr.extract) #labels
        ggdf <- reshape2::melt(ggdf, id.var = "sample_id",
                     value.name = "expression", variable.name = "antigen")
        mm <- match(sample_ids, row.names(assignments_extract))
        ggdf$condition <- assignments_extract$bc[mm]
        ggdf$sample_id <- factor(ggdf$sample_id, levels = levels(factor(ggdf[,"sample_id"]))[order(as.numeric(levels(factor(ggdf[,"sample_id"]) )))])

        if(nrow(ggdf) >= 50000){
          ggdf2 <- ggdf[sample(c(1:nrow(ggdf)), 50000),]
        } else {
          ggdf2 <- ggdf
        }
        ggdf2$sample_id <- rep("total", times = nrow(ggdf2))

        ggdf3 <- rbind(ggdf, ggdf2)

        # Expression plot
        plot1 <-  ggplot2::ggplot(ggdf3, ggplot2::aes(x = expression, color = sample_id)) +
          ggplot2::geom_density(size=1) +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                strip.text = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 5)) +
          ggplot2::scale_color_manual(values = c(color_clusters[1:(length(levels(ggdf3$sample_id))-1)], "grey30")) +
          ggplot2::ggtitle(i)

        plot2 <-  ggplot2::ggplot(ggdf3, ggplot2::aes(x = expression, color = sample_id)) +
          ggplot2::geom_density(size=1, ggplot2::aes(x=expression, y=..scaled..)) +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                strip.text = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 5)) +
          ggplot2::scale_color_manual(values = c(color_clusters[1:(length(levels(ggdf3$sample_id))-1)], "grey30")) +
          ggplot2::ggtitle(paste0(i, " scaled"))

        legend <- cowplot::get_legend(plot1)
        print(cowplot::ggdraw() +
                cowplot::draw_plot(plot1 + ggplot2::theme(legend.position = "none"), 0, .5, .7, .5) +
                cowplot::draw_plot(plot2 + ggplot2::theme(legend.position = "none"), 0, 0, .7, .5) +
                cowplot::draw_plot(legend, .7, .0, .2, 1)
        )

      } # End for each clustering markers

      dev.off()

    } # End Plot Median expr per MEM scores


    # Plot MEM scores per metacl
    {
      pdf.width <- 15
      pdf.height <- 5
      pdf(file=paste0(dir,"/results/MEM/MEM_score_expr_plot_per_metaclusters.pdf"),
          width = pdf.width, height = pdf.height)

      for (i in colnames(MEM.df[,c(3:(ncol(MEM.df)-1))])){

        ## Generate sample IDs corresponding to each cell in the `expr` matrix
        #sample_ids <- rep(gsub(".fcs$","",fs@phenoData@data$name), fsApply(fs, nrow))
        MEM.ids <- MEM.df[order(as.numeric(MEM.df[,"Metacluster"])),i][labels]

        ## ----plot-merker-expression-distribution ---
        expr.extract <- data.frame(matrix(expr[,i], ncol=1))
        colnames(expr.extract) <- i
        ggdf <- data.frame(sample_id = sample_ids, expr.extract) #labels
        ggdf <- reshape2::melt(ggdf, id.var = "sample_id",
                     value.name = "expression", variable.name = "antigen")
        mm <- match(sample_ids, row.names(assignments_extract))
        ggdf$condition <- assignments_extract$bc[mm]
        ggdf$MEMscore <- factor(MEM.ids, levels = levels(factor(MEM.ids))[order(as.numeric(levels(factor(MEM.ids))))])
        ggdf$metacluster <- factor(labels)

        if(nrow(ggdf) >= 50000){
          ggdf2 <- ggdf[sample(c(1:nrow(ggdf)), 50000),]
        } else {
          ggdf2 <- ggdf
        }
        ggdf2$MEMscore <- rep("total", times = nrow(ggdf2))
        ggdf2$metacluster <- rep("total", times = nrow(ggdf2))

        for (j in levels(factor(MEM.ids))[order(as.numeric(levels(factor(MEM.ids))))]){

          ggdf %>%
            .[which(.$MEMscore == j),] %>%
            .[order(.$metacluster),] %>%
            rbind(., ggdf2) ->
            ggdf3

          # Expression plot
          plot1 <-  ggplot2::ggplot(ggdf3, ggplot2::aes(x = expression, color = metacluster)) +
            ggplot2::geom_density(size=1) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                  strip.text = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 5)) +
            ggplot2::scale_color_manual(values = c(color_clusters[1:(length(levels(factor(ggdf3$metacluster)))-1)], "grey30")) +
            ggplot2::ggtitle(paste0("Marqueur ",i, " , MEM score ", j))

          plot2 <-  ggplot2::ggplot(ggdf3, ggplot2::aes(x = expression, color = metacluster)) +
            ggplot2::geom_density(size=1, ggplot2::aes(x=expression, y=..scaled..)) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                  strip.text = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 5)) +
            ggplot2::scale_color_manual(values = c(color_clusters[1:(length(levels(factor(ggdf3$metacluster)))-1)], "grey30")) +
            ggplot2::ggtitle(paste0("Marqueur ",i, " , MEM score ", j, " scaled"))

          legend <- cowplot::get_legend(plot1)
          print(cowplot::ggdraw() +
                  cowplot::draw_plot(plot1 + ggplot2::theme(legend.position = "none"), 0, .5, .7, .5) +
                  cowplot::draw_plot(plot2 + ggplot2::theme(legend.position = "none"), 0, 0, .7, .5) +
                  cowplot::draw_plot(legend, .7, .0, .2, 1)
          )

        } # End for each

      } # End for each clustering markers

      dev.off()

    } # End MEM scores per metacl


    # Plot MEM scores (overlay)
    {
      pdf.width <- 15
      pdf.height <- 15
      pdf(file=paste0(dir,"/results/MEM/MEM_score_expr_plot_overlay.pdf"),
          width = pdf.width, height = pdf.height)

      i <- 1
      for (i in colnames(MEM.df[,c(3:(ncol(MEM.df)-1))])){

        ## Generate sample IDs corresponding to each cell in the `expr` matrix
        #sample_ids <- rep(gsub(".fcs$","",fs@phenoData@data$name), fsApply(fs, nrow))
        MEM.ids <- MEM.df[order(as.numeric(MEM.df[,"Metacluster"])),i][labels]

        ## ----plot-merker-expression-distribution ---
        expr.extract <- data.frame(matrix(expr[,i], ncol=1))
        colnames(expr.extract) <- i
        ggdf <- data.frame(sample_id = sample_ids, expr.extract) #labels
        ggdf <- reshape2::melt(ggdf, id.var = "sample_id",
                     value.name = "expression", variable.name = "antigen")
        mm <- match(sample_ids, row.names(assignments_extract))
        ggdf$condition <- assignments_extract$bc[mm]
        ggdf$MEMscore <- factor(MEM.ids, levels = levels(factor(MEM.ids))[order(as.numeric(levels(factor(MEM.ids))))])
        ggdf$metacluster <- factor(labels)

        if(nrow(ggdf) >= 50000){
          ggdf2 <- ggdf[sample(c(1:nrow(ggdf)), 50000),]
        } else {
          ggdf2 <- ggdf
        }
        ggdf2$MEMscore <- rep("total", times = nrow(ggdf2))
        ggdf2$metacluster <- rep("total", times = nrow(ggdf2))

        ggdf %>%
          .[order(.$MEMscore),] %>%
          rbind(., ggdf2) ->
          ggdf3

        #calculate the density depending on V2
        res <- dlply(ggdf3, .(MEMscore), function(x) stats::density(x$expression))
        dd <- plyr::ldply(res, function(z){
          data.frame(Values = z[["x"]],
                     V1_density = z[["y"]],
                     V1_count = z[["y"]]*z[["n"]])
        })

        #add an offset depending on V2
        d <- stats::density(ggdf3$expression)
        offset <- (max(d$y)/2)
        dd$offest=-as.numeric(dd$MEMscore)*offset # adapt the 0.2 value as you need # now, rely on max density
        dd$V1_density_offest=dd$V1_density+dd$offest

        #and plot
        range <- stats::quantile(ggdf3$expression,c(.01,.99))

        # Expression plot
        plot1 <- ggplot2::ggplot(dd, ggplot2::aes(Values, V1_density_offest, color=MEMscore)) +
          ggplot2::geom_line()+
          ggplot2::geom_ribbon(ggplot2::aes(Values, ymin=offest,ymax=V1_density_offest, fill=MEMscore),alpha=0.3)+
          ggplot2::scale_y_continuous(breaks=NULL) +
          ggplot2::scale_x_continuous(breaks=c(floor(range[1]-1.5):ceiling(range[2]+1.5)),limits=c(range[1]-1.5,range[2]+1.5)) +
          ggplot2::ggtitle(i)

        print(plot1)

      } # End for each clustering markers

      dev.off()

    } # End Plot MEM scores (overlay)

  } # End Plot


  ## Create heatmap of results. For detailed explanation of the arguments, see ?build.heatmaps
  {
    pdf(file = paste0(dir,"/results/MEM/MEM.pdf"), width = 20, height = 10)
    suppressWarnings(MEM::build.heatmaps(MEM_values,
                                         cluster.MEM = "both",
                                         cluster.medians = "both",
                                         cluster.IQRs = "both",
                                         display.thresh = 0,
                                         newWindow.heatmaps=F,
                                         output.files = F))

    MEM::MEM_RMSD(MEM_vals = MEM_values, format=NULL, newWindow.heatmaps = FALSE, output.matrix = FALSE)
    dev.off()
    # PDFs of median and MEM heatmaps will be written to /output files folder in your working directory
  }

} # End plot MEM report
