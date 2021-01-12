#'#########################################################################################
#' CytofWorkFlow QC & statistical functions
#'
#' Adapted from :
#' Nowicka M, Crowell H, Robinson M (2019).
#' cytofWorkflow: CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets.
#' https://f1000research.com/articles/6-748
#' https://github.com/markrobinsonuzh/cytofWorkflow.
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @source https://f1000research.com/articles/6-748
#' @source https://github.com/markrobinsonuzh/cytofWorkflow
#'
#' @import magrittr
#' @importFrom stringr str_replace_all
#' @importFrom matrixStats colQuantiles
#' @importFrom Matrix colSums rowSums
#' @importFrom reshape2 melt dcast
#' @importFrom limma plotMDS
#' @import plyr
#' @importFrom dplyr group_by summarize_all funs
#' @importFrom ggplot2 ggplot aes geom_point geom_density geom_line geom_ribbon facet_wrap facet_grid ggtitle theme theme_bw scale_color_manual scale_colour_gradientn scale_y_continuous scale_x_continuous element_blank element_text element_rect element_line coord_fixed guides position_jitter geom_boxplot stat_summary position_jitterdodge position_dodge
#' @importFrom ggrepel geom_label_repel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom flowCore fsApply flowSet flowFrame exprs
#' @importFrom grid gpar
#'
#' @param dir


CytofWorkFlow <- function(dir = NULL){

  #################
  # CytofWorkFlow #
  #################

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  message("\n")
  message(" ########## CytofWorkFlow ########## ")

  message("Loading data...")
  source(paste0(dir,"/attachments/easyFlow.R"))
  load(paste0(dir,"/Rdata/easyFlow.RData"))
  load(paste0(dir,"/Rdata/som_map.RData"))
  load(paste0(dir,"/Rdata/hierarch.consensus.clustering.RData"))

  set.seed(seed)

  color_clusters <- rep(easyFlow:::color_clusters, ceiling(length(unique(new.parameters[,"som_cluster"]))/length(easyFlow:::color_clusters) ) )

  {
    ###########################################
    ####### Section General CytofWorkFlow #####
    ###########################################

    ## Extract expression
    expr <- data_table[sample_test,]
    expr <- as.matrix(expr[,markers])
    colnames(expr) <- markers

    # Extract cell.clustering
    cell_clustering2m <- metacluster[ID_event] #cell_clustering1m[ID_event]

    # Calculate Quantiles for each marker
    rng <- matrixStats::colQuantiles(expr, probs = c(0.01, 0.99))
    expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    expr01[expr01 < 0] <- 0
    expr01[expr01 > 1] <- 1


    ############################################################################
    # End section General
    #############################################################################


    #############################################################################
    # CytofWorkFlow
    #############################################################################

    ## Create contrasts
    contrast_names <- c(paste0(condition2, "_vs_", condition1))

    suppressWarnings(
      dir.create(paste0(dir, "/results/CytofWorkflow/",contrast_names, "/plots/global/"), recursive = T) +
        dir.create(paste0(dir, "/results/CytofWorkflow/",contrast_names, "/plots/DAC/"), recursive = T) +
        dir.create(paste0(dir, "/results/CytofWorkflow/",contrast_names, "/tables/DAC/"), recursive = T) +
        dir.create(paste0(dir, "/results/CytofWorkflow/",contrast_names, "/plots/DEM/"), recursive = T) +
        dir.create(paste0(dir, "/results/CytofWorkflow/",contrast_names, "/tables/DEM/"), recursive = T)
    )


    ########################################
    #
    ### Statistical models parameters
    #
    ########################################

    ## Model formula without random effects
    model.matrix( ~ bc, data = assignments_extract)

    ## ----diff-FDR-cutoff-----------------------------------------------------
    FDR_cutoff <- 0.05

    ### Differential Abondant Clusters : Binomial regression
    ## ----diff-formula-glmer-binomial-----------------------------------------
    formula_glmer_binomial1 <- y/total ~ condition + (1|sample_id)
    formula_glmer_binomial2 <- y/total ~ condition + (1|patient_id) + (1|sample_id)

    ### Differential Expression of Functionnal and Clustering markers : Linear regression model
    ## ----diff-expr2-formula--------------------------------------------------
    formula_lm <- y ~ condition
    formula_lmer <- y ~ condition + (1|patient_id)


    ##########################
    # Plots
    ##########################

    message("\n")
    message("Printing global plots...")
    #############################
    ####### Section General #####
    #############################
    {
      ## ----plot-merker-expression-distribution ---
      ggdf <- data.frame(sample_id = sample_ids, expr)
      ggdf <- reshape2::melt(ggdf, id.var = "sample_id",
                             value.name = "expression", variable.name = "antigen")
      mm <- match(ggdf$sample_id, row.names(assignments_extract))
      ggdf$condition <- assignments_extract$bc[mm]

      # Expression plot
      pdf.width <- 20
      pdf.height <- length(samples)
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/global/exprs_plot.pdf"),
          width = pdf.width, height = pdf.height)
      plot <-  ggplot2::ggplot(ggdf, ggplot2::aes(x = expression, color = condition,
                                                  group = sample_id)) +
        ggplot2::geom_density() +
        ggplot2::facet_wrap(~ antigen, nrow = 4, scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                       strip.text = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 5)) +
        ggplot2::scale_color_manual(values = color_conditions)
      print(plot)
      dev.off()

      ## ----plot-number-of-cells
      # Barplot showing the number of cells measured for each sample in the dataset.
      # Bars are colored by experimental condition
      # Numbers in the names on the x-axis indicate patient IDs.
      # Numbers on top of the bars indicate the cell counts."----
      cell_table <- table(sample_ids)
      ggdf <- data.frame(sample_id = names(cell_table),
                         cell_counts = as.numeric(cell_table))
      mm <- match(ggdf$sample_id, row.names(assignments_extract))
      ggdf$condition <- assignments_extract$bc[mm] ################################# nb de level a corriger

      # Count plot
      pdf.width <- length(samples)
      pdf.height <- 8
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/global/count_plot.pdf"), width = pdf.width, height = pdf.height)
      plot <-  ggplot2::ggplot(ggdf, ggplot2::aes(x = sample_id, y = cell_counts, fill = condition)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_text(ggplot2::aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::scale_fill_manual(values = color_conditions, drop = FALSE) +
        ggplot2::scale_x_discrete(drop = FALSE)
      print(plot)
      dev.off()

      ## ----plot-mds, fig.cap = "MDS plot for the unstimulated (Ref) and stimulated with BCR/FcR-XL (BCRXL) samples obtained for each of the 8 healthy donors in the PBMC dataset. Calculations are based on the median (arcsinh-transformed) marker expression of 10 lineage markers and 14 functional markers across all cells measured for each sample.  Distances between samples in the plot approximate the typical change in medians. Numbers in the label names indicate patient IDs."----
      # Get the median marker expression per sample
      expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
        dplyr::group_by(sample_id) %>% dplyr::summarize_all(dplyr::funs(median))

      expr_median_sample <- t(expr_median_sample_tbl[, -1])
      colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

      mds <- limma::plotMDS(expr_median_sample, plot = FALSE)

      ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                         sample_id = colnames(expr_median_sample))
      mm <- match(ggdf$sample_id, row.names(assignments_extract))
      ggdf$condition <- assignments_extract$bc[mm]

      # MSD plot
      pdf.width <- length(samples) * 1.2
      pdf.height <- length(samples) * 1.2
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/global/MDS_plot.pdf"),
          width = pdf.width, height = pdf.height)
      plot <- ggplot2::ggplot(ggdf, ggplot2::aes(x = MDS1, y = MDS2, color = condition)) +
        ggplot2::geom_point(size = 2, alpha = 0.8) +
        ggrepel::geom_label_repel(ggplot2::aes(label = sample_id)) +
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual(values = color_conditions) +
        ggplot2::coord_fixed()
      print(plot)
      dev.off()

      ## ----plot-dendogram, fig.cap = "Heatmap of the median (arcsinh-transformed) marker expression of 10 lineage markers and 14 functional markers across all cells measured for each sample in the PBMC dataset. Color-coded with yellow for lower expression and blue for higher expression. The numbers in the heatmap represent the actual expression values. Dendrograms present clustering of samples (columns) and markers (rows) which is based on hierarchical clustering with Euclidean distance metric and average linkage. The two conditions: unstimulated (Ref) and stimulated with BCR/FcR-XL (BCRXL) for each of the 8 healthy donors are presented with a bar colored by experimental condition on top of the heatmap. Numbers in the column label names indicate patient IDs."----
      # Column annotation for the heatmap
      mm <- match(colnames(expr_median_sample), row.names(assignments_extract))
      annotation_col <- data.frame(condition = assignments_extract$bc[mm],
                                   row.names = colnames(expr_median_sample))
      annotation_colors <- list(condition = color_conditions)

      # Colors for the heatmap  ########################################## Normalize color range / exclude exclude.markers
      color <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))(10)) #(100)

      # Heatmap plot
      pdf.width <- length(samples) * 1.2
      pdf.height <- length(markers) * 1.2
      #pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/global/heatmap_plot.pdf"), width = pdf.width, height = pdf.height)
      pheatmap::pheatmap(expr_median_sample, color = color, display_numbers = TRUE,
                         number_color = "black", fontsize_number = 5, annotation_col = annotation_col,
                         annotation_colors = annotation_colors, clustering_method = "average",
                         filename = paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/global/heatmap_plot.pdf"),
                         width = pdf.width, height = pdf.height)
      #print(plot)
      #dev.off()

      ## Calculate the score
      #nrs_sample <- flowCore::fsApply(fs[, clustering.markers.names], easyFlow:::NRS, use.exprs = TRUE)
      nrs_sample <- matrix(data = NA, nrow = length(unique(data_table[, "ID_file"])), ncol = length(clustering.markers))
      for(i in unique(data_table[, "ID_file"])){
        nrs_sample[i,] <- easyFlow:::NRS(data_table[data_table[,"ID_file"] == i, clustering.markers])
      }
      ## Rename parameters names by antigen in nrs_sample tables
      colnames(nrs_sample) <- as.character(clustering.markers)

      rownames(nrs_sample) <- row.names(assignments_extract)
      nrs <- colMeans(nrs_sample, na.rm = TRUE)

      ## Plot the NRS for ordered markers
      names(sort(nrs, decreasing = TRUE)) %>%
        gsub("-",".",.) %>%
        gsub(" ",".",.) ->
        lineage_markers_ord

      nrs_sample <- data.frame(nrs_sample)
      nrs_sample$sample_id <- rownames(nrs_sample)

      ggdf <- reshape2::melt(nrs_sample, id.var = "sample_id",
                             value.name = "nrs", variable.name = "antigen")

      ggdf$antigen <- factor(ggdf$antigen, levels = lineage_markers_ord)
      mm <- match(ggdf$sample_id, row.names(assignments_extract))
      ggdf$condition <- assignments_extract$bc[mm]

      # NRS plot
      pdf.width <- length(markers) * 1.2
      pdf.height <- 7
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/global/NRS_plot.pdf"),
          width = pdf.width, height = pdf.height)
      plot <-  ggplot2::ggplot(ggdf, ggplot2::aes(x = antigen, y = nrs)) +
        ggplot2::geom_point(ggplot2::aes(color = condition), alpha = 0.9,
                            position = ggplot2::position_jitter(width = 0.3, height = 0)) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::scale_color_manual(values = color_conditions)
      print(plot)
      dev.off()

    }
    ############################################################################
    # End section General
    #############################################################################




    ###############################
    ####### Section Diff freq #####
    ###############################

    message("\n")
    message("Computing differentially abondant clusters...")

    ########################################
    #
    # Diff abondance of clusters
    #
    ########################################
    {

      ## Model formula without random effects
      #model.matrix( ~ bc, data = assignments_extract)

      k1 <- c(0, 1)
      K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
      K

      ## ----diff-FDR-cutoff-----------------------------------------------------
      #FDR_cutoff <- 0.05

      ## ----diff-freqs----------------------------------------------------------
      counts_table <- table(as.numeric(cell_clustering2m), sample_ids)
      props_table <- t(t(counts_table) / Matrix::colSums(counts_table)) * 100

      counts <- as.data.frame.matrix(counts_table)
      props <- as.data.frame.matrix(props_table)

      ## ----diff-freqs-plot-props-barplot, fig.cap = "Relative abundance of the 8 PBMC populations in each sample (x-axis), in the PBMC dataset, represented with a barplot. The 8 cell populations are a result of manual merging of the 20 FlowSOM metaclusters."----
      ## ----plot-merker-expression-distribution ---
      ggdf <- data.frame(sample_id = sample_ids, expr)
      ggdf <- reshape2::melt(ggdf, id.var = "sample_id",
                             value.name = "expression", variable.name = "antigen")
      mm <- match(ggdf$sample_id, row.names(assignments_extract))
      ggdf$condition <- assignments_extract$bc[mm]

      ####### Problème dans les colnames si caractères spéciaux présents des samples ids
      ggdf <- reshape2::melt(data.frame(cluster = rownames(props), props),
                             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")

      levels_clusters_merged <- levels(ggdf$cluster) ################################################################################# A revoir dapr?s la variable d'origine, cf. dans code des graphs pour concensus ?
      ggdf$cluster <- factor(ggdf$cluster, levels = levels_clusters_merged)

      ## Add condition info
      tube.names <- row.names(assignments) %>% gsub(" ",".",.) %>% gsub("-",".",.)
      tube.names <- stringr::str_replace_all(tube.names, "[²&é~{}()'[']|è`ç^à=°+$£¤%µù*><,;:!?./§]", ".")
      levels(ggdf$sample_id)[1]
      tube.names[1]

      mm <- match(ggdf$sample_id, tube.names)
      #mm <- match(ggdf$sample_id, row.names(assignments_extract)) ####################################################### A revoir : caract?res g?nants dans nom des tubes (espaces, et "-")
      ggdf$condition <- factor(assignments$bc[mm])

      ## ----diff-freqs-plot-props-boxplot, fig.height = 4, fig.cap = "Relative abundance of the 8 PBMC populations in each sample, in the PBMC dataset, represented with boxplots. Values for the two conditions are indicated with different colors: violet for the unstimulated (Ref) and orange for the stimulated with BCR/FcR-XL (BCRXL) samples. Values for each patient are indicated with different shape. The 8 cell populations are a result of manual merging of the 20 FlowSOM metaclusters."----
      pdf.width <- (length(levels(factor(cell_clustering2m)))*length(samples)/10)
      pdf.height <- length(markers)
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DAC/Differential.abondance.clusters_barplot.pdf"),
          height = pdf.height, width = pdf.width)

      ggdf$patient_id <- factor(assignments$ind[mm])

      plot <- ggplot2::ggplot(ggdf) +
        ggplot2::geom_boxplot(ggplot2::aes(x = condition, y = proportion, color = condition,
                                           fill = condition),  position = ggplot2::position_dodge(), alpha = 0.5,
                              outlier.color = NA) +
        ggplot2::geom_point(ggplot2::aes(x = condition, y = proportion, color = condition,
                                         shape = patient_id), alpha = 0.8, position = ggplot2::position_jitterdodge()) +
        ggplot2::facet_wrap(~ cluster, scales = "free", nrow = 2) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 6)) +
        ggplot2::scale_color_manual(values = color_conditions) +
        ggplot2::scale_fill_manual(values = color_conditions) +
        ggplot2::scale_shape_manual(values = rep(c(0:20),times=1000))
      print(plot)
      dev.off()

      pdf.width <- length(samples)*.6+1
      pdf.height <- length(markers)
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DAC/Differential.abondance.clusters_hist.pdf"),
          height = pdf.height, width = pdf.width)
      plot <- ggplot2::ggplot(ggdf, ggplot2::aes(x = sample_id, y = proportion, fill = cluster)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::facet_wrap(~ condition, scales = "free_x") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::scale_fill_manual(values = color_clusters)
      print(plot)
      dev.off()


      ###############################
      #
      # Set model formula
      #
      ###############################

      ## ----diff-formula-glmer-binomial-----------------------------------------

      #formula_glmer_binomial1 <- y/total ~ condition + (1|sample_id)
      #formula_glmer_binomial2 <- y/total ~ condition + (1|patient_id) + (1|sample_id)

      #################################
      # Format assignments_extract
      ###################################
      assignments_extract$file_name <- row.names(assignments_extract)
      assignments_extract$sample_id <- row.names(assignments_extract)
      assignments_extract$patient_id <- assignments_extract$ind
      assignments_extract$condition <- assignments_extract$bc
      assignments_extract$timepoint <- assignments_extract$tp
      assignments_extract

      if(suppressMessages(inherits(try(easyFlow:::differential_abundance_wrapper(counts, md = assignments_extract,
                                                                                 formula = formula_glmer_binomial1, K = K)
                                       , silent=TRUE), "try-error") == FALSE)) {
        message("with formula_glmer_binomial1...")
        DA_model1 <- TRUE
        ## ----diff-freqs-fit-model-1-----------------------------------------------
        da_out1 <- easyFlow:::differential_abundance_wrapper(counts, md = assignments_extract,
                                                             formula = formula_glmer_binomial1, K = K)
        apply(da_out1$adjp < FDR_cutoff, 2, table)

        ## ----diff-freqs-fit-model-output-1----------------------------------------
        da_output1 <- data.frame(cluster = rownames(props),  props,
                                 da_out1$pvals, da_out1$adjp, row.names = NULL)
        write.csv(x=data.frame(da_output1), file = paste0(dir,"/results/CytofWorkflow/",contrast_names,"/tables/DAC/Differential.abondance.clusters_Binomial.model1.csv"), row.names = TRUE)
      } else {
        DA_model1 <- FALSE
      } # End try 1


      if(suppressMessages(inherits(try(easyFlow:::differential_abundance_wrapper(counts, md = assignments_extract,
                                                                                 formula = formula_glmer_binomial2, K = K)
                                       , silent=TRUE), "try-error")==FALSE)) {
        message("with formula_glmer_binomial2...")
        DA_model2 <- TRUE
        ## ----diff-freqs-fit-model-2-----------------------------------------------
        suppressMessages(da_out2 <- easyFlow:::differential_abundance_wrapper(counts, md = assignments_extract,
                                                                              formula = formula_glmer_binomial2, K = K))
        apply(da_out2$adjp < FDR_cutoff, 2, table)
        ## ----diff-freqs-fit-model-output-2----------------------------------------
        da_output2 <- data.frame(cluster = rownames(props),  props,
                                 da_out2$pvals, da_out2$adjp, row.names = NULL)
        #print(head(da_output2), digits = 2)
        write.csv(x=data.frame(da_output2), file = paste0(dir,"/results/CytofWorkflow/",contrast_names,"/tables/DAC/Differential.abondance.clusters_Binomial.model2.csv"), row.names = TRUE)
      } else {
        DA_model2 <- FALSE
      } # End try 2


      ## ----diff-freqs-asin-sqrt-transformation-2--------------------------------
      ## Apply the arcsine-square-root transformation to the proportions
      asin_table <- asin(sqrt((t(t(counts_table) / Matrix::colSums(counts_table)))))
      asin <- as.data.frame.matrix(asin_table)
      asin_norm <- NULL

      if(DA_model2 == TRUE) {
        ## Get significant clusters and sort them by significance
        sign_clusters <- names(which(sort(da_out2$adjp[, paste0("adjp_", contrast_names)]) < FDR_cutoff))
        ## Get the adjusted p-values for the significant clusters
        sign_adjp <- da_out2$adjp[sign_clusters , paste0("adjp_", contrast_names), drop=FALSE]
        ## Normalize the transformed proportions to mean = 0 and sd = 1
        asin_norm <- easyFlow:::normalization_wrapper(asin[sign_clusters, ])
      }

      if(DA_model1 == TRUE && DA_model2 == FALSE) {
        ## Get significant clusters and sort them by significance
        sign_clusters <- names(which(sort(da_out1$adjp[, paste0("adjp_", contrast_names)]) < FDR_cutoff))
        ## Get the adjusted p-values for the significant clusters
        sign_adjp <- da_out1$adjp[sign_clusters , paste0("adjp_", contrast_names), drop=FALSE]
        ## Normalize the transformed proportions to mean = 0 and sd = 1
        asin_norm <- easyFlow:::normalization_wrapper(asin[sign_clusters, ])
      }

      ## ----diff-freqs-plot-heatmap-with-significant-clusters, fig.height = 3, fig.cap = "Normalized proportions of PBMC cell populations that are significantly differentially abundant between BCR/FcR-XL stimulated and unstimulated condition. The heat represents arcsine-square-root transformed cell frequencies that were subsequently normalized per cluster (rows) to mean of zero and standard deviation of one. The color of the heat varies from blue indicating relative under-representation to orange indicating relative over-representation. Bar at the top of the heatmap indicates the condition the samples (columns) belong to: violet for the unstimulated (Ref) and orange for the stimulated with BCR/FcR-XL (BCRXL) condition. Numbers in the brackets next to the cluster names indicate adjusted p-values. Shown are only the significant clusters for which adjusted p-values < 0.05. Clusters are sorted according to the significance so that a cluster on the top shows the most significant abundance changes between the two conditions."----
      if(!is.null(asin_norm)){
        mm <- match(colnames(asin_norm), assignments_extract$sample_id)
        pdf.width <- length(samples)*.6+1
        pdf.height <- length(sign_clusters)
        pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DAC/Differential.abondance.clusters_heatmap.pdf"),
            height = pdf.height, width = pdf.width)
        easyFlow:::plot_differential_heatmap_wrapper(expr_norm = asin_norm, sign_adjp = sign_adjp,
                                                     condition = assignments_extract$condition[mm], color_conditions = color_conditions)
        dev.off()
      }
    } # End : differential abondance of clusters


    {
      message("\n")
      message("Computing differentially expressed markers...")
      message("Per metaclusters...")
      ## ----diff-expr2-median-expression----------------------------------------
      ## Get median marker expression per sample and cluster
      expr_median_sample_cluster_tbl <- data.frame(expr[, functional_markers],
                                                   sample_id = sample_ids, cluster = as.numeric(cell_clustering2m)) %>%
        dplyr::group_by(sample_id, cluster) %>%
        dplyr::summarize_all(dplyr::funs(median))
      ## Melt
      expr_median_sample_cluster_melt <- reshape2::melt(expr_median_sample_cluster_tbl,
                                                        id.vars = c("sample_id", "cluster"), value.name = "median_expression",
                                                        variable.name = "antigen")
      ## Rearange so the rows represent clusters and markers
      expr_median_sample_cluster <- reshape2::dcast(expr_median_sample_cluster_melt,
                                                    cluster + antigen ~ sample_id,  value.var = "median_expression")
      rownames(expr_median_sample_cluster) <- paste0(expr_median_sample_cluster$cluster,
                                                     "_", expr_median_sample_cluster$antigen)
      ## Eliminate clusters with low frequency
      clusters_keep <- names(which((Matrix::rowSums(counts < 5) == 0)))
      keepLF <- expr_median_sample_cluster$cluster %in% clusters_keep
      expr_median_sample_cluster <- expr_median_sample_cluster[keepLF, ]
      ## Eliminate cases with zero expression in all samples
      keep0 <- Matrix::rowSums(expr_median_sample_cluster[, assignments_extract$sample_id]) > 0
      expr_median_sample_cluster <- expr_median_sample_cluster[keep0, ]

      ## ----diff-expr2-plot-median-expr, fig.height = 8, fig.cap = "Median (arcsinh-transformed) expression of 14 signaling markers (x-axis) across the 8 identified PBMC cell populations (individual panels). Values for the two conditions are indicated with different colors: violet for the unstimulated (Ref) and orange for the stimulated with BCR/FcR-XL (BCRXL) samples. Values for each patient are indicated with different shape. The 8 cell populations are a result of manual merging of the 20 FlowSOM metaclusters."----
      ggdf <- expr_median_sample_cluster_melt[expr_median_sample_cluster_melt$cluster
                                              %in% clusters_keep, ]
      ## Add info about samples
      mm <- match(ggdf$sample_id, assignments_extract$sample_id)
      ggdf$condition <- factor(assignments_extract$condition[mm])
      ggdf$patient_id <- factor(assignments_extract$patient_id[mm])

      pdf.width <- 2*(length(activation.markers))
      pdf.height <- 2*length(levels(factor(clusters_keep)))
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DEM/Differential.expression_activation.markers.pdf"),
          height = pdf.height, width = pdf.width)
      plot <- ggplot2::ggplot(ggdf) +
        ggplot2::geom_boxplot(ggplot2::aes(x = antigen, y = median_expression,
                                           color = condition, fill = condition),
                              position = ggplot2::position_dodge(), alpha = 0.5, outlier.color = NA) +
        ggplot2::geom_point(ggplot2::aes(x = antigen, y = median_expression, color = condition,
                                         shape = patient_id), alpha = 0.8, position = ggplot2::position_jitterdodge(),
                            size = 0.7) +
        ggplot2::facet_wrap(~ cluster, scales = "free_y", ncol=2) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::scale_color_manual(values = color_conditions) +
        ggplot2::scale_fill_manual(values = color_conditions) +
        ggplot2::scale_shape_manual(values = rep(c(0:20),times=1000)) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 2)))
      print(plot)
      dev.off()


      ## ----diff-expr2-formula--------------------------------------------------
      #formula_lm <- y ~ condition
      #formula_lmer <- y ~ condition + (1|patient_id)


      if(suppressMessages(inherits(try(easyFlow:::differential_expression_wrapper(expr_median = expr_median_sample_cluster,
                                                                                  md = assignments_extract, model = "lm", formula = formula_lm, K = K)
                                       , silent=TRUE), "try-error")==FALSE)) {

        message("with formula_lm...")
        DE_model1 <- TRUE
        ## ----diff-expr2-fit-model-1-----------------------------------------------
        de_out1 <- easyFlow:::differential_expression_wrapper(expr_median = expr_median_sample_cluster,
                                                              md = assignments_extract, model = "lm", formula = formula_lm, K = K)
        apply(de_out1$adjp < FDR_cutoff, 2, table)
        ## ----diff-expr2-fit-model-output-1----------------------------------------
        de_output1 <- data.frame(expr_median_sample_cluster,
                                 de_out1$pvals, de_out1$adjp, row.names = NULL)
        write.csv(x=data.frame(de_output1), file = paste0(dir,"/results/CytofWorkflow/",contrast_names,"/tables/DEM/Differential.expression_activation.markers_Linear.regression.model1.csv"), row.names = TRUE)

        ## ----diff-expr2-plot-heatmap-with-significant-markers, fig.height = 8, fig.cap = "Normalized expression of signaling markers in the 8 PBMC populations that are significantly differentially expressed between BCR/FcR-XL stimulated and unstimulated condition. The heat represents median (arcsinh-transformed) marker expression that was subsequently normalized per cluster-marker (rows) to mean of zero and standard deviation of one. The color of the heat varies from blue representing relative under-expression to orange representing relative over-expression. Bar at the top of the heatmap indicates the condition the samples (columns) belong to: violet for the unstimulated (Ref) and orange for the stimulated with BCR/FcR-XL (BCRXL) condition. Numbers in the brackets next to the cluster-marker names indicate adjusted p-values and cluster-marker are sorted so that they block per cluster and within each block, markers on the top show the most significant changes between the two conditions. Shown are only the significant cluster-markers for which adjusted p-values < 0.05."----
        ## Keep the significant markers, sort them by significance and group by cluster
        sign_clusters_markers <- names(which(de_out1$adjp[, paste0("adjp_", contrast_names)] < FDR_cutoff))
        oo <- order(expr_median_sample_cluster[sign_clusters_markers, "cluster"],
                    de_out1$adjp[sign_clusters_markers, paste0("adjp_", contrast_names)])
        sign_clusters_markers <- sign_clusters_markers[oo]

        ## Get the significant adjusted p-values
        sign_adjp <- de_out1$adjp[sign_clusters_markers , paste0("adjp_", contrast_names)]

        ## Normalize expression to mean = 0 and sd = 1
        expr_s <- expr_median_sample_cluster[sign_clusters_markers,assignments_extract$sample_id]
        expr_median_sample_cluster_norm <- easyFlow:::normalization_wrapper(expr_s)

        mm <- match(colnames(expr_median_sample_cluster_norm), assignments_extract$sample_id)

        if(!sum(de_out1$adjp[, paste0("adjp_", contrast_names)] < FDR_cutoff) == 0){

          # if expr_median_sample_cluster_norm not NULL
          pdf.width <- length(samples)*.6+1
          pdf.height <- length(markers)
          pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DEM/Differential.expression_activation.markers_heatmap_lm.pdf"),
              height = pdf.height, width = pdf.width)
          easyFlow:::plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_cluster_norm,
                                                       sign_adjp = sign_adjp, condition = assignments_extract$condition[mm],
                                                       color_conditions = color_conditions)
          dev.off()
        } # End de_out1 sig not NULL

      } # End try 3

      if(suppressMessages(inherits(try(easyFlow:::differential_expression_wrapper(expr_median = expr_median_sample_cluster,
                                                                                  md = assignments_extract, model = "lmer", formula = formula_lmer, K = K)
                                       , silent=TRUE), "try-error")==FALSE)) {
        message("with formula_lmer...")
        DE_model2 <- TRUE
        ## ----diff-expr2-fit-model-2-----------------------------------------------
        de_out2 <- suppressMessages(easyFlow:::differential_expression_wrapper(expr_median = expr_median_sample_cluster,
                                                                               md = assignments_extract, model = "lmer", formula = formula_lmer, K = K))
        apply(de_out2$adjp < FDR_cutoff, 2, table)
        ## ----diff-expr2-fit-model-output-2----------------------------------------
        de_output2 <- data.frame(expr_median_sample_cluster,
                                 de_out2$pvals, de_out2$adjp, row.names = NULL)
        #print(head(de_output2), digits = 2)
        write.csv(x=data.frame(de_output2), file = paste0(dir,"/results/CytofWorkflow/",contrast_names,"/tables/DEM/Differential.expression_activation.markers_Linear.regression.model2.csv"), row.names = TRUE)

        ## ----diff-expr2-plot-heatmap-with-significant-markers, fig.height = 8, fig.cap = "Normalized expression of signaling markers in the 8 PBMC populations that are significantly differentially expressed between BCR/FcR-XL stimulated and unstimulated condition. The heat represents median (arcsinh-transformed) marker expression that was subsequently normalized per cluster-marker (rows) to mean of zero and standard deviation of one. The color of the heat varies from blue representing relative under-expression to orange representing relative over-expression. Bar at the top of the heatmap indicates the condition the samples (columns) belong to: violet for the unstimulated (Ref) and orange for the stimulated with BCR/FcR-XL (BCRXL) condition. Numbers in the brackets next to the cluster-marker names indicate adjusted p-values and cluster-marker are sorted so that they block per cluster and within each block, markers on the top show the most significant changes between the two conditions. Shown are only the significant cluster-markers for which adjusted p-values < 0.05."----
        ## Keep the significant markers, sort them by significance and group by cluster
        sign_clusters_markers <- names(which(de_out2$adjp[, paste0("adjp_", contrast_names)] < FDR_cutoff))
        oo <- order(expr_median_sample_cluster[sign_clusters_markers, "cluster"],
                    de_out2$adjp[sign_clusters_markers, paste0("adjp_", contrast_names)])
        sign_clusters_markers <- sign_clusters_markers[oo]

        ## Get the significant adjusted p-values
        sign_adjp <- de_out2$adjp[sign_clusters_markers , paste0("adjp_", contrast_names)]

        ## Normalize expression to mean = 0 and sd = 1
        expr_s <- expr_median_sample_cluster[sign_clusters_markers,assignments_extract$sample_id]
        expr_median_sample_cluster_norm <- easyFlow:::normalization_wrapper(expr_s)

        mm <- match(colnames(expr_median_sample_cluster_norm), assignments_extract$sample_id)

        if(!sum(de_out2$adjp[, paste0("adjp_", contrast_names)] < FDR_cutoff) == 0){

          pdf.width <- length(samples)*.6+1
          pdf.height <- length(markers)
          pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DEM/Differential.expression_activation.markers_heatmap_lmer.pdf"),
              height = pdf.height, width = pdf.width)
          easyFlow:::plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_cluster_norm,
                                                       sign_adjp = sign_adjp, condition = assignments_extract$condition[mm],
                                                       color_conditions = color_conditions)
          dev.off()

        } # End de_out2 sig not NULL

      } # End try 4

    }
    ########### End : differential expression of activation markers


    {
      message("Based on all events...")
      pdf.width <- 4
      pdf.height <- 2*length(activation.markers)
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DEM/Differential.expr.functionnal.markers_All.clusters.hist.pdf"),
          height = pdf.height, width = pdf.width)

      ## ----diff-expr1-plot-median-expr, fig.height = 8, fig.cap = "Median (arcsinh-transformed) expression of 14 signaling markers calculated from all the cells in a given sample in the PBMC dataset. Values for the two conditions are indicated with different colors: violet for the unstimulated (Ref) and orange for the stimulated with BCR/FcR-XL (BCRXL) samples. Values for each patient are indicated with different shape."----
      ggdf <- reshape2::melt(data.frame(expr_median_sample[functional_markers, ],
                                        antigen = functional_markers), id.vars = "antigen",
                             value.name = "median_expression", variable.name = "sample_id")
      ## Add condition info
      tube.names <- row.names(assignments_extract) %>% gsub(" ",".",.) %>% gsub("-",".",.)
      tube.names <- stringr::str_replace_all(tube.names, "[²&é~{}()'[']|è`ç^à=°+$£¤%µù*><,;:!?./§]", ".")
      tube.names

      mm <- match(ggdf$sample_id, tube.names)
      ggdf$condition <- factor(assignments_extract$condition[mm])
      ggdf$patient_id <- factor(assignments_extract$patient_id[mm])

      plot <- ggplot2::ggplot(ggdf) +
        ggplot2::geom_boxplot(ggplot2::aes(x = condition, y = median_expression, color = condition,
                                           fill = condition),  position = ggplot2::position_dodge(), alpha = 0.5,
                              outlier.color = NA) +
        ggplot2::geom_point(ggplot2::aes(x = condition, y = median_expression, color = condition,
                                         shape = patient_id), alpha = 0.8, position = ggplot2::position_jitterdodge()) +
        ggplot2::facet_wrap(~ antigen, scales = "free", nrow = 5) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
        ggplot2::scale_color_manual(values = color_conditions) +
        ggplot2::scale_fill_manual(values = color_conditions) +
        ggplot2::scale_shape_manual(values = rep(c(0:20),times=1000))
      print(plot)
      dev.off()


      pdf.width <- 2*length(samples)
      pdf.height <- length(functional_markers)
      pdf(file=paste0(dir,"/results/CytofWorkflow/",contrast_names, "/plots/DEM/Differential.expr.functionnal.markers_All.clusters_heatmap.pdf"),
          height = pdf.height, width = pdf.width)

      if(inherits(try(easyFlow:::differential_expression_wrapper(expr_median =
                                                                 expr_median_sample[functional_markers, ],
                                                                 md = assignments_extract, model = "lm", formula = formula_lm, K = K)
                      , silent=TRUE), "try-error")==FALSE) {
        message("with formula_lm...")
        ## ----diff-expr1-fit-model------------------------------------------------
        ## Fit a linear model
        de_out3 <- easyFlow:::differential_expression_wrapper(expr_median =
                                                                expr_median_sample[functional_markers, ],
                                                              md = assignments_extract, model = "lm", formula = formula_lm, K = K)
        apply(de_out3$adjp < FDR_cutoff, 2, table)
        ## ----diff-expr1-fit-model-output-----------------------------------------
        de_output3 <- data.frame(antigen = functional_markers,
                                 expr_median_sample[functional_markers, ], de_out3$pvals, de_out3$adjp)
        #print(head(de_output3), digits=2)
        apply(de_out3$adjp < FDR_cutoff, 2, table)
        write.csv(x=data.frame(de_output3), file = paste0(dir,"/results/CytofWorkflow/",contrast_names,"/tables/DEM/Differential.expr.clustering.markers_All.clusters_Linear.regression.model1.csv"), row.names = TRUE)

        ## ----diff-expr1-plot-heatmap-with-significant-markers
        # "Normalized expression of signaling markers calculated over all the cells in the PBMC dataset
        # that are significantly differentially expressed between BCR/FcR-XL stimulated and unstimulated
        # condition.
        # The heat represents median (arcsinh-transformed) marker expression that was subsequently
        # normalized per marker (rows) to mean of zero and standard deviation of one.
        # The color of the heat varies from blue representing relative under-expression to orange
        # representing relative over-expression. Bar at the top of the heatmap indicates the condition
        # the samples (columns) belong to: violet for the unstimulated (Ref) and orange for
        # the stimulated with BCR/FcR-XL (BCRXL) condition. Numbers in the brackets next to
        # the marker names indicate adjusted p-values and markers are sorted so that markers
        # on the top exhibit the most significant changes between the two conditions.
        # Shown are only the significant markers for which adjusted p-values < 0.05."----

        ## Keep the significant markers and sort them by significance
        sign_markers <- names(which(sort(de_out3$adjp[, paste0("adjp_", contrast_names)]) < FDR_cutoff))
        ## Get the adjusted p-values
        sign_adjp <- de_out3$adjp[sign_markers , paste0("adjp_", contrast_names)]
        ## Normalize expression to mean = 0 and sd = 1
        expr_median_sample_norm <- easyFlow:::normalization_wrapper(expr_median_sample[sign_markers, ])

        if(!sum(de_out3$adjp[, paste0("adjp_", contrast_names)] < FDR_cutoff) == 0){

          mm <- match(colnames(expr_median_sample_norm), assignments_extract$sample_id)
          easyFlow:::plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_norm,
                                                       sign_adjp = sign_adjp,
                                                       condition = assignments_extract$condition[mm],
                                                       color_conditions = color_conditions)
        } # End de_out3 sig not NULL

      } #End Try 5


      if(inherits(try(easyFlow:::differential_expression_wrapper(expr_median =
                                                                 expr_median_sample[functional_markers, ],
                                                                 md = assignments_extract, model = "lmer", formula = formula_lmer, K = K)
                      , silent=TRUE), "try-error")==FALSE) {

        message("with formula_lmer...")
        ## Fit a linear mixed model with patient ID as a random effect
        de_out4 <- easyFlow:::differential_expression_wrapper(expr_median =
                                                                expr_median_sample[functional_markers, ],
                                                              md = assignments_extract, model = "lmer", formula = formula_lmer, K = K)
        apply(de_out4$adjp < FDR_cutoff, 2, table)
        ## ----diff-expr1-fit-model-output-----------------------------------------
        de_output4 <- data.frame(antigen = functional_markers,
                                 expr_median_sample[functional_markers, ], de_out4$pvals, de_out4$adjp)
        #print(head(de_output4), digits=2)
        apply(de_out4$adjp < FDR_cutoff, 2, table)
        write.csv(x=data.frame(de_output4), file = paste0(dir,"/results/CytofWorkflow/",contrast_names,"/tables/DEM/Differential.expr.clustering.markers_All.clusters_Linear.regression.model2.csv"), row.names = TRUE)

        ## ----diff-expr1-plot-heatmap-with-significant-markers
        # "Normalized expression of signaling markers calculated over all the cells in the PBMC dataset
        # that are significantly differentially expressed between BCR/FcR-XL stimulated and unstimulated
        # condition.
        # The heat represents median (arcsinh-transformed) marker expression that was subsequently
        # normalized per marker (rows) to mean of zero and standard deviation of one.
        # The color of the heat varies from blue representing relative under-expression to orange
        # representing relative over-expression. Bar at the top of the heatmap indicates the condition
        # the samples (columns) belong to: violet for the unstimulated (Ref) and orange for
        # the stimulated with BCR/FcR-XL (BCRXL) condition. Numbers in the brackets next to
        # the marker names indicate adjusted p-values and markers are sorted so that markers
        # on the top exhibit the most significant changes between the two conditions.
        # Shown are only the significant markers for which adjusted p-values < 0.05."----

        ## Keep the significant markers and sort them by significance
        sign_markers <- names(which(sort(de_out4$adjp[, paste0("adjp_", contrast_names)]) < FDR_cutoff))
        ## Get the adjusted p-values
        sign_adjp <- de_out4$adjp[sign_markers , paste0("adjp_", contrast_names)]
        ## Normalize expression to mean = 0 and sd = 1
        expr_median_sample_norm <- easyFlow:::normalization_wrapper(expr_median_sample[sign_markers, ])

        mm <- match(colnames(expr_median_sample_norm), assignments_extract$sample_id)

        if(!sum(de_out4$adjp[, paste0("adjp_", contrast_names)] < FDR_cutoff) == 0){

          easyFlow:::plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_norm,
                                                       sign_adjp = sign_adjp,
                                                       condition = assignments_extract$condition[mm],
                                                       color_conditions = color_conditions)
        } # End de_out4 sig not NULL

      } # End try6
      dev.off()
    }

    #'###########################################################################
    #' End section Diff expr All markers
    #############################################################################

  } # End CytofWorkFlow

}

