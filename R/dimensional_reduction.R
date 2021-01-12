#' Dimensionality reduction
#'
#' @importFrom Rtsne Rtsne
#' @importFrom FNN KL.divergence
#' @importFrom FNN get.knnx
#' @importFrom caret createDataPartition
#'
#' @param dir a character string specifying the project directory
#'
#' @export

dimensional_reduction <- function(dir = NULL){

  {
    if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

    message("\n")
    message(" ########## Running dimensionality reduction ########## ")

    message("Loading data...")
    source(paste0(dir,"/attachments/easyFlow.R"))
    load(paste0(dir,"/Rdata/easyFlow.RData"))
    load(paste0(dir,"/Rdata/som_map.RData"))
    source(FItSNE.dir, chdir=T)

    subsample_data_table <- NULL
    set.seed(seed)

    data_table <- data_table[, clustering.markers]
  }

  # Initialisation matrix on clusters
  if(init.matrix == TRUE){

    message("Computing initialisation matrix...")

    {
      # Extract Median expr for each clusters
      median.expr <- data.frame()
      median <- data.frame()

      median.expr <- easyFlow:::helper_cluster_medians(data_table, som.cluster)
      colnames(median.expr) <- c(colnames(data_table)) #"seed", "som.cluster",

    } # End Extract Median expr for each clusters

    {
      tsne_out1 <- Rtsne::Rtsne(as.matrix(median.expr),
                                dims = 2, #n.dims,
                                perplexity = 6,#1.5*log(nrow(data_table))/ncol(data_table),#perplexity.value,
                                theta = .3, #theta.value,
                                pca = FALSE,
                                pca_center = FALSE,
                                pca_scale = FALSE,
                                max_iter = 5000, #n.iter,
                                mom_switch_iter = 500, #250
                                stop_lying_iter = 500,
                                momentum = 0.5,
                                final_momentum = 0.8,
                                eta = 200,
                                exaggeration_factor = 12
      )

    } # end tnse on clusters

  } # End initialisation matrix


  cluster_id <- sort(unique(som.cluster))


  ###### tsne_prediction #############
  if(!tsne_prediction == "none"){

    # Avoid subsample dataset to be larger than all events
    if(tsne_pred_n_events > length(som.cluster)){
      tsne_pred_n_events <- length(som.cluster)
    }

    # Cluster stratified method
    if(tsne_pred_sampling == "CS"){  # cluster_stratified or DCP deep_cluster_profile
      # Stratified sampling
      sample_pct <- tsne_pred_n_events / length(som.cluster)
      ID <- caret::createDataPartition(y = som.cluster, p= sample_pct, list = FALSE)

      # subsample data_table
      subsample_data_table <- data_table[ID,]
      subsample_som.cluster <- som.cluster[ID]

      mm <- match(som.cluster, cluster_id)
      mm

      if(init.matrix == TRUE){
        init <- cbind(as.numeric(tsne_out1$Y[,1][mm]),
                      as.numeric(tsne_out1$Y[,2][mm]))
        colnames(init) <- c("tsne1", "tsne2")
        init <- init[ID,]
      } else {
        init <- NULL
      }

    } # End CS

    # Deep cluster profile method
    if(tsne_pred_sampling == "DCP"){

      {
        # Select event_per_clus IDs in each clusters (100 clusters x 500 events = 50.000events)
        message("Extracting reference clusters (deep clusters profile)...")
        ID <- NULL
        i <- sort(unique(som.cluster))[1]
        event_per_clus <- floor(tsne_pred_n_events / length(as.numeric(levels(som.cluster))))
        message(paste0("clusters : ", length(as.numeric(levels(som.cluster))) ))
        message(paste0("max events per cluster : ", event_per_clus ))
        ID_cluster <- sort(unique(som.cluster))

        for(i in ID_cluster){
          if(sum(as.numeric(som.cluster) == i) > 0){
            clustering.map.id <- which(as.numeric(som.cluster) == i)
          } else {
            clustering.map.id <- NULL
          }

          if(is.null(ID)){

            if(length(clustering.map.id) > event_per_clus){
              ID <- sample(clustering.map.id, event_per_clus, replace = FALSE)
            }
            if(length(clustering.map.id) <= event_per_clus && length(clustering.map.id) > 0){
              ID <- clustering.map.id
            }

          } else {

            if(length(clustering.map.id) > event_per_clus){
              ID <- c(ID, sample(clustering.map.id, event_per_clus, replace = TRUE))
            }
            if(length(clustering.map.id) <= event_per_clus && length(clustering.map.id) > 0){
              ID <- c(ID, clustering.map.id)
            }

          }
        } # End for cluster

        # Cluster Map
        #subsample_data_table <- data_table[sort(as.numeric(na.omit(ID))),]
        #subsample_som.cluster <- som.cluster[sort(as.numeric(na.omit(ID)))]
        subsample_data_table <- data_table[ID,]
        subsample_som.cluster <- som.cluster[ID]


        if(init.matrix == TRUE){
          ## New clustering
          mm <- match(som.cluster, cluster_id)
          mm

          init <- cbind(as.numeric(tsne_out1$Y[,1][mm]),
                        as.numeric(tsne_out1$Y[,2][mm]))
          colnames(init) <- c("tsne1", "tsne2")
          init <- init[ID]
        } else {
          init <- NULL
        }

        # Check nrow for each ref cluster
        #for(lev in levels(factor(som.cluster))){print(paste0("som.cluster ", lev, " nrow ", sum(subsample_som.cluster==lev)))}

      } # End reference clusters

    } #end DCP

    {
      message("Computing tSNE map...")

      tsne_out2 <- fftRtsne(as.matrix(subsample_data_table),
                            dims = 2, #n.dims,
                            perplexity = 0,#perplexity.value,
                            perplexity_list = perplexity.value,
                            theta = theta.value,
                            max_iter = n.iter,
                            fft_not_bh = tsne_algo, #FIt-SNE = TRUE
                            # if theta is nonzero, this determins whether to
                            # use FIt-SNE or Barnes Hut approximation. Default is FIt-SNE.
                            # set to be True for FIt-SNE
                            ann_not_vptree = TRUE, #approximate nearest neighbors,
                            # use vp-trees (as in bhtsne) or approximate nearest neighbors (default).
                            # set to be True for approximate nearest neighbors
                            initialization =  init,
                            momentum = 0.5,
                            final_momentum = 0.8,
                            exaggeration_factor = 12,
                            no_momentum_during_exag = 1,
                            # Set to 0 to use momentum and other optimization tricks. 1 to do plain,vanilla
                            # gradient descent (useful for testing large exaggeration coefficients)
                            stop_early_exag_iter = 500, #n.iter/4,
                            # When to switch off early exaggeration. Default 250.
                            start_late_exag_iter = -1,
                            # When to start late exaggeration. set to -1 to not use late exaggeration. Default -1.
                            # late_exag_coeff - Late exaggeration coefficient.
                            # Set to -1 to not use late exaggeration. Default -1
                            df = freedom.value
                            # Degree of freedom of t-distribution, must be greater than 0.
                            # Values smaller than 1 correspond to heavier tails, which can often
                            # resolve substructure in the embedding. See Kobak et al. (2019) for
                            # details. Default is 1.0
      );
      colnames(tsne_out2) <- c("tsne1", "tsne2")

    } # end FItSNE


    {
      # tSNE prediction

      ##############################
      # Train the model with 5 cvf #
      ##############################

      suppressWarnings(rm(list = c("intrain", "X_train", "X_test", "y_train", "y_test", "nn", "dist_min", "dist_max", "dist_k", "coeff", "x_init", "y_init", "pred_tsne", "rmse_val")))

      n_cvf <- 3
      k <- 2
      rmse_val <- data.frame(matrix(data = NA, nrow = 10, ncol = 4))
      i <- 1

      if(tsne_prediction == "k1init" || tsne_prediction == "k1final"){
        k_min <- 1
        k_max <- 1
      }
      if(tsne_prediction == "complete") {
        k_min <- 2
        k_max <- 5
      }

      for(k in k_min:k_max){

        message(paste0("k = ", k))

        for(n_cvf in c(1:5)){

          # Stratified sampling
          intrain <- caret::createDataPartition(y = som.cluster[ID], p= 0.9, list = FALSE)
          X_train <- subsample_data_table[intrain, ]
          X_test <- subsample_data_table[-intrain, ]

          y_train <- tsne_out2[intrain,]
          y_test <- tsne_out2[-intrain,]

          # Get nearest neighbors
          knn2 <- FNN::get.knnx(X_train,
                                X_test,
                                k=k, algorithm=c("kd_tree") ) # "kd_tree", "cover_tree", "CR", "brute"))

          # Coefficients weigthing
          coeff <- data.frame(matrix(data = NA, nrow = nrow(X_test), ncol = k))
          y_init <- x_init <- data.frame(matrix(data = NA, nrow = nrow(X_test), ncol = k))

          dist_min <- apply(knn2$nn.dist, 1, function(x) min(x))
          dist_max <- apply(knn2$nn.dist, 1, function(x) max(x))

          for(i in c(1:k)){

            # Compute weighted tsne coord
            coeff[ , i] <- abs( 1 - ( ( ( knn2$nn.dist[ , i] + 1) - ( dist_min + 1) ) / ( dist_max + 1) ) )

            # Predicted tsne coordinates
            x_init[ , i] <- coeff[ , i] * y_train[knn2$nn.index[ , i] , 1]
            y_init[ , i] <- coeff[ , i] * y_train[knn2$nn.index[ , i] , 2]

          } # end for k

          pred_tsne <- matrix(c(rowMeans(x_init),rowMeans(y_init)), ncol=2)

          rmse_val[k, n_cvf] <- caret::RMSE(pred_tsne, y_test)

        } # end for n_cvf

      } # end for k in k_min : k_max

      rmse_val[ , n_cvf+1] <- rowMeans(rmse_val[ , c(1:3)])
      rmse_val


      ##################################
      # Apply the model to the dataset #
      ##################################
      if(tsne_prediction == "k1init" || tsne_prediction == "k1final"){

        # Best k value
        best_k <- 1
        message(paste0("k value : ", best_k))
        message(paste0("Mean RMSE : ", sort(rmse_val[ , n_cvf+1])[1]))

        # Get nearest neighbors
        knn <- FNN::get.knnx(subsample_data_table,
                             data_table,
                             k=best_k, algorithm=c("kd_tree") ) # "kd_tree", "cover_tree", "CR", "brute"))

        # Predicted tsne coordinates
        x_init <- tsne_out2[knn$nn.index[ , 1] , 1]
        y_init <- tsne_out2[knn$nn.index[ , 1] , 2]

        tsne.coord <- cbind(x_init,y_init)
        init <- tsne.coord

      }

      if(tsne_prediction == "complete") {

        # Best k value
        best_k <- order(rmse_val[ , n_cvf+1])[1]
        message(paste0("Best k value : ", best_k))
        message(paste0("Mean RMSE : ", sort(rmse_val[ , n_cvf+1])[1]))

        best_k <- 2
        message(paste0("But k value kept : ", best_k))
        message(paste0("Mean RMSE : ", sort(rmse_val[ , n_cvf+1])[1]))

        # Get nearest neighbors
        knn <- FNN::get.knnx(subsample_data_table,
                             data_table,
                             k=best_k, algorithm=c("kd_tree") ) # "kd_tree", "cover_tree", "CR", "brute"))

        # Coefficients weigthing
        coeff <- data.frame(matrix(data = NA, nrow = nrow(data_table), ncol = best_k))
        y_init <- x_init <- data.frame(matrix(data = NA, nrow = nrow(data_table), ncol = best_k))

        dist_min <- apply(knn$nn.dist, 1, function(x) min(x))
        dist_max <- apply(knn$nn.dist, 1, function(x) max(x))

        # Compute weighted tsne coord
        for(i in c(1:best_k)){
          coeff[ , i] <- abs( 1 - ( ( ( knn$nn.dist[ , i] + 1) - ( dist_min + 1) ) / ( dist_max + 1) ) )

          # Predicted tsne coordinates
          x_init[ , i] <- coeff[ , i] * tsne_out2[knn$nn.index[ , i] , 1]
          y_init[ , i] <- coeff[ , i] * tsne_out2[knn$nn.index[ , i] , 2]

        }

        tsne.coord <- cbind(rowMeans(x_init),rowMeans(y_init))

      } # end tsne_prediction

    } # end prediction on knn

  } # end tsne prediction


  # Compute FItSNE on all data
  message("Compute tSNE on all the dataset...")

  if(tsne_prediction == "none" || tsne_prediction == "k1init"){

    if(init.matrix == TRUE && tsne_prediction == "none"){
      exaggeration_factor <- 3

      ## New clustering
      mm <- match(som.cluster, cluster_id)
      mm

      init <- cbind(as.numeric(tsne_out1$Y[,1][mm]),
                    as.numeric(tsne_out1$Y[,2][mm]))
      colnames(init) <- c("tsne1", "tsne2")

      perplexity.value2 <- perplexity.value
      stop_early_exag_iter <- ceiling(n.iter/4)
      n.iter2 <- n.iter
      learning_rate <- 200

    }

    if(init.matrix == FALSE && tsne_prediction == "none"){
      exaggeration_factor <- 12
      init <- NULL
      perplexity.value2 <- perplexity.value
      stop_early_exag_iter <- ceiling(n.iter/4)
      n.iter2 <- n.iter
      learning_rate <- 200
    }

    if(tsne_prediction == "k1init"){
      exaggeration_factor <- 3
      perplexity.value2 <- 50
      stop_early_exag_iter <- 1000
      n.iter2 <- 5000
      learning_rate <- 1000
    }

    tsne.coord <- fftRtsne(as.matrix(data_table),
                           dims = 2, #n.dims,
                           perplexity = 0,#perplexity.value,
                           perplexity_list = perplexity.value2,
                           theta = theta.value,
                           max_iter = n.iter2,
                           fft_not_bh = tsne_algo, #FIt-SNE = TRUE
                           # if theta is nonzero, this determins whether to
                           # use FIt-SNE or Barnes Hut approximation. Default is FIt-SNE.
                           # set to be True for FIt-SNE
                           ann_not_vptree = TRUE, #approximate nearest neighbors,
                           # use vp-trees (as in bhtsne) or approximate nearest neighbors (default).
                           # set to be True for approximate nearest neighbors
                           initialization =  init,
                           momentum = 0.5,
                           final_momentum = 0.8,
                           learning_rate = learning_rate,
                           exaggeration_factor = exaggeration_factor,
                           no_momentum_during_exag = 1,
                           # Set to 0 to use momentum and other optimization tricks. 1 to do plain,vanilla
                           # gradient descent (useful for testing large exaggeration coefficients)
                           stop_early_exag_iter = stop_early_exag_iter,
                           # When to switch off early exaggeration. Default 250.
                           start_late_exag_iter = -1,
                           # When to start late exaggeration. set to -1 to not use late exaggeration. Default -1.
                           # late_exag_coeff - Late exaggeration coefficient.
                           # Set to -1 to not use late exaggeration. Default -1
                           df = freedom.value
                           # Degree of freedom of t-distribution, must be greater than 0.
                           # Values smaller than 1 correspond to heavier tails, which can often
                           # resolve substructure in the embedding. See Kobak et al. (2019) for
                           # details. Default is 1.0
    );

  } # end FItSNE on all data

  {
    #message("Computing final Kullback–Leibler divergence...")
    #pca1 <- stats::prcomp(data_table, rank. = 2)
    #pca2 <- stats::prcomp(tsne.coord, rank. = 2)
    #KLD <- FNN::KL.divergence(X = pca1$x,
    #                          Y = pca2$x,
    #                          k = 10,
    #                          algorithm=c("kd_tree"))

    #message(paste0("Final KL Divergence (k = 10, distances are PCA-based) = ", KLD[10]))


    message(paste0("Saving tSNE results..."))
    clustering.markers.tsne <- clustering.markers
    suppressWarnings(dir.create(paste0(dir,"/Rdata/")))
    save(list = c("clustering.markers.tsne",
                  "tsne.coord", "seed", #"KLD",
                  "perplexity.value", "theta.value", "n.iter", "freedom.value",
                  "init.matrix", "tsne_prediction", "tsne_pred_sampling", "tsne_pred_n_events", "tsne_algo"
    ),
    file = paste0(dir,"/Rdata/tsne.RData"),
    compress = "gzip"
    )

  }

} # end dim red
