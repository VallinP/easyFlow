##########################################################################################
#' Metaclustering using Hierarchical consensus clustering
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#' SSE and findElbow function are adapted from FlowSOM R package
#' Sofie Van Gassen et al., Cytometry Part A, 2015
#' https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625
#' http://bioconductor.org/packages/release/bioc/html/FlowSOM.html
#########################################################################################
#'
#' @source https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625
#' @source http://bioconductor.org/packages/release/bioc/html/FlowSOM.html
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom NbClust NbClust
#'
#' @param dir a character string specifing the project directory
#'
#' @export

metaClustering_HC <- function(dir = NULL){

  if(is.null(dir)){dir <- easyFlow:::choose_directory("Select the directory containing the project")}

  message("\n")
  message(" ########## Running metaclustering ########## ")

  message("Loading data...")
  source(paste0(dir,"/attachments/easyFlow.R"))
  load(paste0(dir,"/Rdata/easyFlow.RData"))
  load(paste0(dir,"/Rdata/som_map.RData"))

  message(paste0("Hierarchical consensus clustering. maxK = ", n.metaclusters, " metaclusters"))

  som.cluster <- factor(x=som.res$unit.classif, levels=as.numeric(levels(factor(som.res$unit.classif))))

  # Compute CCP for n.metaclusters (1 to n.metametaclusters)
  SOM.metaclustering <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(t(som.res$codes[[1]]),
                                                                                    maxK = n.metaclusters,
                                                                                    reps = 100,
                                                                                    pItem = 0.9, pFeature = 1,
                                                                                    title = tempdir(), plot = "pdf",
                                                                                    verbose = FALSE,
                                                                                    innerLinkage = HCC_method, #"average", #ward.D, simple, average, complete
                                                                                    finalLinkage = HCC_method, #"average", #ward.D
                                                                                    clusterAlg = "hc",
                                                                                    distance =  "euclidean", #"pearson",
                                                                                    seed = seed))
  #SOM.metaclustering[[n.metaclusters]]["consensusClass"][[1]]


  # For each consensus, clusters are associated to their respective metacluster
  # MetaClustering_consensus ? (cf. vignette)
  message(paste0("Formating metaclusters results..."))
  metacl <- data.frame(matrix(nrow = nrow(som.res$codes[[1]]), ncol = (n.metaclusters)))

  min.metaclusters <- 2
  for(i in n.metaclusters:min.metaclusters){
    metacl[,i] <- SOM.metaclustering[[i]]$consensusClass
  } #end for n metaclusters


  # Test unimodality of metaclusters
  unimodal.res <- NULL
  if(FALSE){
    message(paste0("Testing unimodality of metaclusters..."))

    unimodal.res <- as.list(rep(TRUE, n.metaclusters))

    for(i in rev(2:n.metaclusters)){

    classif <- metacl[,i][som.res$unit.classif]
    message(paste0("With ", i, " metaclusters selected"))

    unimodal.res[[i]] <- easyFlow:::SPADEVizR_unimodality(df = som.res$codes[[1]],
                                                          cluster = classif,
                                                          clusters = NULL,
                                                          uniform.test = "unimodality", #"unimodality", "spread", "both"
                                                          th.pvalue = 0.05,
                                                          th.IQR = 2,
                                                          #density.PDFfile = NULL #"qcUniformClusters_density.pdf",
                                                          #density.PDFfile.dim = c(17, 10),
                                                          #heatmap.PDFfile = NULL #"qcUniformClusters_heatmap.pdf"
                                                          tile.color = "black",
                                                          verbose = TRUE)
    }

    unimodal.res2 <- NULL
    for(i in rev(2:n.metaclusters)){
      if( (i != which(is.null(unimodal.res2[i])) ) || (i != which(is.na(unimodal.res2[i])) ) ){
      unimodal.res2[i] <- sum(unimodal.res[[i]]$accuracy.matrix$uniform.cluster == TRUE ) / n.metaclusters
    } else {
      unimodal.res2[i] <- FALSE
    }

    } # end n.metaclusters

    max_unimodal <- max(unimodal.res2[!is.na(unimodal.res2)])
    max_unimodal_ix <- which(unimodal.res2 == max_unimodal)
    unimodal_ix <- min(max_unimodal_ix)

    n.metaclusters <- unimodal_ix

    message("For ", n.metaclusters, " metaclusters, ", 100*unimodal.res2[n.metaclusters], "% of metaclusters are unimodal.")

  } # end HC_unimodal

  elbow_point_res <- NULL
  if(use_elbow_point == TRUE){
    # From FlowSOM R package
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


    if(n.metaclusters > 50){max <- 50} else {max <- n.metaclusters}
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
    n.metaclusters <- findElbow(elbow_point_res)

    plot(1:max, elbow_point_res, type = "b", xlab = "Number of Clusters",
         ylab = "Within groups sum of squares")
    abline(v = (n.metaclusters + .01), col = "blue")

  }

  if(HC_unimodal == TRUE){
    message(paste0("Testing unimodality of metaclusters..."))

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

    elbow_point <- FlowSOM:::findElbow(res)
    elbow_point2 <- which(rev(res) == res[elbow_point]) #(n.metaclusters-1) - elbow_point
    plot(unimodal.perc)
    abline(v = (elbow_point2 - .01), col = "red")
    message("Based on elbow point of unimodal percentages")
    message("Found ", elbow_point2, " metaclusters")
    message("Unimodal percentage : ", round(res[elbow_point], 2), "%")

    nbclus <- NbClust::NbClust(res, method = "ward.D", min.nc = 2, max.nc = 2, index = c("silhouette"))
    silh.score <- min(which(nbclus$Best.partition == 2))
    silh.score2 <- which(rev(res) == res[max(which(nbclus$Best.partition == 1))]) #(n.metaclusters-1) - silh.score
    message("Based on sihouette score of unimodal percentages")
    message("Found ", silh.score2, " metaclusters")
    message("Unimodal percentage : ", round(res[silh.score], 2), "%")
    abline(v = (silh.score2 + .01), col = "blue")

    n.metaclusters <- elbow_point2

  } # end HC_unimodal

  metacluster <- metacl[,n.metaclusters][som.res$unit.classif]

  # Save Hierarchical Consensus Clustering
  message("Saving metaclustering results...")
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

}
