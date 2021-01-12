#'#########################################################################################
#' Function to choose FItSNE install directory.
#'
#' Patrice Vallin, easyFlow, Sept 2019.
#'########################################################################################
#'
#' @param dir a character string specifying the project directory
#'
#' @return FItSNE.dir a character string containing the path to FItSNE folder
#'

choose_FItSNE.dir <- function(dir){

  FItSNE.dir <- NULL

  if(!file.exists(paste0(dir, "/attachments/FItSNE.path.txt"))){
    caption.value <- "Select FItSNE directory"
    message(caption.value)
    FItSNE.dir <- easyFlow:::choose_directory(caption.value)
    suppressWarnings(dir.create(paste0(dir, "/attachments/")))
    write.table(x=paste0(FItSNE.dir,"/fast_tsne.R"),file=paste0(dir, "/attachments/FItSNE.path.txt"))
  }

  if(suppressWarnings(file.exists(paste0(dir, "/attachments/FItSNE.path.txt"))) &&
     suppressMessages(class(try(source(as.character(read.table(paste0(dir, "/attachments/FItSNE.path.txt"))[,1])), silent = T)) == "try-error")){
    FItSNE.dir <- "error"
    message("FItSNE directory detected.")
    message("An error occured during FItSNE initialisation. Please choose a new directory.")
  }

  if(suppressWarnings(file.exists(paste0(dir, "/attachments/FItSNE.path.txt"))) &&
     !suppressMessages(class(try(source(as.character(read.table(paste0(dir, "/attachments/FItSNE.path.txt"))[,1])), silent = T)) == "try-error")){
    message("FItSNE directory detected.")
    FItSNE.dir <- as.character(read.table(paste0(dir, "/attachments/FItSNE.path.txt"))[,1])
  }

  if(FItSNE.dir == "error"){

    caption.value <- "Select FItSNE directory"
    message(caption.value)
    FItSNE.dir <- easyFlow:::choose_directory(caption.value)

    if(suppressMessages(class(try(source(paste0(FItSNE.dir,"/fast_tsne.R")), silent = T)) == "try-error")){
      stop("An error occured during FItSNE initialisation. Try to reinstall it.")
    }
    suppressWarnings(dir.create(paste0(dir, "/attachments/")))
    write.table(x=paste0(FItSNE.dir,"/fast_tsne.R"),file=paste0(dir, "/attachments/FItSNE.path.txt"))
  }
  return(FItSNE.dir)
}
