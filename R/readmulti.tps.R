#' Read and combine multiple tps files
#'
#' Read multiple tps files to obtain landmark coordinates and combine them into a single array
#'
#' This is a wrapper of \code{\link{readland.tps}} to allow reading landmark coordinates, in 2D or 3D, from several tps files , and compiling them into an array for proceeding with GM procedures.
#' 
#' The arguments specID and negNA of \code{\link{readland.tps}} can be directly set through this function.
#' 
#'
#' @param filelist A vector containing the file paths to all the tps files to be compiled
#' @param ... other arguments to be passed to \code{\link{readland.tps}}
#' @author Antigoni Kaliontzopoulou
#' @author Michael Collyer
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the total number of specimens across all tps files included in the folder read. 
#' @export

readmulti.tps <- function(filelist, ... ){
  tps.list <- filelist
  readland.args <- list(...)
  if(is.null(readland.args$specID)) readland.args$specID <- "None"
  if(is.null(readland.args$negNA)) readland.args$negNA <- FALSE
  if(is.null(readland.args$readcurves)) 
    readland.args$readcurves <- FALSE
  if(is.null(readland.args$warnmsg)) readland.args$warnmsg <- FALSE
  
  file.ext <- substr(tps.list, nchar(tps.list)-3, nchar(tps.list))
  if(!all(file.ext%in%c(".tps", ".TPS"))) 
    stop("\nFile list includes files in a format other than tps, please ammend",
         call. = FALSE)
  
  readland.args$file <- tps.list[[1]]
  all.lms <- lapply(1:length(tps.list), function(j){
    readland.args$file <- tps.list[[j]]
    do.call(.readland.tps, readland.args)
  })
  
  specnames <- unlist(lapply(all.lms, names))
  if(readland.args$specID == "None") specnames <- 1:length(specnames)
  if(length(unique(specnames)) < length(specnames)) 
    warning("\n\nInput files seem to include repeated specimen names.\n", 
            call. = FALSE)
  
  lm.tab<- do.call(cbind, lapply(all.lms, function(x) sapply(x, dim)))
  lm.check <- apply(lm.tab, 1, unique)
  bad.lm <- is.list(lm.check)
  
  if(bad.lm) {
    lmo <- t(lm.tab)
    colnames(lmo) <- c("p", "k")
  } else {
    lmo <- unlist(all.lms, recursive = FALSE)
    names(lmo) <- specnames
    lmo <- simplify2array(lmo)
  }

  if(bad.lm) {
    cat("\nInput tps files include either different numbers of landmarks\n")
    cat("or different dimensions (2D or 3D) for some landmarks.\n\n")
    cat("A table of dimensions is returned rather than an array, so the issue can be investigated.\n\n")
    print(lmo)
  }
  
  return(lmo)
}
