#' Read and combine multiple tps files
#'
#' Read multiple tps files to obtain landmark coordinates and combine them into a single array
#'
#' This is a wrapper of \code{\link{readland.tps}} to allow reading landmark coordinates, in 2D or 3D, from several tps files , and compiling them into an array for proceeding with GM procedures.
#' 
#' The arguments specID and negNA of \code{\link{readland.tps}} can be directly set through this function.
#' Note that if specID is set to either "None" or "ID", a check for duplicate specimen names is not possible and specimens will be numbered with 1:N, where N the total number of specimens.  
#' 
#'
#' @param filelist A vector containing the file paths to all the tps files to be compiled
#' @param ... other arguments to be passed to \code{\link{readland.tps}}
#' @author Antigoni Kaliontzopoulou
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the total number of specimens across all tps files included in the folder read. 
#' @export

readmulti.tps <- function(filelist, ... ){
  tps.list <- filelist
  readland.args <- list(...)
  if(is.null(readland.args$specID)) readland.args$specID <- "None"
  
  file.ext <- substr(tps.list, nchar(tps.list)-3, nchar(tps.list))
  if(!all(file.ext%in%c(".tps", ".TPS"))) 
    stop("File list includes files in a format other than tps, please ammend")
  
  dt.dims <- sapply(1:length(tps.list), function(x){
    dim(readland.tps(tps.list[x], warnmsg = F, ...))
  }, simplify = T)
  p1 <- dt.dims[1, 1]; k1 <- dt.dims[2, 1]; n1 <- dt.dims[3, 1]
  
  if(any(dt.dims[1,]!=p1)) stop("Input tps files include different numbers of landmarks, please correct")
  
  if(any(dt.dims[2,]!=k1)) stop("Input tps files include landmarks in different dimensions (2D and 3D), please correct")
  
  all.lms <- NULL
  for(f in 1:length(tps.list)){       
    lms <- two.d.array(readland.tps(tps.list[f], ...))
    all.lms <- rbind(all.lms, lms)
  }
  all.lms <- arrayspecs(all.lms, p1, k1)
  
  if(any(table(dimnames(all.lms)[3])!=1)) {
    if(readland.args$specID != "imageID") {
      dimnames(all.lms)[[3]] <- 1:dim(all.lms)[3]
    } else {
      warning("Input files seem to include repeated specimen names")
    }
  }
  
  return(all.lms)
}
