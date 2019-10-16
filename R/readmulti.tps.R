#' Read and combine multiple tps files
#'
#' Read multiple *.tps files from a directory to obtain landmark coordinates and combine them into a single array
#'
#' This is a wrapper of \code{\link{readland.tps}} to allow reading landmark coordinates, in 2D or 3D, from several tps files in a specific folder, and compiling them into an array for proceeding with GM procedures.
#' 
#' The arguments specID and negNA of \code{\link{readland.tps}} can be directly set through this function. 
#' 
#'
#' @param folder A folder containing all the *.tps files to be compiled
#' @param ... other arguments to be passed to \code{\link{readland.tps}}
#' @author Antigoni Kaliontzopoulou
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the total number of specimens across all tps files included 
#'   in the folder read. 

readmulti.tps <- function(folder, ... ){
  tps.list <- list.files(folder, pattern = ".tps")
  
  if(length(tps.list)==0) stop("No tps files found in the indicated folder")
  
  dt.dims <- sapply(1:length(tps.list), function(x){
    dim(readland.tps(paste(folder, tps.list[x], sep="")))
  }, simplify = T)
  p1 <- dt.dims[1, 1]; k1 <- dt.dims[2, 1]; n1 <- dt.dims[3, 1]
  
  if(any(dt.dims[1,]!=p1)) stop("Input tps files include different numbers of landmarks, please correct")
  
  if(any(dt.dims[2,]!=k1)) stop("Input tps files include landmarks in different dimensions (2D and 3D), please correct")
  
  all.lms <- NULL
  for(f in 1:length(tps.list)){       
    lms <- two.d.array(readland.tps(paste(folder, tps.list[f], sep=""), specID = ...))
    all.lms <- rbind(all.lms, lms)
  }
  all.lms <- arrayspecs(all.lms, p1, k1)
  return(all.lms)
}
