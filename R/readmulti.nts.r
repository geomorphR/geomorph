#' Read and combine multiple nts files
#'
#' Read multiple nts (or .dta) files to obtain landmark coordinates and combine them into a single array
#'
#' This is a wrapper of \code{\link{readland.nts}} to allow reading landmark coordinates, in 2D or 3D, 
#' from several nts (or .dta) files, and compiling them into an array for proceeding with GM procedures.
#' 
#' See \code{\link{readland.nts}} for adequately formatting NTS files and requirements that need to be met 
#' for that (and therefore this) function to work correctly.
#' 
#'
#' @param filelist A vector containing the file paths to all the nts files to be compiled
#' @author Antigoni Kaliontzopoulou
#' @return Function returns a 3D array (p x k x n), where p is the number of landmark points, k is 
#'   the number of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension 
#'   of this array contains names for each specimen, which are retrieved from the nts specimen labels, 
#'   if those are available. Specimens in nts files without labels are named as filename_# where # are consecutive
#'   numbers.
#' @references  Rohlf, F. J. 2012 NTSYSpc: Numerical taxonomy and multivariate analysis system. Version 
#'   2.2. Exeter Software, New York.
#' @export


readmulti.nts <- function(filelist){
  nts.list <- filelist
  
  file.ext <- substr(nts.list, nchar(nts.list)-3, nchar(nts.list))
  if(!all(file.ext%in%c(".nts", ".NTS", ".dta", ".DTA"))) 
    stop("File list includes files in a format other than nts or dta, please ammend")
  
  dt.dims <- sapply(1:length(nts.list), function(x){
    dim(readland.nts(nts.list[x]))
  }, simplify = T)
  p1 <- dt.dims[1, 1]; k1 <- dt.dims[2, 1]; n1 <- dt.dims[3, 1]
  
  if(any(dt.dims[1,]!=p1)) stop("Input tps files include different numbers of landmarks, please correct")
  if(any(dt.dims[2,]!=k1)) stop("Input tps files include landmarks in different dimensions (2D and 3D), please correct")
  
  all.lms <- NULL
  for(f in 1:length(nts.list)){       
    lms <- two.d.array(readland.nts(nts.list[f]))
    if(is.null(rownames(lms))) {
      rownames(lms) <- paste(nts.list[f], 1:nrow(lms), sep = "_")
    }
    all.lms <- rbind(all.lms, lms)
  }
  all.lms <- arrayspecs(all.lms, p1, k1)
  
  return(all.lms)
}
