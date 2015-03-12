#' Identify specimen closest to the mean of a set of aligned specimens
#' 
#' A function to identify which specimen lies closest to the estimated mean 
#' shape for a set of aligned specimens.
#' 
#' Function takes an array of aligned specimens (such as made by \code{\link{gpagen}}, 
#' calculates the distance of each to the estimated mean shape, and returns the name and 
#' address of the closest specimen. This function can be used
#' to identify the specimen to be used by \code{\link{warpRefMesh}}. 
#' 
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of aligned specimens
#' @export
#' @seealso  \code{\link{warpRefMesh}}
#' @keywords utilities
#' @author Emma Sherratt
#' @return Function returns the name and address of the specimen closest to the mean of the set of 
#' aligned specimens.
findMeanSpec <- function(A){
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  ref <- mshape(A)
  x <- two.d.array(A)
  x <- rbind(x, as.vector(t(ref)))
  dists <- as.matrix(dist(x))[,(nrow(x))]
  spec <- which(dists == min(dists[dists>0]))
  return(spec)
}
