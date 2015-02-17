#' Convert (p x k x n) data array into 2D data matrix
#'
#' Convert a three-dimensional array of landmark coordinates into a 2-dimensional matrix 
#'
#' This function converts a (p x k x n) array of landmark coordinates into a 2-dimensional 
#'  matrix. The latter format of the shape data is useful for performing subsequent statistical 
#'  analyses in R (e.g., PCA, MANOVA, PLS, etc.). Row labels are preserved if included in 
#'  the original array. 
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @keywords utilities
#' @export
#' @author Dean Adams and Emma Sherratt
#' @return Function returns a two-dimensional matrix of dimension (n x [p x k]), where rows 
#'   represent specimens and columns represent variables. 
#' @examples
#' data(plethodon) 
#' plethodon$land    #original data in the form of 3D array
#' 
#' two.d.array(plethodon$land)   # Convert to a 2D data matrix
two.d.array<-function(A){  
  pxk <- dim(A)[1]*dim(A)[2]
  n <- dim(A)[3]
  tmp <- aperm(A, c(3,2,1))
  dim(tmp) <- c(n,pxk)
  rownames(tmp)<-dimnames(A)[[3]]  
  return(tmp)
}