#' Convert (p x k x n) data array into 2D data matrix
#'
#' Convert a three-dimensional array of landmark coordinates into a two-dimensional matrix 
#'
#' This function converts a (p x k x n) array of landmark coordinates into a two-dimensional 
#'  matrix (n x [p x k]). The latter format of the shape data is useful for performing subsequent statistical 
#'  analyses in R (e.g., PCA, MANOVA, PLS, etc.). Row labels are preserved if included in 
#'  the original array. 
#'
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param sep An optional argument for variable labeling, combining landmark labels (e.g., 1, 2, 3, ...)
#' and partial dimension labels (e.g., "x", "y", and "z"), much like the \code{\link{paste}} function.  
#' The default is sep = ".", but this can be changed to any separator.  
#' One should make sure to match separators with \code{\link{arrayspecs}} if switching between matrices and arrays.
#' 
#' @keywords utilities
#' @export
#' @author Dean Adams and Emma Sherratt
#' @return Function returns a two-dimensional matrix of dimension (n x [p x k]), where rows 
#'   represent specimens and columns represent variables.
#' @seealso \code{\link{arrayspecs}} 
#' @examples
#' \dontrun{
#' data(plethodon) 
#' plethodon$land    #original data in the form of 3D array
#' 
#' two.d.array(plethodon$land)   # Convert to a 2D data matrix
#' }
two.d.array<-function(A, sep = "."){  
  pxk <- dim(A)[1]*dim(A)[2]
  n <- dim(A)[3]
  tmp <- aperm(A, c(3,2,1))
  dim(tmp) <- c(n,pxk)
  rownames(tmp)<-dimnames(A)[[3]] 
  colnames(tmp)<-as.vector(t(outer(dimnames(A)[[1]], 
                                   dimnames(A)[[2]], 
                                   FUN = paste, sep=sep)))
  return(tmp)
}
