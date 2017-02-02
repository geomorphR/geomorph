#' Convert landmark data matrix into array (p x k x n)
#'
#' Convert a matrix of landmark coordinates into a three-dimensional array 
#'
#' This function converts a matrix of landmark coordinates into a 3D array (p x k x n), 
#'  which is the required input format for many functions in geomorph. 
#'  The input matrix can be arranged such that the coordinates
#'  of each landmark are found on a separate row, or that each row contains all landmark 
#'  coordinates for a single specimen.
#'
#' @param A A matrix containing landmark coordinates for a set of specimens
#' @param p Number of landmarks
#' @param k Number of dimensions (2 or 3)
#' @export
#' @keywords utilities
#' @author Dean Adams & Mike Collyer
#' @seealso \code{\link{two.d.array}}
#' @return Function returns a 3D array (p x k x n), where p is the number of landmark points, k is the
#' number of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of
#' this array contains names for each specimen if specified in the original input matrix.
#' @examples 
#' x<-matrix(rnorm(18),nrow=3)  # Random triangles (all coordinates on same row for each triangle)
#' arrayspecs(x,3,2) 
#'  
#' x2<-matrix(rnorm(18),ncol=2) # Random triangles (each landmark on its own row)
#' arrayspecs(x2,3,2)
arrayspecs<-function(A,p,k){  
  if(!is.matrix(A) && !is.data.frame(A)) stop("A must be a data frame or matrix")
  dnames <- dimnames(A)
  n <- length(unlist(A))/(p * k)
  if(k < 2 ) stop("One-dimensional data cannot be used")
  if(n != NROW(A)) stop("Matrix dimensions do not match input")
  specimens <- aperm(array(t(A), c(k,p,n)), c(2,1,3)) 
  dimnames(specimens)[[3]] <- dnames[[1]]
  col.names <- dnames[[2]]
  if(!is.null(col.names)){
    a <- strsplit(col.names, split = "[.]")
    if(length(unlist(a)) == 2*k*p){
      b <- simplify2array(a)
      rn <- unique(b[1,])
      cn <- unique(b[2,])
    }  else b <- rn <- cn <- NULL
    dnames2 <- list(rn=rn, cn = cn)
  }
  dimnames(specimens)[1:2] <- dnames2
  return(specimens)
}