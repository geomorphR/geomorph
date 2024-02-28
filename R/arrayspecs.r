#' Convert landmark data matrix into array (p x k x n)
#'
#' Convert a matrix of landmark coordinates into a three-dimensional array 
#'
#' This function converts a matrix of landmark coordinates into a 3D array (p x 
#' k x n), 
#'  which is the required input format for many functions in geomorph. 
#'  The input matrix can be arranged such that the coordinates
#'  of each landmark are found on a separate row, or that each row contains all 
#'  landmark 
#'  coordinates for a single specimen.
#'
#' @param A A matrix containing landmark coordinates for a set of specimens
#' @param p Number of landmarks
#' @param k Number of dimensions (2 or 3)
#' @param sep An optional argument to attempt to separate variable names into 
#' landmark dimension
#' and landmark number variables.  For example, X.1, Y.1, Z.1, X.2, Y.2, Z.2, 
#' ..., can be separated 
#' with sep = ".", such that rows of landmark configurations are labeled 1, 2, 
#' 3, ..., and columns
#' are labeled X, Y, Z. Note, for variables, X1, Y1, Z1, X2, Y2, Z2, ...,
#' where no separator is evident, use sep = "".  Any illogical separation 
#' argument will result in unlabeled
#' variables.  If sep = NULL (the default), unlabeled variables are forced.  
#' This is a good idea if the original matrix
#' has landmarks out of order (as the the landmark labels might not sort as 
#' expected).
#' @export
#' @keywords utilities
#' @author Dean Adams & Mike Collyer
#' @seealso \code{\link{two.d.array}}
#' @return Function returns a 3D array (p x k x n), where p is the number of 
#' landmark points, k is the
#' number of landmark dimensions (2 or 3), and n is the number of specimens. 
#' The third dimension of
#' this array contains names for each specimen if specified in the original 
#' input matrix.
#' @examples 
#' \dontrun{
#' 
#' x <- matrix(rnorm(18), nrow = 3)  # Random triangles (all coordinates on same 
#'  # row for each triangle)
#' arrayspecs(x, 3, 2) 
#'  
#' x2 <- matrix(rnorm(18), ncol = 2) # Random triangles (each landmark on its 
#' # own row)
#' arrayspecs(x2, 3, 2)
#' }

arrayspecs<-function(A, p, k, sep = NULL){  
  if(!is.matrix(A) && !is.data.frame(A)) stop("A must be a data frame or matrix")
  dnames <- dimnames(A)
  n <- length(unlist(A))/(p * k)
  if(k < 2 ) stop("One-dimensional data cannot be used")
  if(all(is.na(match(c(n, n*p), NROW(A))))) stop("Matrix dimensions do not match input")
  specimens <- aperm(array(t(A), c(k,p,n)), c(2,1,3)) 
  dimnames(specimens)[[3]] <- dnames[[1]]
  col.names <- dnames[[2]]
  if(is.null(sep)) col.names <- NULL
  if(!is.null(col.names)){
    if(sep == ""){
      options(warn = -1)
      split.lab <-sort(unique(unlist(strsplit(col.names, split=""))))
      split.lab <- split.lab[length(split.lab)]
      a <- strsplit(col.names, split = split.lab)
      a <- a[sapply(a, length) == 2]
      a <- na.omit(as.numeric(unlist(a)))
      if(length(unlist(a)) == p){
        split.no <- as.character(a[length(a)])
        b <- strsplit(col.names, split = split.no)
        no.char <- sapply(b, nchar)
        dim.tag <- which(no.char == min(no.char))
        b <- unlist(b[dim.tag])
        rn <- a
        cn <- b
      } else rn <- cn <- NULL
      options(warn = 0)
    } else{
      sep = paste("[", noquote(sep), "]")
      a <- strsplit(col.names, split = sep)
      if(length(unlist(a)) == 2*k*p){
        b <- simplify2array(a)
        rn <- unique(b[1,])
        cn <- unique(b[2,])
        if(length(rn) < length(cn)) {tmp <- cn; cn <- rn; rn <- tmp}
      }
      else rn <- cn <- NULL
    }
    dnames2 <- list(rn=rn, cn = cn)
  }
  if(!is.null(sep)) dimnames(specimens)[1:2] <- dnames2
  return(specimens)
}
