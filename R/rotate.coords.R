#' Rotate or flip landmark or coordinate configurations
#'
#' Function to rotate or flip 2D landmark or coordinate configurations to re-orientate them, as desired. 
#'
#' This function will allow a user to rotate or flip one or several configurations of raw landmarks or
#' Procrustes residuals from GPA, as desired.  This function only works with two-dimensional configurations.
#' This is a useful tool if importing coordinates or performing GPA produces undesired orientations (such as flipping
#' configurations upside down, when aligning them to their PCs).  It is not as useful with 3D coordinates, as the plotting tools 
#' for 3D coordinates already have built-in rotation capabilities.
#' 
#' The function returns an object of the same class as input, after having rotated or flipped
#' the coordinates of each specimen.  If multiple steps are required, the function
#' can be used in a recursive fashion.
#' 
#' @param A One of either an array of landmark coordinates (p, 2, n dimensions for n specimens and p 2D points),
#' a class gpagen object, or a p x 2 matrix for a single specimen.  
#' @param type The type of rotation or flip to be performed.  Specimens can be flipped with respect to x or y axes,
#' or rotated clockwise (C) or counter-clockwise (CC).  
#' @param index An index to indicate which specimens should be rotated or flipped.  If NULL (default)
#' all specimens are rotated.  A binary index (0 = do not rotate; 1 = rotate) as a vector with the same length as 
#' the number of specimens will direct which specimens are manipulated.  A factor with two levels or a numeric vector with 
#' two any two values can also be used.  
#' The second level (or larger value) will be the level manipulated.  For example, a factor with levels = c("a", "b") 
#' will rotate specimens matching level "b".  This function might be useful for reflecting specimens with bilateral structures.
#' @export
#' @keywords utilities
#' @author Michael Collyer
#' 
#' @examples
#' \dontrun{
#' data(plethodon)
#' Y.gpa <- gpagen(plethodon$land)
#' plot(Y.gpa)
#' Y.gpa2 <- rotate.coords(Y.gpa, "flipX")
#' plot(Y.gpa2)
#' Y.gpa3 <- rotate.coords(Y.gpa2, "rotateCC")
#' plot(Y.gpa3)
#' 
#' spec1 <- Y.gpa$coords[,,1]
#' plot(spec1, asp = 1)
#' spec1 <- rotate.coords(spec1, "flipY")
#' plot(spec1, asp = 1)
#' 
#' specs1to3 <- Y.gpa$coords[,,1:3]
#' plotAllSpecimens(specs1to3)
#' specs1to3 <- rotate.coords(specs1to3, "rotateC")
#' plotAllSpecimens(specs1to3)
#' }
rotate.coords <- function(A, type = c("flipX", "flipY", "rotateC", "rotateCC"),
                          index = NULL) {
  if(inherits(A, "gpagen")) {
    obj <- "gpagen" 
    n <- length(A$Csize)
  } else

  if(is.array(A)) obj <- "array" else 
      stop("A must be an array or a class gpagen object.")
    
    if(obj == "array") {
      dims <- dim(A)
      check <- intersect(length(dims), c(2, 3))
      if(length(check) == 0)
        stop("Either a matrix with two columns or an array with three dimensions is required.")
      if(length(dims) == 2) {
        obj <- "matrix"
        n <- 1
      } else n <- dims[[3]]
    }
    
   type <- match.arg(type)
   if(type == "flipX") rot <- matrix(c(-1, 0, 0, 1), 2, 2)
   if(type == "flipY") rot <- matrix(c(1, 0, 0, -1), 2, 2)
   if(type == "rotateC") rot <- matrix(c(0, 1, -1, 0), 2, 2) 
   if(type == "rotateCC") rot <- matrix(c(0, -1, 1, 0), 2, 2) 
   
   if(!is.null(index)){
     if(is.factor(index)) index <- as.numeric(index) - 1
     if(length(index) != n)
       stop("\nThe index must be the same length as the number of specimens.\n")
     ui <- unique(index)
     sdiff <- setdiff(ui, c(0, 1))
     if(length(sdiff) > 1) {
       if(length(ui) != 2)
         stop("\nThe index is illogical.  It must have two levels only.\n")
       index <- floor(index/max(index))
       ui <- unique(index)
       sdiff <- setdiff(ui, c(0, 1))
       if(length(sdiff) > 1) 
         stop("\nThe index cannot be coerced into a binary vector.\n")
     }
   } else index <- rep(1, n)
   
   if(obj == "gpagen") {
     id <- dimnames(A$coords)
     Y <- A$coords
     Y <- lapply(1:(dim(Y)[[3]]), function(j) as.matrix(Y[,,j]))
     Y <- lapply(1:length(Y), function(j){
       if(index[j] == 1) Y[[j]] %*% rot else Y[[j]]
     })
     A$coords <- simplify2array(Y)
     dimnames(A$coords) <- id
   }
   
   if(obj == "array") {
     id <- dimnames(A)
     Y <- A
     Y <- lapply(1:(dim(Y)[[3]]), function(j) as.matrix(Y[,,j]))
     Y <- lapply(1:length(Y), function(j){
       if(index[j] == 1) Y[[j]] %*% rot else Y[[j]]
     })
     A <- simplify2array(Y)
     dimnames(A) <- id
   }
   
   if(obj == "matrix") {
     id <- dimnames(A)   
     A <- A %*% rot
     dimnames(A) <- id
   }
   return(A)
}
