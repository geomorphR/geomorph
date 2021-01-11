#' Rotate a subset of 2D landmarks to common articulation angle
#'
#' A function for rotating a subset of landmarks so that the articulation angle between subsets is constant
#' 
#' This function standardizes the angle between two subsets of landmarks for a set of specimens. The approach assumes a simple
#'  hinge-point articulation between the two subsets, and rotates all specimens such that the angle between landmark subsets 
#'  is equal across specimens (see Adams 1999).  As a default, the mean angle is used, though the user may specify an additional amount by which 
#'  this may be augmented. To quantify the angle, users may specify a single landmark in each subset as angle endpoints, or may specify a set of landmarks.
#'  If the latter, the centroid of those points is used. 
#' 
#' Presently, the function is only implemented for two-dimensional landmark data. 
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param art.pt A number specifying which landmark is the articulation point between the two landmark subsets
#' @param angle.pts.1 A vector or single value specifying the angle point of one subset.  If more that one value
#' is provided, the centroid of the landmarks described by the vector will be used; a single value
#' identifies a specific landmark to use.  
#' @param angle.pts.2 A vector or single value specifying the angle point of the second subset. This could be 
#' the entire set of points of an articulated structure to be rotated.
#' @param rot.pts A vector containing numbers specifying which landmarks are in the subset to be rotated.  If NULL,
#' it is assumed that the points to be rotated are the same as those in angle.pts.2.
#' @param angle An optional value specifying the additional amount by which the rotation should be augmented (in radians).
#' It might be essential to use a negative angle if centroids from multiple points are used for angle points.  It should be 
#' clear if this is the case, upon plotting results.
#' @param degrees A logical value specifying whether the additional rotation angle is expressed in degrees or radians (radians is default)
#' @author Dean Adams and Michael Collyer
#' @keywords utilities
#' @return Function returns a (p x k x n) array of landmark coordinates. 
#' @export
#' @references Adams, D. C. 1999. Methods for shape analysis of landmark data from articulated structures. Evolutionary Ecology Research. 1:959-970.
#' @examples
#' #Example using Plethodon
#' #Articulation point is landmark 1, rotate mandibular landmarks (2-5) 
#' # relative to cranium
#'
#' data(plethspecies) 
#' # Using specific points:
#' newLM1 <- fixed.angle(plethspecies$land,
#' art.pt=1, angle.pts.1 = 5, 
#' angle.pts.2 = 6, rot.pts = c(2,3,4,5))
#' Y.gpa1 <- gpagen(newLM1)
#' plot(Y.gpa1, mean = FALSE)
#' 
#' # Using centroids from subsets
#' newLM2 <- fixed.angle(plethspecies$land,art.pt=1, 
#' angle.pts.1 = c(1, 6:11), 
#' angle.pts.2 = 2:5, 
#' rot.pts = NULL, angle = 20, 
#' degrees = TRUE) # rotated points same as second partition
#' Y.gpa2 <- gpagen(newLM2)
#' plot(Y.gpa2, mean = FALSE)
#' 
fixed.angle<-function(A, art.pt=NULL, angle.pts.1, angle.pts.2, 
                      rot.pts=NULL, angle=0, degrees = FALSE){
  if (length(dim(A)) != 3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  }
 
   n <- dim(A)[3];   k <- dim(A)[2];  p <- dim(A)[1]
  
   if(is.null(art.pt)) 
    stop("An articulation point must be provided.")
   
  if(is.null(angle.pts.1) || is.null(angle.pts.2))
    stop("At least one point must be provided for set of angle points.")
  
  if(is.null(rot.pts)) rot.pts <- angle.pts.2
  
   if(length(angle.pts.1) > p-1 || length(angle.pts.1) > p-1)
    stop("Too many points from which to calculate angle points")
  
   if(length(intersect(angle.pts.1, angle.pts.2) > 0))
    stop("Points in each subset must be unique (one or more points are shared between subsets)")
  
   if(k != 2)
    stop("Method presently implemented for two-dimensional data only.")
  
   if(degrees == FALSE & angle > pi*2)
    stop("Angle was specified in degrees (use degrees = T).")
  
   if(degrees == FALSE & angle < -pi*2)
    stop("Angle was specified in degrees (use degrees = T).")
  
  centroid1 <- ifelse(length(angle.pts.1) > 1, TRUE, FALSE)  
  centroid2 <- ifelse(length(angle.pts.2) > 1, TRUE, FALSE)  
  angl<-array(NA, dim=n)
  
  angle.pts <- function(y, ap1, ap2, artp){
    ap <- y[artp,]
    y <- y - matrix(ap, p, k, byrow = TRUE)
    ay1 <- y[ap1,]
    ay2 <- y[ap2,]
    if(centroid1) ay1 <- colMeans(ay1)
    if(centroid2) ay2 <- colMeans(ay2)
    list(y = y, ap = ap, ay1 = ay1, ay2 = ay2)
  }
  
  A.list <- lapply(1:n, function(j){
    a <- A[,,j]
    angle.pts(a, angle.pts.1, angle.pts.2, art.pt)
  })
  
  options(warn = -1)
  angl<- sapply(1:n, function(j){
    a <- A.list[[j]]
    acos(crossprod(a$ay1/sqrt(sum(a$ay1^2)),
                   a$ay2/sqrt(sum(a$ay2^2))))
  })
  options(warn = 0)
  
  if(any(is.na(angl))) {
    cat("\nThere is a problem with the choice of articulated subsets.",
        "\nThe vectors defining articulation angles switch positions in",
        "\nat least one specimen.\n")
    cat("\nThese specimens appear to have this issue:\n")
    print(which(is.na(angl)))
    cat("\n")
    stop("\nThe specimens above will have to be redigitized or removed,", 
         "\nor different angle points will have to be chosen.", call. = FALSE)
  }
  
  if(degrees)  angle = angle*pi/180
  dev.angle <- (angl - mean(angl)) + angle
  
  check <- mean(A[angle.pts.1,1,1])
  if(check < 0){dev.angle <- -1 * dev.angle}
  
  for (i in 1:n){   
    r <- matrix(c(cos(dev.angle[i]),-sin(dev.angle[i]),
                 sin(dev.angle[i]),cos(dev.angle[i])),2, 2)
    y <- A.list[[i]]$y
    A[,,i] <- y
    A[rot.pts,,i] = y[rot.pts,] %*% r
  }  
  return(A)
}
