#' Rotate a subset of 2D landmarks to common articulation angle
#'
#' A function for rotating a subset of landmarks so that the articulation angle between subsets is constant
#' 
#' This function standardizes the angle between two subsets of landmarks for a set of specimens. The approach assumes a simple
#'  hinge-point articulation between the two subsets, and rotates all specimens such that the angle between landmark subsets 
#'  is equal across specimens (see Adams 1999).  As a default, the mean angle is used, though the user may specify an additional amount by which 
#'  this may be augmented. 
#' 
#' Presently, the function is only implemented for two-dimensional landmark data. 
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @param art.pt A number specifying which landmark is the articulation point between the two landmark subsets
#' @param angle.pts A vector containing numbers specifying which two points used to define the angle (one per subset)
#' @param rot.pts A vector containing numbers specifying which landmarks are in the subset to be rotated
#' @param angle An optional value specifying the additional amount by which the rotation should be augmented (in radians)
#' @param degrees A logical value specifying whether the additional rotation angle is expressed in degrees or radians (radians is default)
#' @author Dean Adams
#' @keywords utilities
#' @return Function returns a (p x k x n) array of landmark coordinates. 
#' @export
#' @references Adams, D. C. 1999. Methods for shape analysis of landmark data from articulated structures. Evolutionary Ecology Research. 1:959-970.
#' @examples
#' #Example using Plethodon
#' #Articulation point is landmark 1, rotate mandibular landmarks (2-5) relative to cranium
#'
#' data(plethspecies) 
#' fixed.angle(plethspecies$land,art.pt=1,angle.pts=c(5,6),rot.pts=c(2,3,4,5))
fixed.angle<-function(A,art.pt=NULL,angle.pts=NULL,rot.pts=NULL,angle=0,degrees = FALSE){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  }
  n<-dim(A)[3];   k<-dim(A)[2];  p<-dim(A)[1]
  if (k!=2){
    stop("Method presently implemented for two-dimensional data only.")}
  if (degrees == FALSE & angle>pi*2){
    stop("Angle was specified in degrees (use degrees = T).")}
  if (degrees == FALSE & angle< -pi*2){
    stop("Angle was specified in degrees (use degrees = T).")}
  angl<-array(NA,dim=n)
  for (i in 1:n){   
    A[,,i]<-t(t(A[,,i])-A[art.pt,,i])  
    angl[i]<-  acos((A[angle.pts[1],,i]/sqrt(sum(A[angle.pts[1],,i]^2)))%*%
                      (A[angle.pts[2],,i]/sqrt(sum(A[angle.pts[2],,i]^2))))
  }  
  if(degrees ==TRUE)  angle = angle*pi/180
  dev.angle<- (angl-mean(angl))+angle
  if(A[angle.pts[1],1,1]<0){dev.angle<- -1*dev.angle}
  for (i in 1:n){   
    r = matrix(c(cos(dev.angle[i]),-sin(dev.angle[i]),
                 sin(dev.angle[i]),cos(dev.angle[i])),2) 
    A[rot.pts,,i] = A[rot.pts,,i]%*%r 
  }  
  return(A)
}