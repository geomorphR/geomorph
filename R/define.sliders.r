#' Select points to "slide" along curves
#'
#' An interactive function to define which landmarks will "slide" along curves.
#'
#' Function takes a matrix of digitized landmark coordinates, such as made by \code{\link{digitize2d}},
#'  and helps user choose which landmarks will be treated as "sliders" in Generalized Procrustes analysis
#'  \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks 
#'  (see Bookstein 1997 for algorithm details). 
#'  Each sliding semilandmark ("sliders") will slide between two designated points, along a line 
#'  tangent to the specified curvature.
#'  
#' Defining landmarks is an interactive procedure (see below for 2D routines). The procedure is overlapping. 
#' For example: there are 5 landmarks (1:5), 1 and 5 are landmarks and 2,3,4 are sliders,
#' the user must select '1' '2' '3', and then '2' '3' '4', and then '3' '4' '5'.
#' 
#' \subsection{Selection in 2D}{ 
#' Choosing which landmarks will be sliders involves landmark selection using a mouse in the plot window. 
#' To define the sliders, for each sliding landmark along the curve in the format 'before-slider-after',
#' using the LEFT mouse button (or regular button for Mac users), click on the hollow circle to choose the landmark
#' in the following order:
#' \enumerate{
#'  \item Click to choose the first landmark between which semi-landmark will "slide",
#'  \item Click to choose sliding landmark,
#'  \item Click to choose the last landmark between which semi-landmark will "slide",
#'  Selected landmarks will be filled in and lines are drawn connecting the three landmarks, 
#'  and will highlight the sliding semilandmark in red and the flanking landmarks in blue. 
#' } }
#' 
#'  \subsection{Notes for geomorph 4.1}{ 
#'  Starting with geomorph version 4.1, interactive module selection in 3D is no longer
#'  supported, as RGL is not supported under MacOS Tahoe.  Only AUTO mode is available.
#'  }
#'  
#' \subsection{AUTO mode}{ 
#' The input 'landmarks' can be simply a vector of numbers corresponding to the "sliders" (semilandmarks) in the order they appear along a curve on the specimen. This can be made by c() or seq() or any other reasonable method.
#'  
#'  If the sliders form a closed curve, then the function assumes that the first and last landmarks in the 'landmarks' vector are THE SAME are fixed (not sliders). e.g. if landmark 1 is a fixed landmark, and 2, 3 and 4 are semilandmarks, then sliders = c(1,2,3,4,1).
#   If the sliders form an open curve, then the function assumes the first and last landmarks are DIFFERENT and are fixed (not sliders), e.g. if landmark 1 and 5 are fixed landmarks, and 2, 3 and 4 are semilandmarks, then sliders = c(1,2,3,4,5).
#' }
#' 
#' @param landmarks A matrix containing landmark coordinates of landmarks and semilandmarks, OR A vector containing a sequence of numbers corresponding to the landmarks in the order they appear along the curve (for AUTO mode)
#' @param nsliders Number of landmarks to be semilandmarks that slide along curves
#' @param surfsliders (3D only) If 'landmarks' contains "surface sliders",
#'  these should be given as a vector or use surfsliders = T, and function looks for "surfslide.csv" in working directory.
#' @param write.file A logical value indicating whether the matrix is written to file as .csv.
#' @return Function returns a 'nsliders-x-3' matrix containing the landmark address of the curve sliders, indicating the landmarks between which the slider landmarks will "slide". If write.file = T the matrix is also written to working directory as "curveslide.csv". Matrix (or "curveslide.csv") is designed for use by \code{\link{gpagen}} during GPA.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{digitize2d}}, \code{\link{gpagen}}, \code{\link{digit.curves}}
#' @author Emma Sherratt & Dean Adams 
#' @examples  
#' ## (not run)
#' ## Examples of AUTO mode 
#'  # data(scallops)
#'  ## 1 curve of sliding semilandmark
#'  # Define sliders for scallopdata
#'  #sliders = define.sliders(c(5:16,1))
#' 
#'  ## 2 curves of sliding semilandmarks
#'  # Define sliders for 10 landmarks, where LMs 1, 5, and 10 fixed
#'  # 2, 3, and 4 are along a curve between 1 and 5
#'  # and 6, 7, 8, and 9 are along a curve between 5 and 10.
#'  #sliders = rbind(define.sliders(1:5), define.sliders(5:10)) 
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.

define.sliders <- function(landmarks, 
                         nsliders, 
                         surfsliders = NULL, 
                         write.file = TRUE){
  checkmat <- is.matrix(landmarks)
  if (checkmat==FALSE) { 
    if(length(dim(landmarks) == 3)){ spec <- as.matrix(landmarks[,,1]) }
    # Auto Mode
    if(is.null(dim(landmarks))){
      nsliders <- length(landmarks)
      CV <- matrix(NA, ncol=3, nrow=nsliders-2)
      for (i in 1:(nsliders-2)){
        CV[i,] <- landmarks[1:3]
        landmarks <- landmarks[-1] }
      colnames(CV)<-c("before","slide","after")
      if(write.file == TRUE){write.table(CV,file="curveslide.csv",row.names=FALSE,col.names=TRUE,sep=",")}
      return(CV)
    }
  }
  if (checkmat == TRUE) { spec <- landmarks }
  checkdim <- dim(spec)[2]
  # 2D interactive routine
  if (checkdim==2) {
    plot(spec[,1],spec[,2],cex=1,pch=21,bg="white",xlim=range(spec[,1]),ylim=range(spec[,2]),asp=1)
    text(spec[,1],spec[,2],label=paste(1:dim(spec)[1]),adj=.5,pos=4)
    selected<-matrix(NA,ncol=3,nrow=nsliders)
    select<-NULL
    for(i in 1:nsliders){
      for(j in 1:3){
        select<-identify(spec,n=1,plot=FALSE,cex=5,pch=25)
        selected[i,j]<-select
        if(j==2){
          points(spec[select,][1],spec[select,][2],cex=1.5,pch=19,col="red")
          arrows(spec[selected[i,j],][1],spec[selected[i,j],][2],spec[selected[i,j-1],][1],spec[selected[i,j-1],][2]
                 ,col="red",lwd=2,length=.15)
        } else {
          points(spec[select,][1],spec[select,][2],cex=1.1,pch=19,col="blue")
        }
        if(j==3){
          arrows(spec[selected[i,j],][1],spec[selected[i,j],][2],spec[selected[i,j-1],][1],spec[selected[i,j-1],][2]
                 ,col="red",lwd=2,length=.15,code=1)
        } else { NA
        } 
      }
    }
    colnames(selected)<-c("before","slide","after")
    if(write.file == TRUE){write.table(selected,file="curveslide.csv",row.names=FALSE,col.names=TRUE,sep=",")}
    return(selected)
  }
}
