#' Select points to "slide" along two-dimensional curves.
#'
#' An interactive function to define which landmarks will "slide" along two-dimensional (2D) curves.
#'
#' Function takes a matrix of digitized landmark coordinates, such as made by \code{\link{digitize2d}},
#'  and helps user choose which landmarks will be treated as "curve sliders" in Generalized Procrustes analysis
#'  \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks 
#'  (see Bookstein 1997 for algorithm details). 
#'  Each sliding semilandmark ("sliders") will slide between two designated points, along a line 
#'  tangent to the specified curvature.
#' 
#' \subsection{Selection}{ 
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
#' }
#' This procedure is overlapping, so for example, a curve defined by a sequence of semilandmarks,
#' the user must select the 2nd point of the first three to be the 1st for the next
#' e.g. 1 2 3 then 2 3 4, then 3 4 5 etc.
#' }
#' @param spec Name of specimen, as an object matrix containing 2D landmark coordinates
#' @param nsliders Number of landmarks to be semilandmarks that slide along curves
#' @return Function returns a 'nsliders-x-3' matrix containing the landmark address of the curve sliders, indicating the 
#' landmarks between which the slider landmarks will "slide". The matrix is also written to working directory
#' as "curveslide.csv". Matrix (or "curveslide.csv") is designed for use by \code{\link{gpagen}} during GPA.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{digitize2d}}, \code{\link{gpagen}}
#' @author Dean Adams, Erik Otarola-Castillo, Emma Sherratt
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.
define.sliders.2d<-function(spec, nsliders){
  spec.name <- deparse(substitute(spec))
  checkmat <- is.matrix(spec)
  if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(spec)[2]
  if (checkdim==3) {stop("Input must be a p-x-k matrix of two-dimensional landmark coordinates") }
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
  write.table(selected,file="curveslide.csv",row.names=FALSE,col.names=TRUE,sep=",")
  return(selected)
}