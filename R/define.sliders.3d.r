#' Select points to "slide" along three-dimensional curves.
#'
#' An interactive function to define which digitized landmarks of an '*.nts' file
#' will "slide" along three-dimensional (3D) curves.
#' 
#' Function takes a matrix of digitized landmark coordinates, such as made by \code{\link{digit.fixed}}, or \code{\link{buildtemplate}}
#'  and helps user choose which landmarks will be treated as "curve sliders" in Generalized Procrustes analysis
#'  \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks 
#'  (see Bookstein 1997 for algorithm details). 
#'  Each sliding semilandmark ("sliders") will slide between two designated points, along a line 
#'  tangent to the specified curvature.
#' 
#' \subsection{Selection}{ 
#' Choosing which landmarks will be sliders involves landmark selection using a mouse in the rgl plot window. 
#' With a standard 3-button (PC) buildtemplate uses:
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to choose points to be defined as sliders,
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. Macs using platform 
#' specific single button mice, XQuartz must be configured: go to Preferences > Input > tick "Emulate three button mouse":
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing COMMAND key to select points to be defined as sliders,
#'  \item press button while pressing OPTION key to adjust mesh perspective.
#'  \item the mouse SCROLLER or trackpad two finger scroll is used to zoom in an out.
#'  }
#' 
#' To define the sliders, for each sliding landmark along the curve in the format 'before-slider-after':
#' \enumerate{
#'  \item Click to choose the first landmark between which semi-landmark will "slide",
#'  \item Click to choose sliding landmark,
#'  \item Click to choose the last landmark between which semi-landmark will "slide",
#' Screen will show lines connecting the three landmarks, and will highlight the sliding semilandmark in red. 
#' }
#' This procedure is overlapping, so for example a curve defined by a sequence of semilandmarks,
#' the user must select the 2nd point of the first three to be the 1st for the next
#' e.g. 1 2 3 then 2 3 4, etc.
#' }
#' @param spec Name of specimen, as an object matrix containing 3D landmark coordinates
#' @param nsliders Number of landmarks to be semilandmarks that slide along curves
#' @param surfsliders Logical If spec contains landmarks that are "surface sliders",
#' made by \code{\link{buildtemplate}}, "surfslide.csv" should be in working directory 
#' @export
#' @seealso  \code{\link{digit.fixed}}, \code{\link{digitsurface}}, \code{\link{gpagen}}
#' @keywords utilities
#' @author Erik Otarola-Castillo & Emma Sherratt
#' @return Function returns a 'curves x 3' matrix containing the landmark address of the curve sliders, indicating the 
#' points between which the selected point will "slide". Written to the working directory is this matrix as "curveslide.csv".
#' Matrix (or "curveslide.csv") is designed for use by \code{\link{gpagen}} during GPA.
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.
define.sliders.3d<-function(spec, nsliders,surfsliders=FALSE)    {
  spec.name <- deparse(substitute(spec))
  rownames(spec) <- c(1:nrow(spec)) 
  checkmat <- is.matrix(spec)
  if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(spec)[2]
  if (checkdim==2) {stop("Input must be a p-x-k matrix of three-dimensional landmark coordinates") }
  if (surfsliders == TRUE){
    surf <- as.matrix(read.csv("surfslide.csv", header=T))
    spec <- spec[-surf,]
  }
  n <- dim(spec)[1]
  index <- as.numeric(rownames(spec))
  clear3d();ids <- plot3d(spec,size=5,col="black",xlab="x",ylab="y",zlab="z",aspect=FALSE)
  text3d(spec, texts=index,cex=1,adj=c(2,1))  
  curveslide<-NULL 
  for (i in 1:nsliders)      {
    selected<-NULL
    for (j in 1:3)      	{
      keep<-NULL
      keep <- selectpoints3d(ids["data"], value= FALSE, button = "right")[2]
      points3d(spec[keep,1],spec[keep,2],spec[keep,3],size=10,color="red",add=TRUE)
      selected.temp<-cbind(index[keep],spec[,1][keep],spec[,2][keep],spec[,3][keep])
      selected<-rbind(selected,selected.temp)
      if (j==2) {
        points3d(selected[2,2],selected[2,3],selected[2,4],size=10,color="red",add=TRUE) 
      } else {
        points3d(selected[,2],selected[,3],selected[,4],size=7,color="blue",add=TRUE)
      }  
      if (j>=2) {
        lines3d(selected[c(j-1,j),2],selected[c(j-1,j),3],selected[c(j-1,j),4],size=10,color="red",add=TRUE)
      } 
    }
    cat("semi-landmark", paste(selected[2, 1]), "slides between landmarks", 
        paste(selected[1, 1]), "and", paste(selected[3, 1]), "\n")
    curveslide<-rbind(curveslide,selected[,1])
  }
  
  colnames(curveslide)<-c("before","slide","after")
  write.table(curveslide,file="curveslide.csv",row.names=FALSE,col.names=TRUE,sep=",")
  return(curveslide)  
}