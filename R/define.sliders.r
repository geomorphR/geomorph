#' Select points to "slide" along curves
#'
#' An interactive function to define which landmarks will "slide" along two-dimensional (2D) or
#' three-dimensional (3D) curves.
#'
#' Function takes a matrix of digitized landmark coordinates, such as made by \code{\link{digitize2d}} or \code{\link{digit.fixed}},
#'  and helps user choose which landmarks will be treated as "sliders" in Generalized Procrustes analysis
#'  \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks 
#'  (see Bookstein 1997 for algorithm details). 
#'  Each sliding semilandmark ("sliders") will slide between two designated points, along a line 
#'  tangent to the specified curvature.
#'  
#' Defining landmarks is an interactive procedure (see below for 2D and 3D routines). The procedure is overlapping. 
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
#' \subsection{Selection in 3D}{ 
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
#' } }
#' 
#' \subsection{AUTO mode}{ 
#' The input 'landmarks' can be simply a vector of numbers corresponding to the "sliders" (semilandmarks) in the order they appear along a curve on the specimen. This can be made by c() or seq() or any other reasonable method.
#'  
#'  If the sliders form a closed curve, then the function assumes that the first and last landmarks in the 'landmarks' vector are THE SAME are fixed (not sliders). e.g. if landmark 1 is a fixed landmark, and 2, 3 and 4 are semilandmarks, then sliders = c(1,2,3,4,1).
#   If the sliders form an open curve, then the function assumes the first and last landmarks are DIFFERENT and are fixed (not sliders), e.g. if landmark 1 and 5 are fixed landmarks, and 2, 3 and 4 are semilandmarks, then sliders = c(1,2,3,4,5).
#' }
#' 
#' @param landmarks A matrix containing 2D or 3D landmark coordinates of landmarks and semilandmarks, OR A vector containing a sequence of numbers corresponding to the landmarks in the order they appear along the curve (for AUTO mode)
#' @param nsliders Number of landmarks to be semilandmarks that slide along curves
#' @param surfsliders (3D only) If 'landmarks' contains "surface sliders",
#' e.g. made by \code{\link{buildtemplate}}, these should be given as a vector or use surfsliders = T, and function looks for "surfslide.csv" in working directory.
#' @param write.file A logical value indicating whether the matrix is written to file as .csv.
#' @return Function returns a 'nsliders-x-3' matrix containing the landmark address of the curve sliders, indicating the landmarks between which the slider landmarks will "slide". If write.file = T the matrix is also written to working directory as "curveslide.csv". Matrix (or "curveslide.csv") is designed for use by \code{\link{gpagen}} during GPA.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{digitize2d}}, \code{\link{digit.fixed}}, \code{\link{gpagen}}, \code{\link{digit.curves}}
#' @author Emma Sherratt, Dean Adams, Erik Otarola-Castillo 
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.

define.sliders<-function(landmarks, nsliders, surfsliders=NULL, write.file = TRUE) {
  checkmat <- is.matrix(landmarks)
  if (checkmat==FALSE) { 
    if(length(dim(landmarks)) == 3){stop("'landmarks' should be the shape matrix of a single specimen") }
    nsliders <- length(landmarks)
    CV <- matrix(NA, ncol=3, nrow=nsliders-2)
    for (i in 1:(nsliders-2)){
      CV[i,] <- landmarks[1:3]
      landmarks <- landmarks[-1] }
    colnames(CV)<-c("before","slide","after")
    if(write.file == TRUE){write.table(CV,file="curveslide.csv",row.names=FALSE,col.names=TRUE,sep=",")}
    return(CV)
  }
  if (checkmat == TRUE) { spec <- landmarks }
  checkdim <- dim(spec)[2]
  # 2D routine
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
  # 3D routine
  if (checkdim==3) {
    rownames(spec) <- c(1:nrow(spec)) 
    if(!is.null(surfsliders)){
    if(is.logical(surfsliders)){
      surf <- as.matrix(read.csv("surfslide.csv", header=T))
      spec <- spec[-surf,]} else     
    spec <- spec[-surfsliders,]
    }
    n <- dim(spec)[1]
    index <- as.numeric(rownames(spec))
    clear3d();ids <- plot3d(spec,size=5,col="black",xlab="x",ylab="y",zlab="z",aspect=FALSE)
    text3d(spec, texts=index,cex=1,adj=c(2,1))  
    curveslide<-NULL 
    for (i in 1:nsliders)      {
      selected<-NULL
      for (j in 1:3)        {
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
    if(write.file == TRUE){write.table(curveslide,file="curveslide.csv",row.names=FALSE,col.names=TRUE,sep=",")}
    return(curveslide)  
}
}