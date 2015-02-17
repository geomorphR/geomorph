#' Plot shape differences between a reference and target specimen
#'
#' Function plots shape differences between a reference and target specimen
#'
#' The function generates a plot of the shape differences of a target specimen relative to a reference 
#'  specimen. The option {mag} allows the user to indicates the degree of magnification to be used when 
#'  displaying the shape difference. The function will plot either two- or three-dimensional data. For
#'  two-dimensional data and thin-plate spline deformation plots, the user may also supply boundary curves 
#'  of the object, which will be deformed from 
#'  the reference to the target specimen using the thin-plate spline. Such curves are often useful in describing 
#'  the biological shape differences expressed in the landmark coordinates.  Note that to utilize this option, 
#'  a boundary curve from a representative specimen must first be warped to the reference specimen using
#'   \code{\link{warpRefOutline}}.
#'   
#'  Four distinct methods for plots are available:
#'  \enumerate{
#'  \item {TPS} a thin-plate spline deformation grid is generated. For 3D data, 
#'  this method will generate thin-plate spline deformations in the x-y and x-z planes. A boundary curve 
#'  will also be deformed if provided by the user.
#'  \item {vector}: a plot showing the vector displacements between corresponding landmarks in the reference 
#'  and target specimen is shown. 
#'  \item {points} a plot is displayed with the landmarks in the target (black) 
#'  overlaying those of the reference (gray). Additionally, if a matrix of links is provided, the 
#'  landmarks of the mean shape will be connected by lines.  The link matrix is an M x 2 matrix, where 
#'  M is the desired number of links. Each row of the link matrix designates the two landmarks to be 
#'  connected by that link. 
#'  \item {surface} a mesh3d surface is warped using thin-plate spline (for 3D data only). 
#'  Requires mesh3d object in option {mesh}, made using \code{\link{warpRefMesh}}. 
#'  }
#'  This function combines numerous plotting functions found in Claude (2008).
#'
#' @param M1 Matrix of landmark coordinates for the first (reference) specimen
#' @param M2 Matrix of landmark coordinates for the second (target) specimen
#' @param mesh A mesh3d object for use with {method="surface"}
#' @param outline An x,y curve or curves warped to the reference (2D only)
#' @param method Method used to visualize shape difference; see below for details
#' @param mag The desired magnification to be used when visualizing the shape difference (e.g., mag=2)
#' @param links An optional matrix defining for links between landmarks
#' @param label A logical value indicating whether landmark numbers will be plotted 
#' @param ... Parameters to be passed to \code{\link{plot}}, \code{\link{plot3d}} or \code{\link{shade3d}}
#' @return If using {method="surface"}, function will return the warped mesh3d object.
#' @keywords visualization
#' @export
#' @author Dean Adams & Emma Sherratt
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York. 
#' @seealso  \code{\link{define.links}}
#' @seealso  \code{\link{warpRefMesh}}
#' @seealso  \code{\link{warpRefOutline}}
#' @examples
#' 
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' ref<-mshape(Y.gpa$coords)
#' plotRefToTarget(ref,Y.gpa$coords[,,39])
#' plotRefToTarget(ref,Y.gpa$coords[,,39],mag=2,outline=plethodon$outline)   #magnify difference by 2X
#' plotRefToTarget(ref,Y.gpa$coords[,,39],method="vector",mag=3)
#' plotRefToTarget(ref,Y.gpa$coords[,,39],method="points",outline=plethodon$outline)
#' 
#'
#' # Three dimensional data
#' 
#' data(scallops)
#' Y.gpa<-gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
#' ref<-mshape(Y.gpa$coords)
#' plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
plotRefToTarget<-function(M1,M2,mesh= NULL,outline=NULL,method=c("TPS","vector","points","surface"),
                          mag=1.0,links=NULL, label=FALSE,...){
  method <- match.arg(method)
  if(any(is.na(M1))==T){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(any(is.na(M2))==T){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  k<-dim(M1)[2]
  mag<-(mag-1)
  M2<-M2+(M2-M1)*mag
  limits = function(x,s){ 
    r = range(x)
    rc=scale(r,scale=F)
    l=mean(r)+s*rc
  }
  if(k==2){
    if(method=="TPS"){
      tps(M1,M2,20)
      if(label == TRUE){text(M2, label = paste(1:dim(M2)[1]), adj = 0.5, pos = 1, cex=0.8)}
      if(!is.null(outline)){
        curve.warp <- tps2d(outline, M1, M2)
      }
      if(!is.null(outline)){
        points(curve.warp,pch=19, cex=0.1) 
      }
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=2)
        }
      }
    }
    if(method=="vector"){
      plot(M1,asp=1,pch=21,bg="gray",xlim=limits(M1[,1],1.25),ylim=limits(M1[,2],1.25),cex=1,xlab="x",ylab="y",...)
      if(label == TRUE){text(M1, label = paste(1:dim(M1)[1]), adj = 0.5, pos = 1, cex=0.8)}
      arrows(M1[,1],M1[,2],M2[,1],M2[,2],length=0.075,lwd=2)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1,lty=2)
        }
      }
    }
    if(method=="points"){
      plot(M1,asp=1,pch=21,bg="gray",xlim=limits(M1[,1],1.25),ylim=limits(M1[,2],1.25),cex=1,xlab="x",ylab="y",...)
      if(label == TRUE){text(M1, label = paste(1:dim(M1)[1]), adj = 0.5, pos = 1, col = "gray", cex=0.8)}
      points(M2,asp=1,pch=21,bg="black",cex=1)
      if(!is.null(outline)){
        curve.warp <- tps2d(outline, M1, M2)
      }
      if(!is.null(outline)){
        points(outline,pch=19, cex=0.1,col="gray") 
        points(curve.warp,pch=19, cex=0.1,col="black") 
      }
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M1[links[i,1],1],M1[links[i,1],2],M1[links[i,2],1],M1[links[i,2],2],lwd=1,lty=2, col= "gray")
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1,lty=2, col= "black")
        }
      }
    }
    if(method=="surface"){
      stop("Surface plotting for 3D landmarks only.")
    }      
  }
  if(k==3){
    if(method=="TPS"){
      old.par <- par(no.readonly = TRUE)
      layout(matrix(c(1,2),1,2))
      par(mar=c(1,1,1,1))
      tps(M1[,1:2],M2[,1:2],20)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1)
        }
      }
      title("X,Y tps grid")
      b<-c(1,3)
      tps(M1[,b],M2[,b],20)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],3],M2[links[i,2],1],M2[links[i,2],3],lwd=1)
        }
      }
      title("Y,Z tps grid")
      layout(1)
      on.exit(par(old.par))
    }
    if(method=="vector"){
      plot3d(M1,type="s",col="gray",size=1.25,aspect=FALSE,...)
      if(label == TRUE){text3d(M1, texts = paste(1:dim(M1)[1]), adj = 1.3, pos = 4, cex=0.8)}
      for (i in 1:nrow(M1)){
        segments3d(rbind(M1[i,],M2[i,]),lwd=2)
      }
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),lwd=.75)
        }
      }
    }
    if(method=="points"){
      plot3d(M1,type="s",col="gray",size=1.25,aspect=FALSE,...)
      if(label == TRUE){text3d(M1, texts = paste(1:dim(M1)[1]), adj = 1.3, pos = 4, cex=0.8)}
      points3d(M2,color="black",size=5)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),lwd=.75)
        }
      }
    }
    if (method == "surface") {
      if(is.null(mesh)==T){
        stop("Surface plotting requires a template mesh3d object (see 'warpRefMesh').")
      }
      warp.PLY <- mesh
      vb <- as.matrix(t(mesh$vb)[,-4])
      warp <- tps2d3d(vb, M1, M2)
      warp.PLY$vb <- rbind(t(warp), 1)
      open3d(); shade3d(warp.PLY, ...)
      if(label == TRUE){text3d(M1, texts = paste(1:dim(M1)[1]), adj = 1.3, pos = 4, cex=0.8)}
      return(warp.PLY)
    }
  }
}