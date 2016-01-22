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
#' @param gridPars An optional object made by \code{\link{gridPar}}
#' @param useRefPts An option (logical value) to use reference configuration points rather than target configuration points (when {method = "TPS"})
#' @param ... Additional parameters not covered by \code{\link{gridPar}} to be passed to \code{\link{plot}}, \code{\link{plot3d}} or \code{\link{shade3d}}
#' @return If using {method="surface"}, function will return the warped mesh3d object.
#' @keywords visualization
#' @export
#' @author Dean Adams, Emma Sherratt & Michael Collyer
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York. 
#' @seealso  \code{\link{gridPar}}
#' @seealso  \code{\link{define.links}}
#' @seealso  \code{\link{warpRefMesh}}
#' @seealso  \code{\link{warpRefOutline}}
#' @examples
#' # Two dimensional data
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' ref<-mshape(Y.gpa$coords)
#' plotRefToTarget(ref,Y.gpa$coords[,,39])
#' plotRefToTarget(ref,Y.gpa$coords[,,39],mag=2,outline=plethodon$outline)   #magnify by 2X
#' plotRefToTarget(ref,Y.gpa$coords[,,39],method="vector",mag=3)
#' plotRefToTarget(ref,Y.gpa$coords[,,39],method="points",outline=plethodon$outline)
#' plotRefToTarget(ref,Y.gpa$coords[,,39],gridPars=gridPar(pt.bg = "green", pt.size = 1),
#' method="vector",mag=3)
#'
#' # Three dimensional data
#' # data(scallops)
#' # Y.gpa<-gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
#' # ref<-mshape(Y.gpa$coords)
#' # plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
#' # scallinks <- matrix(c(1,rep(2:16, each=2),1), nrow=16, byrow=TRUE)
#' # plotRefToTarget(ref,Y.gpa$coords[,,1],gridPars=gridPar(tar.pt.bg = "blue", tar.link.col="blue",
#' # tar.link.lwd=2), method="points", links = scallinks)
#' 
plotRefToTarget<-function(M1,M2,mesh= NULL,outline=NULL,method=c("TPS","vector","points","surface"),
                          mag=1.0,links=NULL, label=FALSE, gridPars = NULL, useRefPts=FALSE,...){
  method <- match.arg(method)
  if(any(is.na(M1))==TRUE){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(any(is.na(M2))==TRUE){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(is.null(gridPars)) gP = gridPar() else gP=gridPars
  k<-dim(M1)[2]
  mag<-(mag-1)
  M2<-M2+(M2-M1)*mag
  limits = function(x,s){ 
    r = range(x)
    rc=scale(r,scale=FALSE)
    l=mean(r)+s*rc
  }
  if(k==2){
    if(method=="TPS"){
      tps(M1,M2,gP$n.col.cell, sz=gP$tar.pt.size, pt.bg=gP$tar.pt.bg, grid.col=gP$grid.col, 
          grid.lwd=gP$grid.lwd, grid.lty=gP$grid.lty, refpts=useRefPts)
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      if(label == TRUE){text(M2, label=paste(1:dim(M2)[1]),adj=gP$txt.adj,
                             pos=gP$txt.pos,cex=gP$txt.cex,col=gP$txt.col)}
      if(!is.null(outline)){
        curve.warp <- tps2d(outline, M1, M2)
        points(curve.warp,pch=19, cex=gP$tar.out.cex, col=gP$tar.out.col) 
      }
      if(useRefPts==FALSE) points(M2,pch=21,cex=gP$tar.pt.size, bg=gP$tar.pt.bg) else points(M1,pch=21,cex=gP$pt.size, bg=gP$pt.bg)
    }
    if(method=="vector"){
      plot(M1,asp=1,type="n",xlab="x",ylab="y",xlim=limits(M1[,1],1.25),
           ylim=limits(M1[,2],1.25),...)
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      if(label == TRUE){text(M1, label=paste(1:dim(M1)[1]),adj=gP$txt.adj,
                             pos=gP$txt.pos,cex=gP$txt.cex,col=gP$txt.col)}
      arrows(M1[,1],M1[,2],M2[,1],M2[,2],length=0.075,lwd=2)
      points(M1,pch=21,bg=gP$pt.bg,cex=gP$pt.size)
    }
    if(method=="points"){
      plot(M1,asp=1,pch=21,type="n",xlim=limits(M1[,1],1.25),
           ylim=limits(M1[,2],1.25),xlab="x",ylab="y",...)
      if(label == TRUE){text(M1, label=paste(1:dim(M1)[1]),adj=gP$txt.adj,
                             pos=gP$txt.pos,cex=gP$txt.cex,col=gP$txt.col)}
      if(!is.null(outline)){
        curve.warp <- tps2d(outline, M1, M2)
        points(outline,pch=19, cex=gP$out.cex, col=gP$out.col) 
        points(curve.warp,pch=19, cex=gP$tar.out.cex, col=gP$tar.out.col) 
      }
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty,nrow(links))[1:nrow(links)]
        tarlinkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        tarlinklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M1[links[i,1],1],M1[links[i,1],2],
                   M1[links[i,2],1],M1[links[i,2],2],col=linkcol[i],
                   lty=linklty[i],lwd=linklwd[i])
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],
                   M2[links[i,2],2],col=tarlinkcol[i],
                   lty=tarlinklty[i],lwd=tarlinklwd[i])
        }
      }
      points(M2,pch=21,bg=gP$tar.pt.bg,cex=gP$tar.pt.size)
      points(M1,pch=21,bg=gP$pt.bg,cex=gP$pt.size)
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
      tps(M1[,1:2],M2[,1:2],gP$n.col.cell, sz=gP$tar.pt.size, pt.bg=gP$tar.pt.bg, grid.col=gP$grid.col, 
          grid.lwd=gP$grid.lwd, grid.lty=gP$grid.lty, refpts=useRefPts)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          linkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
          linklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
          linklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
          segments(M2[links[i,1],1],M2[links[i,1],2],
                   M2[links[i,2],1],M2[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      title("X,Y tps grid")
      b<-c(1,3)
      tps(M1[,b],M2[,b],gP$n.col.cell, sz=gP$tar.pt.size, pt.bg=gP$tar.pt.bg, grid.col=gP$grid.col, 
          grid.lwd=gP$grid.lwd, grid.lty=gP$grid.lty, refpts=useRefPts)
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],3],
                   M2[links[i,2],1],M2[links[i,2],3],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      title("Y,Z tps grid")
      layout(1)
      on.exit(par(old.par))
    }
    if(method=="vector"){
      plot3d(M1,type="s",col=gP$pt.bg,size=gP$pt.size,aspect=FALSE,...)
      if(label == TRUE){text3d(M1, texts = paste(1:dim(M1)[1]), adj=(gP$txt.adj+gP$pt.size),
                               pos=(gP$txt.pos+gP$pt.size),cex=gP$txt.cex,col=gP$txt.col)}
      for (i in 1:nrow(M1)){
        segments3d(rbind(M1[i,],M2[i,]),lwd=2)
      }
      if(is.null(links)==FALSE){
        tarlinkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        tarlinklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),
                     col=tarlinkcol[i],lty=tarlinklty[i],lwd=tarlinklwd[i])
        }
      }
    }
    if(method=="points"){
      plot3d(M1,type="s",col=gP$pt.bg,size=gP$pt.size,aspect=FALSE,
             xlab="x", ylab="y",zlab="z",...)
      plot3d(M2,type="s", col=gP$tar.pt.bg,size=gP$tar.pt.size, add=TRUE)
      if(label == TRUE){text3d(M1, texts = paste(1:dim(M1)[1]), adj=(gP$txt.adj+gP$pt.size),
                               pos=(gP$txt.pos+gP$pt.size),cex=gP$txt.cex,col=gP$txt.col)}
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty,nrow(links))[1:nrow(links)]
        tarlinkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        tarlinklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments3d(rbind(M1[links[i,1],],M1[links[i,2],]),
                     col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),
                     col=tarlinkcol[i],lty=tarlinklty[i],lwd=tarlinklwd[i])
        }
      }
    }
    if (method == "surface") {
      if(is.null(mesh)==TRUE){
        stop("Surface plotting requires a template mesh3d object (see 'warpRefMesh').")
      }
      warp.PLY <- mesh
      vb <- as.matrix(t(mesh$vb)[,-4])
      warp <- tps2d3d(vb, M1, M2)
      warp.PLY$vb <- rbind(t(warp), 1)
      open3d(); shade3d(warp.PLY, ...)
      return(warp.PLY)
    }
  }
}
