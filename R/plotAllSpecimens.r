#' Plot landmark coordinates for all specimens
#'
#' Function plots landmark coordinates for a set of specimens
#'
#' The function creates a plot of the landmark coordinates for all specimens. This is useful for examining 
#'  patterns of shape variation after GPA. If "mean=TRUE", the mean shape will be calculated and added to the plot.
#'  Additionally, if a matrix of links is provided, the landmarks of the mean shape will be connected by lines.  
#'  The link matrix is an m x 2 matrix, where m is the desired number of links. Each row of the link matrix 
#'  designates the two landmarks to be connected by that link. The function will plot either two- or 
#'  three-dimensional data (e.g. see \code{\link{define.links}}).
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param mean A logical value indicating whether the mean shape should be included in the plot
#' @param links An optional matrix defining for links between landmarks (only if mean=TRUE)
#' @param label A logical value indicating whether landmark numbers will be plotted (only if mean=TRUE)
#' @param plot.param A list of plotting parameters for the points (pt.bg, pt.cex), mean (mean.bg, mean.cex), links (link.col, link.lwd, link.lty) and landmark labels (txt.cex, txt.adj, txt.pos, txt.col)
#' @export
#' @keywords visualization
#' @author Dean Adams
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'
#' plotAllSpecimens(Y.gpa$coords,links=plethodon$links)
plotAllSpecimens<-function(A,mean=TRUE,links=NULL,label=FALSE,plot.param = list()){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  k<-dim(A)[2]
  if(mean==TRUE){
    mn<-mshape(A)
  }
  
  p.p <- plot.param
  if(is.null(p.p$pt.bg)) p.p$pt.bg="gray" ; if(is.null(p.p$pt.cex)) p.p$pt.cex=1 ; 
  if(is.null(p.p$mean.bg)) p.p$mean.bg="black" ; if(is.null(p.p$mean.cex)) p.p$mean.cex=2
  if(is.null(p.p$link.col)) p.p$l.col="black" ; if(is.null(p.p$link.lwd)) p.p$link.lwd=2
  if(is.null(p.p$link.lty)) p.p$link.lty=1 ; if(is.null(p.p$txt.adj)) p.p$txt.adj=c(-.1,-.1)
  if(is.null(p.p$txt.col)) p.p$txt.col="black" ; if(is.null(p.p$txt.cex)) p.p$txt.cex=0.8
  if(is.null(p.p$txt.pos)) p.p$txt.pos=1
  
  if(k==2){
    plot(A[,1,],A[,2,],asp=1, pch=21,bg=p.p$pt.bg,cex=p.p$pt.cex*1,xlab="x",ylab="y") 
    if(mean==TRUE){ 
      if(is.null(links)==FALSE){
        linkcol <- rep(p.p$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(p.p$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(p.p$link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(mn[links[i,1],1],mn[links[i,1],2],mn[links[i,2],1],mn[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      points(mn,pch=21,bg=p.p$mean.bg,cex=p.p$mean.cex)
      if(label == TRUE){text(mn, label=paste(1:dim(mn)[1]),adj=(p.p$txt.adj+p.p$mean.cex),
                             pos=p.p$txt.pos,cex=p.p$txt.cex,col=p.p$txt.col)}
    }
  }
  if(k==3){
    A3d<-NULL
    for (i in 1:dim(A)[[3]]){
      A3d<-rbind(A3d,A[,,i])
    }
    plot3d(A3d,type="s",col=p.p$pt.bg,xlab="x",ylab="y",zlab="z",size=p.p$pt.cex*1.5,aspect=FALSE)
    if(mean==TRUE){ 
      if(is.null(links)==FALSE){
        linkcol <- rep(p.p$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(p.p$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(p.p$link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments3d(rbind(mn[links[i,1],],mn[links[i,2],]),
                     col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      points3d(mn,color=p.p$mean.bg,size=p.p$mean.cex*2)
      if(label == TRUE){text3d(mn, texts = paste(1:dim(mn)[1]), adj=(p.p$txt.adj+p.p$mean.cex),
                               pos=p.p$txt.pos,cex=p.p$txt.cex,col=p.p$txt.col)}
    }
  }
}