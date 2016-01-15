#' Plot landmark coordinates for all specimens
#'
#' Function plots landmark coordinates for a set of specimens
#'
#' The function creates a plot of the landmark coordinates for all specimens. This is useful for examining 
#'  patterns of shape variation after GPA. If "mean=TRUE", the mean shape will be calculated and added to the plot.
#'  Additionally, if a matrix of links is provided, the landmarks of the mean shape will be connected by lines.  
#'  The link matrix is an m x 2 matrix, where m is the desired number of links. Each row of the link matrix 
#'  designates the two landmarks to be connected by that link. The function will plot either two- or 
#'  three-dimensional data.
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param mean A logical value indicating whether the mean shape should be included in the plot
#' @param links An optional matrix defining for links between landmarks (only if mean=TRUE)
#' @param label A logical value indicating whether landmark numbers will be plotted (only if mean=TRUE)
#' @param pt.bg An optional value defining the background color of points (single value or vector of values)
#' @param pt.cex An optional value defining the the size of the points (single value or vector of values)
#' @param mean.bg An optional value defining the background color of the points for all specimens
#' @param mean.cex An optional value defining the size of the points representing the average specimen
#' @param link.col An optional value defining color of links (single value or vector of values)
#' @param link.lwd An optional value defining line weight of links (single value or vector of values)
#' @param link.lty An optional value defining line type of links (single value or vector of values)
#' @param txt.adj The adjustment value of the landmark label (one or two values, as in base R \code{\link{text}}) 
#' @param txt.pos The position of the landmark label (single numerical value, as in base R \code{\link{text}}) 
#' @param txt.cex The size of the landmark label text (single numerical value, as in base R \code{\link{text}})
#' @param txt.col The color of the landmark label text (single numerical value, as in base R \code{\link{text}})
#' @export
#' @keywords visualization
#' @author Dean Adams
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'
#' plotAllSpecimens(Y.gpa$coords,links=plethodon$links)
plotAllSpecimens<-function(A,mean=TRUE,links=NULL,label=FALSE,
                           pt.bg="gray",
                           pt.cex=1,
                           mean.bg="black",
                           mean.cex=2,
                           link.col="black",
                           link.lwd = 2,
                           link.lty = 1,
                           txt.adj = 0.5,
                           txt.pos = 1, 
                           txt.cex = 0.8,
                           txt.col = "black"){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  k<-dim(A)[2]
  if(mean==TRUE){
    mn<-mshape(A)
  }
  if(k==2){
    plot(A[,1,],A[,2,],asp=1, pch=21,bg=pt.bg,cex=pt.cex*1,xlab="x",ylab="y") 
    if(mean==TRUE){ 
      if(is.null(links)==FALSE){
        linkcol <- rep(link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(mn[links[i,1],1],mn[links[i,1],2],mn[links[i,2],1],mn[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      points(mn,pch=21,bg=mean.bg,cex=mean.cex)
      if(label == TRUE){text(mn, label=paste(1:dim(mn)[1]),adj=(txt.adj+mean.cex),
                             pos=txt.pos,cex=txt.cex,col=txt.col)}
    }
  }
  if(k==3){
    A3d<-NULL
    for (i in 1:dim(A)[[3]]){
      A3d<-rbind(A3d,A[,,i])
    }
    plot3d(A3d,type="s",col=pt.bg,xlab="x",ylab="y",zlab="z",size=pt.cex*1.5,aspect=FALSE)
    if(mean==TRUE){ 
      if(is.null(links)==FALSE){
        linkcol <- rep(link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments3d(rbind(mn[links[i,1],],mn[links[i,2],]),
                     col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      points3d(mn,color=mean.bg,size=mean.cex*2)
      if(label == TRUE){text3d(mn, texts = paste(1:dim(mn)[1]), adj=(txt.adj+mean.cex),
                               pos=txt.pos,cex=txt.cex,col=txt.col)}
    }
  }
}