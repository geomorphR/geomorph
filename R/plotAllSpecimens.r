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
#' @param links An optional matrix defining for links between landmarks
#' @param pointscale An optional value defining the size of the points for all specimens
#' @param meansize An optional value defining the size of the points representing the average specimen
#' @export
#' @keywords visualization
#' @author Dean Adams
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'
#' plotAllSpecimens(Y.gpa$coords,links=plethodon$links)
plotAllSpecimens<-function(A,mean=TRUE,links=NULL,pointscale=1,meansize=2){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  k<-dim(A)[2]
  if(mean==TRUE){
    mn<-mshape(A)
  }
  if(k==2){
    plot(A[,1,],A[,2,],asp=1, pch=21,bg="gray",xlab="x",ylab="y",cex=pointscale*1) 
    if(mean==TRUE){ 
      points(mn,pch=21,bg="black",cex=meansize) 
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(mn[links[i,1],1],mn[links[i,1],2],mn[links[i,2],1],mn[links[i,2],2],lwd=2)
        }
      }
    }
  }
  if(k==3){
    A3d<-NULL
    for (i in 1:dim(A)[[3]]){
      A3d<-rbind(A3d,A[,,i])
    }
    plot3d(A3d,type="s",col="gray",xlab="x",ylab="y",zlab="z",size=pointscale*1.5,aspect=FALSE)
    if(mean==TRUE){ 
      points3d(mn,color="black",size=meansize*2)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(mn[links[i,1],],mn[links[i,2],]),lwd=3)
        }
      }
    }
  }
}