#' Plot specimens in tangent space
#'
#' Function plots a set of Procrustes shape variables in tangent space along their principal component axes
#'
#' The function performs a principal components analysis of shape variation and plots two 
#'  dimensions of tangent space for a set of Procrustes shape variables (default is PC1 vs. PC2). 
#'  The percent variation along each PC-axis is returned. Additionally (and optionally, {warpgrids=T}), 
#'  deformation grids can be requested, which display the shape of specimens at the ends 
#'  of the range of variability along PC1. If groups are provided, specimens from 
#'  each group are plotted using distinct colors based on the order in which the groups are found in the dataset, 
#'  and using R's standard color palette: black, red, green, blue, cyan, magenta, yellow, and gray. NOTE: to change
#'  the colors of the groups, simply substitute a vector of the desired colors for each specimen (see example below).
#'  
#'  \subsection{Notes for geomorph 3.1.0 and subsequent versions}{ 
#'  The function \code{\link{gm.prcomp}} can also be used to generate principal components plots of 
#'  tangent space, and will yield plots identical to those of the current function. 
#'  }
#'  
#'  NOTE: previous versions of plotTangentSpace had option 'verbose' to return the PC scores and PC shapes. 
#'  From version 3.0.2 this is automatic when assigned to an object.
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of specimens 
#' @param warpgrids A logical value indicating whether deformation grids for shapes along X-axis should be displayed
#' @param mesh A mesh3d object to be warped to represent shape deformation along X-axis (when {warpgrids=TRUE})
#' as described in \code{\link{plotRefToTarget}}.
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param label An optional vector indicating labels for each specimen are to be displayed 
#' (or if TRUE, numerical addresses are given)
#' @param groups An optional factor vector specifying group identity for each specimen (see example)
#' @param legend A logical value for whether to add a legend to the plot (only when groups are assigned).
#' @param ... Arguments passed on to \code{\link{prcomp}}.  By default, \code{\link{plotTangentSpace}}
#' will attempt to remove redundant axes (eigen values effectively 0).  To override this, adjust the 
#' argument, tol, from \code{\link{prcomp}}.
#' @return If user assigns function to object, returned is a list of the following components:
#' \item{pc.summary}{A table summarizing the percent variation explained by each pc axis, equivalent to summary of \code{\link{prcomp}}.}
#' \item{pc.scores}{The set of principal component scores for all specimens.}
#' \item{pc.shapes}{A list with the shape coordinates of the extreme ends of all PC axes, e.g. $PC1min}
#' \item{sdev}{The standard deviations of the principal components (i.e., the square roots of the eigenvalues of the 
#' covariance/correlation matrix, as per \code{\link{prcomp}}.}
#' \item{rotation}{The matrix of variable loadings, as per \code{\link{prcomp}}.}
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @seealso \code{\link{gm.prcomp}}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' 
#' gp <- interaction(plethodon$species, plethodon$site) # group must be a factor
#' plotTangentSpace(Y.gpa$coords, groups = gp) 
#' 
#' ## To save and use output
#' PCA <- plotTangentSpace(Y.gpa$coords, groups = gp, legend=TRUE) 
#' summary(PCA)
#' PCA$pc.shapes
#' PCA$rotation
#' 
#' ##To change colors of groups
#' col.gp <- rainbow(length(levels(gp))) 
#'    names(col.gp) <- levels(gp)
#' col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must NOT be a factor
#' plotTangentSpace(Y.gpa$coords, groups = col.gp)
#' 
#' ## To plot residual shapes from an allometry regression (note: must add mean back in!) 
#' plotTangentSpace(arrayspecs(resid(lm(two.d.array(Y.gpa$coords)~log(Y.gpa$Csize)))+
#'          predict(lm(two.d.array(Y.gpa$coords)~1)),12,2))

plotTangentSpace<-function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, label = NULL, 
                            groups=NULL, legend=FALSE, ...){
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }
  dots <- list(...)
  retx <- dots$retx
  if(is.null(retx)) retx <- TRUE
  scale. <- dots$scale.
  if(is.null(scale.)) scale. <- FALSE
  center <- dots$center
  if(is.null(center)) center <- TRUE
  tol <- dots$tol
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  if(is.null(tol)){
    d <- prcomp(x)$sdev^2
    cd <-cumsum(d)/sum(d)
    cd <- length(which(cd < 1)) 
    if(length(cd) < length(d)) cd <- cd + 1
    if(length(d) > 2) tol <- max(c(d[cd]/d[1],0.005)) else tol <- 0
  }
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, tol = tol)
  pcdata <- pc.res$x
  if (warpgrids == FALSE) {
    if(legend==TRUE){ layout(t(matrix(c(1, 1, 2, 1, 1, 1, 1, 1, 1), 3,3))) }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
         ylab = paste("PC ",axis2))
    if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    if (length(label!=0)) {
      if(isTRUE(label)){text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) }
      else{text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.7)) }
    }
    if(!is.null(groups) && legend==TRUE){
        plot.new(); 
        if(is.factor(groups)){legend(0.5,1, legend=unique(groups), pch=19, bty="n", col=unique(groups))
        } else {legend(0.5,1, legend=unique(names(groups)), pch=19, bty="n", col=unique(groups)) }
    }
  }
  shapes <- shape.names <- NULL
  for(i in 1:ncol(pcdata)){
    pcaxis.min <- min(pcdata[, i]) ; pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min ; pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + as.vector(t(ref))
    shapes <- rbind(shapes,pc.min, pc.max)
    shape.names <- c(shape.names,paste("PC",i,"min", sep=""),paste("PC",i,"max", sep=""))
  }
  shapes <- arrayspecs(shapes,p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[,,x])
  names(shapes) <- shape.names
  if (warpgrids == TRUE) {
    if (k == 2) {
      layout(t(matrix(c(2, 1, 4, 1, 1, 1, 1, 1, 3), 3,3)))
    }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
         ylab = paste("PC ",axis2))
      if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    if (length(label!=0)) {
      if(isTRUE(label)){text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) }
      else{text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.1)) }
    }
    shape.min <- shapes[[which(names(shapes) == paste("PC",axis1,"min", sep=""))]]
    shape.max <- shapes[[which(names(shapes) == paste("PC",axis1,"max", sep=""))]]
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[,axis2])), min(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[,axis2])), max(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      tps(ref, shape.min, 20)
      tps(ref, shape.max, 20)
    }
    if(!is.null(groups) && legend==TRUE){
      plot.new(); 
      if(is.factor(groups)){legend(0.5,1, legend=unique(groups), pch=19, bty="n", col=unique(groups))
        } else {legend(0.5,1, legend=unique(names(groups)), pch=19, bty="n", col=unique(groups)) }
      }
    if (k == 3) {
      if (is.null(mesh)==TRUE){
        open3d() ; mfrow3d(1, 2) 
        plot3d(shape.min, type = "s", col = "gray", main = paste("PC ", axis1," negative"),size = 1.25, aspect = FALSE,xlab="",ylab="",zlab="",box=FALSE, axes=FALSE)
        plot3d(shape.max, type = "s", col = "gray", main = paste("PC ", axis1," positive"), size = 1.25, aspect = FALSE,xlab="",ylab="",zlab="",box=FALSE, axes=FALSE)
        }
      if(is.null(mesh)==FALSE){
        open3d() ; mfrow3d(1, 2) 
        cat(paste("\nWarping mesh to negative end of axis ", axis1, "\n", sep=""))
        plotRefToTarget(ref, shape.min, mesh, method = "surface")
        title3d(main=paste("PC ", axis1," negative"))
        next3d()
        cat(paste("\nWarping mesh to positive end of axis ", axis1, "\n", sep=""))
        plotRefToTarget(ref, shape.max, mesh, method = "surface")
        title3d(main=paste("PC ", axis1," positive"))
        }
      }
    layout(1)
    }
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata, pc.shapes = shapes, 
              sdev = pc.res$sdev, rotation = pc.res$rotation)
  class(out) = "plotTangentSpace"
  out
}
