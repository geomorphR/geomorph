#' Plot specimens in tangent space
#'
#' Function plots a set of Procrustes-aligned specimens in tangent space along their principal axes
#'
#' The function performs a principal components analysis of shape variation and plots two 
#'  dimensions of tangent space for a set of Procrustes-aligned specimens (default is PC1 vs. PC2). 
#'  The percent variation along each PC-axis is returned. Additionally (and optionally, {warpgrids=T}), 
#'  deformation grids can be requested, which display the shape of specimens at the ends 
#'  of the range of variability along PC1. If groups are provided, specimens from 
#'  each group are plotted using distinct colors based on the order in which the groups are found in the dataset, 
#'  and using R's standard color palette: black, red, green, blue, cyan, magenta, yellow, and gray. NOTE: to change
#'  the colors of the groups, simply substitute a vector of the desired colors for each specimen (see example below).
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of aligned specimens 
#' @param warpgrids A logical value indicating whether deformation grids for shapes along X-axis should be displayed
#' @param mesh A mesh3d object to be warped to represent shape deformation along X-axis (when {warpgrids=TRUE})
#' as described in \code{\link{plotRefToTarget}}.
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param label An optional vector indicating labels for each specimen are to be displayed 
#' (or if TRUE, numerical addresses are given)
#' @param groups An optional factor vector specifying group identity for each specimen
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @return Function returns a table summarizing the percent variation explained by each
#'  pc axis (equivalent to summary() of \code{\link{prcomp}}) (default, verbose = FALSE).
#'  If verbose=TRUE, function returns a list containing the following components: 
#' \item{pc.summary}{A PC summary table as above}
#' \item{pc.scores}{The set of principal component scores for all specimens.} 
#' \item{pc.shapes}{A list with four components of the shape coordinates of the extreme ends of axis1 and axis2 
#' is returned, which can be used by \code{\link{warpRefMesh}}.}
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' 
#' gp <- as.factor(paste(plethodon$species, plethodon$site)) # group must be a factor
#' plotTangentSpace(Y.gpa$coords, groups = gp) 
#' 
#' ## To plot residual shapes from an allometry regression (note: must add mean back in!) 
#' plotTangentSpace(arrayspecs(resid(lm(two.d.array(Y.gpa$coords)~Y.gpa$Csize))+
#'          predict(lm(two.d.array(Y.gpa$coords)~1)),12,2))
#'  
#' ##To change colors of groups
#' col.gp <- rainbow(length(levels(gp))) 
#'    names(col.gp) <- levels(gp)
#' col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must not be a factor
#' 
#' plotTangentSpace(Y.gpa$coords, groups = col.gp)
#' 
#' ## To make a quick legend [not run]
#' #plot.new(); legend(0,1, legend= levels(gp), pch=19, col = levels(as.factor(col.gp)))
plotTangentSpace<-function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, label = NULL, 
                            groups=NULL, verbose=FALSE){
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  pc.res <- prcomp(x)
  pcdata <- pc.res$x
  if (warpgrids == FALSE) {
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
         ylab = paste("PC ",axis2))
    if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    if (length(label!=0)) {
      if(isTRUE(label)){text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) }
      else{text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.7)) }
    }
  }
  pcaxis.min.1 <- min(pcdata[, axis1])
  pcaxis.max.1 <- max(pcdata[, axis1])
  pc.min.1 <- pc.max.1 <- rep(0, dim(pcdata)[2])
  pc.min.1[axis1] <- pcaxis.min.1
  pc.max.1[axis1] <- pcaxis.max.1
  shape.min.1 <- arrayspecs(as.matrix(pc.min.1 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shape.max.1 <- arrayspecs(as.matrix(pc.max.1 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  pcaxis.min.2 <- min(pcdata[, axis2])
  pcaxis.max.2 <- max(pcdata[, axis2])
  pc.min.2 <- pc.max.2 <- rep(0, dim(pcdata)[2])
  pc.min.2[axis2] <- pcaxis.min.2
  pc.max.2[axis2] <- pcaxis.max.2
  shape.min.2 <- arrayspecs(as.matrix(pc.min.2 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shape.max.2 <- arrayspecs(as.matrix(pc.max.2 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shapes <- list(shape.min.1, shape.max.1, shape.min.2, shape.max.2)
  names(shapes) <- c(paste("PC",axis1,"min", sep=""),paste("PC",axis1,"max", sep=""),
                     paste("PC",axis2,"min", sep=""),paste("PC",axis2,"max", sep=""))
  if (warpgrids == TRUE) {
    if (k == 2) {
      layout(t(matrix(c(2, 1, 1, 1, 1, 1, 1, 1, 3), 3,3)))
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
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[,axis2])), min(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[,axis2])), max(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      tps(ref, shape.min.1, 20)
      tps(ref, shape.max.1, 20)
    }
    if (k == 3) {
      if (is.null(mesh)==TRUE){
        open3d()
        plot3d(shape.min.1, type = "s", col = "gray", main = paste("PC ", axis1," negative"),size = 1.25, aspect = FALSE)
        open3d()
        plot3d(shape.max.1, type = "s", col = "gray", main = paste("PC ", axis1," positive"), size = 1.25, aspect = FALSE)
        }
      if(is.null(mesh)==FALSE){
        print("Warping mesh to axis 1 minima and maxima...")
        plotRefToTarget(ref, shape.min.1, mesh, method = "surface")
        title3d(main=paste("PC ", axis1," negative"))
        plotRefToTarget(ref, shape.max.1, mesh, method = "surface")
        title3d(main=paste("PC ", axis1," positive"))
        }
      }
    layout(1)
    }
  if(verbose==TRUE){ return(list(pc.summary = summary(pc.res), pc.scores = pcdata, pc.shapes= shapes)) }
  if(verbose==FALSE){ return(pc.summary = summary(pc.res)) }
}