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
#' \subsection{Notes for geomorph 3.3.0 and subsequent versions}{
#'  This function is defunct and replaced with \code{\link{gm.prcomp}}.  The functions,
#'  \code{\link{summary.gm.prcomp}}, \code{\link{plot.gm.prcomp}}, and \code{\link{picknplot.shape}} 
#'  allow summaries, ammendable plots, and warpgrids, for points anywhere in a plot, respectively, can
#'  be used to accomplish any task of plotTangentSpace, with greater flexibility.
#'  }
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

plotTangentSpace<-function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, label = NULL, 
                            groups=NULL, legend=FALSE, ...){
  .Defunct("gm.prcomp", package = "geomorph", 
              msg = "plotTangentSpace has been removed from geomorph.
              A combination of gm.prcomp, summary.prcomp, and plot.gm.prcomp has much greater flexibility.
              One can also use picknplot.shape to generate warpgrids, anywhere in a plot.")
}
