#' Calculate semilandmarks along a curve
#'
#' A function that "digitizes curves" by calculating equidistant two-dimensional or three-dimensional semilandmarks along a curve. These landmarks will be treated as "sliders" in Generalized Procrustes analysis \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks (see Bookstein 1997 for algorithm details). Each sliding semilandmark ("sliders") will slide between two designated points, along a line tangent to the specified curvature, as specified by \code{\link{define.sliders}}. 
#'
#' The function is based upon tpsDig2 'resample curve by length' for 2D data by James Rohlf (Rohlf 2015).
#' The start of the curve is a fixed landmark on the curve that is equivalent (homologous) in each specimen in the sample (and will be treated as a fixed point during Procrustes Superimposition using \code{\link{gpagen}}). Then nPoints are calculated along the curve at equidistant points from the start to the end. 
#' 
#' 'curve' is a p-x-k matrix of 2D or 3D coordinates for a set of ordered points defining a curve. This can be the pixels of an outline calculated in ImageJ (save xy coordinates) or any other reasonable way of obtaining ordered coordinates along a curve (including sampling by hand using 
#' \code{\link{digit.fixed}} or \code{\link{digitize2d}} - but note that there should be more points defining the curve than nPoints in order to accurately calculate the semilandmarks).
#'
#' If 'closed = T', the function returns the coordinates of the 'start' landmark plus nPoints. If 'closed = F', the function returns the coordinates of the 'start' landmark, plus nPoints and the end of the curve. 
#' 
#' If unsure if the points defining the curve are ordered, then plot and color them using the rainbow function, e.g. plot(curve, pch=19, cex=0.1, col=rainbow(nrow(outline))), and it should be easy to visualize.
#'
#' @param start A numeric vector of x,y,(z) coordinates for the landmark defining the start of the curve (can be simply first point on open outline: curve[1,])
#' @param curve A matrix (p x k) of 2D or 3D coordinates for a set of ordered points defining a curve
#' @param nPoints Numeric how many semilandmarks to place equidistantly along the curve (not counting beginning and end points) 
#' @param closed Logical Whether the curve is closed (TRUE) or open (FALSE)
#' @return Function returns a matrix of coordinates for nPoints equally spaced semilandmarks sampled along the curve (plus start and end if 'closed = F', or only including start if 'closed = T') 
#' @seealso \code{\link{digit.fixed}} \code{\link{digitize2d}}
#' @export
#' @keywords digitizing
#' @author Emma Sherratt and Michael Collyer
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.
#' @references Rohlf, F.J., 2015. The tps series of software. Hystrix 26(1):9-12.

digit.curves <- function(start, curve, nPoints, closed=TRUE){
  nPoints <- nPoints+2
  if(!is.matrix(curve)) stop("Input must be a p-x-k matrix of curve coordinates")
  nCurvePoints = NROW(curve)
  if(nCurvePoints < 2) stop("curve matrix does not have enough points to estimate any interior points")
  if(nPoints > (nCurvePoints - 1)) {
    if((nCurvePoints - 1) == 1) nPoints = 1
    if((nCurvePoints - 1) > 1) nPoints = nCurvePoints - 2
    cat("\nWarning: because the number of desired points exceeds the number of curve points,")
    cat("\nthe number of points will be truncated to", nPoints, "\n\n")
  }
  start <- as.numeric(start)
  if(!setequal(start, curve[1,])) curve <- rbind(start, curve)
  if(closed) curve <- rbind(curve, curve[1,])
  res <- evenPts(curve, nPoints)
  if(closed) res <- res[-NROW(res),]
  res

}