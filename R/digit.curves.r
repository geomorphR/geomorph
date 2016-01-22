#' Calculate semilandmarks along a curve
#'
#' A function to calculate equidistant two-dimensional and three-dimensional semilandmarks along a curve. These landmarks will be treated as "sliders" in Generalized Procrustes analysis \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks (see Bookstein 1997 for algorithm details). Each sliding semilandmark ("sliders") will slide between two designated points, along a line tangent to the specified curvature, as specified by \code{\link{define.sliders}}. 
#'
#' The function is based upon tpsDig2 'resample curve by length' for 2D data by James Rohlf.
#' The start of the curve is a fixed landmark on the curve that is equivalent (homologous) in each specimen in the sample (and will be treated as a fixed point during Procrustes Superimpoistion using \code{\link{gpagen}}). Then nPoints are calculated along the curve at equidistant points from the start to the end. 
#' 
#' 'curve' is a p-x-k matrix of 2D or 3D coordinates for a set of ordered points defining a curve. This can be the pixels of an outline calculated in ImageJ (save xy coordinates) or any other reasonable way of obtaining ordered coordinates along a curve (including sampling by hand using 
#' \code{\link{digit.fixed}} or \code{\link{digitize2d}} - but note that there should be more points defining the curve than nPoints in order to accurately calculate the semilandmarks).
#'
#' If 'closed = T', the function returns the coordinates of the 'start' landmark plus nPoints. If 'closed = F', the function returns the coordinates of the 'start' landmark, plus nPoints and the end of the curve. 
#'
#' @param start A vector of x,y,(z) coordinates for the fixed landmark defining the start of the curve
#' @param curve A p-x-k matrix of 2D or 3D coordinates for a set of ordered points defining a curve
#' @param nPoints Numeric how many semilandmarks to place equidistantly along the curve 
#' @param closed Logical Whether the curve is closed (TRUE) or open (FALSE)
#' @return Function returns a matrix of coordinates for nPoints equally spaced semilandmarks sampled along the curve (plus start and end if 'closed = F', or only including start if 'closed = T') 
#' @seealso \code{\link{digit.fixed}} \code{\link{digitize2d}}
#' @export
#' @keywords digitizing
#' @author Emma Sherratt
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.

digit.curves <- function(start, curve, nPoints, closed=T){
  nPoints=nPoints+2
  checkmat <- is.matrix(curve)
  if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(curve)[2]
  nCurvePoints = nrow(curve)
  if (checkdim==2) {  newPoints <- matrix(NA, ncol=2, nrow = nPoints)
                      start <- as.numeric(start)
                      start <- which.min(sqrt((start[1]-curve[,1])^2+
                                              (start[2]-curve[,2])^2))
                    }
  
  if (checkdim==3) {  newPoints <- matrix(NA, ncol=3, nrow = nPoints) 
                      start <- as.numeric(start)
                      start <- which.min(sqrt((start[1]-curve[,1])^2+
                                              (start[2]-curve[,2])^2+
                                              (start[3]-curve[,3])^2))
                      }
  newPoints[1,] <- curve[start,]
  if(start!=1){curve <- rbind(curve[start:nCurvePoints,],
                               curve[1:(start-1),])} 
  if(closed==F){newPoints[nPoints,] <- curve[nrow(curve),]}
  if(closed==T){curve <- rbind(curve, curve[1,])
                nCurvePoints <- nCurvePoints+1
                newPoints[nPoints,] <- curve[nCurvePoints,]}
  B <- rep(0, nCurvePoints) 
  for(i in 1:(nCurvePoints-1)){ 
    if (checkdim==2) {Interval<-sqrt((curve[i,1]-curve[i+1,1])^2 
                                       + (curve[i,2]-curve[i+1,2])^2)} 
    if (checkdim==3) {Interval<-sqrt((curve[i,1]-curve[i+1,1])^2 
                                       + (curve[i,2]-curve[i+1,2])^2
                                       + (curve[i,3]-curve[i+1,3])^2)}
  B[i+1]<-B[i]+Interval} 
  TotalLength <- B[nCurvePoints]
  j = 2
  for (i in 2:(nPoints-1)){
    NextLength <- TotalLength*(i - 1) / (nPoints - 1) 
    while(B[j-1] < NextLength) {j=j+1} 
    xy0 <- curve[j - 2,] 
    xy <- curve[j - 1,] 
    CurrInterval <- B[j - 1] - B[j - 2] 
    if (CurrInterval > 0){p <- (NextLength - B[j - 2]) / CurrInterval } else p <- 0
    newPoints[i,1] <- round((1 - p) * xy0[1] + p * xy[1], digits=4) 
    newPoints[i,2] <- round((1 - p) * xy0[2] + p * xy[2], digits=4) 
    if (checkdim==3) {newPoints[i,3] <- round((1 - p) * xy0[3] + p * xy[3], digits=4)}
  }
  if (closed==T){return(newPoints[1:(nPoints-1),])}
  if (closed==F){return(newPoints[1:nPoints,])}
}