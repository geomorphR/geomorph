#' Creates a 2D outline warped to the mean shape
#' 
#' A function to take an outline (defined by many points) and use thin-plate spline method to warp the outline
#' into the estimated mean shape for a set of aligned specimens.
#'
#' Function takes an outline (defined by many points) with a set of fixed landmark coordinates and uses the thin-plate spline method (Bookstein 1989) 
#' to warp the outline into the shape defined by a second set of landmark coordinates, usually those of the 
#' mean shape for a set of aligned specimens. It is highly recommended that the mean shape is used as the 
#' reference for warping (see Rohlf 1998). The workflow is as follows:
#' \enumerate{
#' \item {Calculate the mean shape using \code{\link{mshape}}}
#' \item{Choose an actual specimen to use for the warping. The specimen used as the template for this warping 
#' is recommended as one most similar in shape to the average of the sample, but can be any reasonable 
#' specimen - do this by eye, or use \code{\link{findMeanSpec}} }
#' \item{Warp this specimen into the mean shape using \code{\link{warpRefOutline}} }
#' \item{Use this average outline where it asks for a outline= in the analysis functions and visualization functions  }
#' }
#' The returned outline object is for use in geomorph
#' functions where shape deformations are plotted (\code{\link{plotTangentSpace}}, 
#' \code{\link{two.b.pls}}, \code{\link{bilat.symmetry}}, and \code{\link{plotRefToTarget}}). 
#' 
#' @param file A .txt or .csv file of the outline point coordinates, or a .TPS file with OUTLINES= or CURVES= elements
#' @param coord A p x k matrix of 2D fixed landmark coordinates
#' @param ref A p x k matrix of 2D coordinates made by \code{\link{mshape}}
#' @export
#' @seealso \code{\link{findMeanSpec}}
#' @keywords utilities
#' @keywords visualization
#' @author Emma Sherratt
#' @return Function returns an outline object
#' @references  Bookstein, F. L. 1989 Principal Warps: Thin-Plate Splines and the Decomposition
#' of Deformations. IEEE Transactions on Pattern Analysis and Machine Intelligence 11(6):567-585.
#' @references  Rohlf, F. J. 1998. On Applications of Geometric Morphometrics to Studies of Ontogeny and Phylogeny. Systematic Biology. 47:147-158.
warpRefOutline <- function(file, coord, ref){
  checkmat <- is.matrix(coord)
  if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(coord)[2]
  if (checkdim==3) {stop("Input must be a p-x-k matrix of two-dimensional landmark coordinates") }
  if(grepl(".txt", file, ignore.case=TRUE) == TRUE) {outline <- as.matrix(read.table(file, header = F))[,1:2]
                                                     npoints <- as.vector(nrow(outline)) }
  if(grepl(".csv", file, ignore.case=TRUE) == TRUE) {outline <- as.matrix(read.csv(file, header = F))[,1:2]
                                                     npoints <- as.vector(nrow(outline))}
  if(grepl(".tps", file, ignore.case=TRUE) == TRUE) {
    tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
    curves <- grep("POINTS=", tpsfile)
    npoints <- as.numeric(sub("POINTS=", "", tpsfile[curves]))
    outline = NULL
    for (i in 1:length(curves)){
      tmp <- tpsfile[(curves[i]+1): ((curves[i])+npoints[i])]
      outline <- rbind(outline, matrix(as.numeric(unlist(strsplit(tmp, split = " +"))), ncol = 2, byrow = T)) } 
    imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE", tpsfile)]))
    if (identical(imscale, numeric(0)) == TRUE) {imscale = 1}
    outline <- outline * imscale
  }
  if(sum(range(outline[,1])) < sum(range(outline[,2]))){ layout(t(c(1,2))) }
  if(sum(range(outline[,1])) > sum(range(outline[,2]))){ layout(c(1,2)) }
  plot(outline, pch=19, cex=0.3, main = "Imported outline", asp=T, xlab="x", ylab="y")
  points(coord, pch=19, col="red")
  text(coord, labels = c(1:nrow(coord)), adj=2)
  coord.sc <- scale(coord, scale=F)
  sc.mat <- matrix(rep(1,nrow(outline)), ncol=1) %*% apply(coord,2,mean)
  outline <- outline - sc.mat
  warp <- tps2d(outline, coord.sc, ref)
    plot(warp, pch=19, cex=0.3, main = "Warped outline", asp=T, xlab="x", ylab="y")
    points(ref, pch=19, cex=0.8, col= "red")
    text(ref, labels = c(1:nrow(ref)), adj=2)
  layout(1)
return(list(outline=warp, npoints = npoints))
}
