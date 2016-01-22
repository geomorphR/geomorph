#' Plot 3D specimen, fixed landmarks and surface semilandmarks
#'
#' A function to plot three-dimensional (3D) specimen along with its landmarks.
#' 
#' Function to plot 3D specimens along with their digitized "fixed" landmarks and semilandmarks
#' "surface sliders" and "curve sliders". If specimen is a 3D surface (class shape3d/mesh3d) mesh is plotted.
#' For visualization purposes, 3D coordinate data collected using \code{\link{digit.fixed}} or 
#' \code{\link{digitsurface}} and \code{\link{buildtemplate}} prior to build 1.1-6 were centered by default. 
#' Therefore use this function with {centered=TRUE}. Data collected outside geomorph should be read using
#' {centered=FALSE}. The function assumes the fixed landmarks are listed at the beginning of
#' the coordinate matrix ({digitspec}).
#'
#' @param spec An object of class shape3d/mesh3d, or matrix of 3D vertex coordinates.
#' @param digitspec Name of data matrix containing 3D fixed and/or surface sliding coordinates.
#' @param fixed Numeric The number of fixed template landmarks (listed first in {digitspec})
#' @param ptsize Numeric Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes                                       
#' @param centered Logical Whether the data matrix is in the surface mesh coordinate system ({centered=FALSE}) or
#' if the data were collected after the mesh was centered ({centered=TRUE})- see details.
#' @param ... additional parameters which will be passed to \code{\link{plot3d}}.
#' @export
#' @keywords visualization
#' @seealso \code{\link{warpRefMesh}}
#' @seealso \code{\link{read.ply}}
#' @examples
#' # data(scallopPLY)
#' # ply <- scallopPLY$ply
#' # digitdat <- scallopPLY$coords
#' # plotspec(spec=ply,digitspec=digitdat,fixed=16, centered =TRUE)
#' @author Erik Otarola-Castillo & Emma Sherratt
plotspec <- function (spec, digitspec, fixed=NULL, ptsize = 1, centered = FALSE, ...) 
{
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    if(centered == TRUE){
      specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    if(centered == FALSE){specimen <- as.matrix(t(spec$vb)[,-4])}
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" }
  } 
  else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } 
  else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    if(centered == TRUE){ specimen <- scale(spec, scale = FALSE)}
    if(centered == FALSE){specimen <- spec}
  } 
  else { stop ("File is not matrix in form: vertices by xyz")} 
  if (is.null(dim(digitspec)) || dim(digitspec)[2] != 3) {
    stop("Digitized file is not xyz matrix in form: p x k")}
  plot3d(specimen[, 1], specimen[, 2], specimen[, 3], size = ptsize, aspect = FALSE, ...)
    if (!is.null(mesh)) { shade3d(mesh, add=TRUE) }
  if(!is.null(fixed)) {points3d(digitspec[1:fixed, ], aspect = FALSE, size = 10, col = "red")
  points3d(digitspec[(fixed + 1):nrow(digitspec), ], aspect = F, size = 10, col = "green")}
  else { points3d(digitspec, aspect = F, size = 10, col = "green") }
}
