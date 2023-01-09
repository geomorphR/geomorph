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
#' This function is a wrapper for several functions in the \code{\link{rgl-package}} package.  
#' Although there is some allowance for arguments to be passed to \code{\link{rgl-package}} functions,
#' some override of rgl-plot3d arguments is required.  Errors that result from trying
#' to pass rgl-plot3d or rgl-points3d arguments should inspire the user
#' to find solutions with \code{\link{rgl-package}} core functions.
#'
#' @param spec An object of class shape3d/mesh3d, or matrix of 3D vertex coordinates.
#' @param digitspec Name of data matrix containing 3D fixed and/or surface sliding coordinates.
#' @param fixed Numeric The number of fixed template landmarks (listed first in {digitspec})
#' @param fixed.pt.col The color for plotting fixed template landmarks (if any)
#' @param fixed.pt.size The size for plotting fixed template landmarks (if any)
#' @param mesh.ptsize Numeric Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes                                       
#' @param centered Logical Whether the data matrix is in the surface mesh coordinate system ({centered = FALSE}) or
#' if the data were collected after the mesh was centered ({centered = TRUE})- see details.
#' @param ... additional parameters which will be passed to rgl-plot3d or 
#' rgl-points3d.
#' @export
#' @keywords visualization
#' @seealso \code{\link{warpRefMesh}}
#' @seealso \code{\link{read.ply}}
#' @seealso \code{\link{rgl-package}} (used in 3D plotting)
#' @examples
#' 
#' data(scallopPLY)
#' ply <- scallopPLY$ply
#' digitdat <- scallopPLY$coords
#' plotspec(spec = ply, digitspec = digitdat, fixed = 16, 
#' centered = TRUE, fixed.pt.col = "red", 
#' fixed.pt.size = 15, col = "blue", size = 5)
#' @author Erik Otarola-Castillo, Emma Sherratt, Antigoni Kaliontzopoulou, & Michael Collyer

plotspec <- function (spec, digitspec, fixed = NULL, fixed.pt.col = "red", fixed.pt.size = 10, 
                      mesh.ptsize = 1, centered = FALSE, ...) {
  dots <- list(...)
  if(!is.null(dots$type)) {
    cat("Warning: plot3d argument 'type' cannot be varied with this function.
         \nIf you wish to change the plot3d type, you should use the plot3d function, directly.\n\n")
    dots[names(dots) == "type"] <- NULL
  }
    
  plot3d.args <- list(x = 0, y = 0, z = 0, 
                      xlab = NULL, ylab = NULL, zlab = NULL,
                      col = 1, size = mesh.ptsize,
                      lwd = 1, radius = 0, add = FALSE, aspect = FALSE,
                      xlim = NULL, ylim = NULL, zlim = NULL, 
                      forceClipregion = FALSE, box = TRUE,
                      axes = TRUE, main = NULL, sub = NULL, top = TRUE,
                      expand = 1.03)
  
  m.p <- match(names(dots), names(plot3d.args))
  
  if(any(is.na(m.p)))
    stop("Some of the arguments in ... are arguments that cannot be passed to plot3d or points 3d.\n",
         call. = FALSE)
  
  plot3d.args[m.p] <- dots
  plot3d.args$type <- "n"
  
  if (inherits(spec, "shape3d") || inherits(spec, "mesh3d")){
    
    if(centered){
      specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    
    if(!centered) specimen <- as.matrix(t(spec$vb)[,-4])
    
    mesh <- spec
    if (is.null(mesh$material)) mesh$material$color <- "gray" 
    if (is.null(mesh$material$color)) mesh$material$color <- "gray" 
  } 
  
  else if (!inherits(spec, "matrix")) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
    
  } 
  
  else if (inherits(spec, "matrix") && dim(spec)[2] == 3) {
     specimen <- if(centered) scale(spec, scale = FALSE) else
       spec
     
  } else stop ("File is not matrix in form: vertices by xyz.\n", call. = FALSE)
  
  if (is.null(dim(digitspec)) || dim(digitspec)[2] != 3) {
    stop("Digitized file is not xyz matrix in form: p x k.\n", call. = FALSE)}
  
  plot3d.args$x <- specimen[,1]
  plot3d.args$y <- specimen[,2]
  plot3d.args$z <- specimen[,3]
  
  if(!is.null(colnames(specimen))) {
    dn <- colnames(specimen)
    if(is.null(plot3d.args$xlab)) plot3d.args$xlab <- dn[1]
    if(is.null(plot3d.args$ylab)) plot3d.args$ylab <- dn[2]
    if(is.null(plot3d.args$zlab)) plot3d.args$zlab <- dn[3]
  }
  
  if(is.null(plot3d.args$xlab)) plot3d.args$xlab <- "x"
  if(is.null(plot3d.args$ylab)) plot3d.args$ylab <- "y"
  if(is.null(plot3d.args$zlab)) plot3d.args$zlab <- "z"
  
  do.call(plot3d, plot3d.args)

  plot3d.args <- plot3d.args[c("x", "y", "z", "col", "size")]

  do.call(points3d, plot3d.args)
  
  if(!is.null(mesh)) shade3d(mesh, meshColor = "legacy", add = TRUE)
  
  if(!is.null(fixed)) {
    p2.args <- xyz.coords(digitspec[1:fixed,])
    p2.args$size = fixed.pt.size
    p2.args$col = fixed.pt.col
    do.call(points3d, p2.args)
    
    if(nrow(digitspec) > fixed) {
      p2.args$x <- digitspec[(fixed + 1):nrow(digitspec), 1]
      p2.args$y <- digitspec[(fixed + 1):nrow(digitspec), 2]
      p2.args$z <- digitspec[(fixed + 1):nrow(digitspec), 3]
      p2.args$size <- plot3d.args$size
      p2.args$col <- plot3d.args$col
      do.call(points3d, p2.args)
    } 
  }

}
