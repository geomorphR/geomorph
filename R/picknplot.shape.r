#' Pick points in geomorph scatterplots to visualize shape variation
#'
#' Function plots the shape corresponding to a clicked point in the area of a geomorph plot
#'
#' This function makes the scatterplots created by generic \code{\link{plot}} functions for geomorph analytical objects interactive, 
#' which allows one to visualize shape variation by selecting one or more points in morphospace. The function uses \code{\link{shape.predictor}} 
#' to estimate the shape corresponding to the selected point(s) based on the prediction underlying the scatterplot, and it plots the estimated 
#' shape as compared to the consensus landmark configuration using \code{\link{plotRefToTarget}}. The user is then prompted as to whether the plotted
#' shape is to be saved as a png file, in which case the name of the file needs to be provided (without quotation marks).
#' Interactive plots are at present available for plots produced by \code{\link{plot.gm.prcomp}}.
#' 
#' @param x a geomorph plot object of class plot.gm.prcomp, plot.procD.lm, plot.trajectory.analysis or plot.pls 
#' @param ... other arguments passed to \code{\link{plotRefToTarget}}
#' @return A list with the following components:
#' \item{points}{A list with the xy coordinates of the selected points.}
#' \item{shapes}{A list with the corresponding estimated shapes.}
#' @keywords visualization
#' @export
#' @author Antigoni Kaliontzopoulou & Emma Sherratt
#' @seealso  \code{\link{shape.predictor}}, \code{\link{plotRefToTarget}}
#' 
#' @examples
#' 2d
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)
#' pleth.pca <- gm.prcomp(Y.gpa$coords)
#' pleth.pca.plot <- plot(pleth.pca)
#' picknplot.shape(pleth.pca.plot) 
#' # May change arguments for plotRefToTarget
#' picknplot.shape(plot(pleth.pca), method = "vector", mag = 3)
#' 
#' # 2d with phylomorphospace
#' data(plethspecies) 
#' Y.gpa <- gpagen(plethspecies$land)
#' gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
#' pleth.phylomorpho <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
#' pleth.phylo.plot <- plot(pleth.phylomorpho, phylo = TRUE, cex = 2, pch = 22, bg = gps, 
#'      phylo.par = list(edge.color = "blue", edge.width = 2, edge.lty = 2,
#'      node.pch = 22, node.bg = "black")) 
#' picknplot.shape(pleth.phylo.plot, method = "points", mag = 5)
#' 
#' # 3d
#' data(scallops)
#' Y3d.gpa <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
#' scallops.pca <- gm.prcomp(Y3d.gpa$coords)
#' scallops.pca.plot <- plot(scallops.pca)
#' picknplot.shape(scallops.pca.plot) #May change method for plotRefToTarget

picknplot.shape <- function(x, ...){
  if(!class(x)%in%c("plot.gm.prcomp")){
    stop("Class of plot object not compatible with picknplot.shape. \n Please see the help file for allowed plot objects")
  }
  
  plot.args <- list(...)
  if(is.null(plot.args$method)) {plot.args$method <- "TPS"}
  eval(x$call)
  continue <- "y"
  p = 1
  picked.pts <- list()
  picked.shapes <- list()
  while(continue == "y"){
    cat("Pick a point in the shape space", "\n")
    picked.pts[[p]] <- unlist(locator(n = 1, type = "p", pch = 20, col = "red", cex = 1))
    cat("Picked point coordinates are:", "\n")
    cat(picked.pts[[p]], "\n")
    
    y <- get(paste(as.list(x$call)$x))
    A <- y$Y
    picked.shapes[[p]] <- shape.predictor(A, x = x$points[1:dim(A)[3],], pred1 = picked.pts[[p]])$pred1 
    if (dim(A)[2]==2) {
      plot.args$M1 <- cbind(mshape(A), 0)
      plot.args$M2 <- cbind(picked.shapes[[p]], 0)
      class(plot.args$M2) <- "predshape.k2"
      view3d(phi = 0, fov = 30, interactive = FALSE) 
      do.call(plotRefToTarget, plot.args)
    }
    if (dim(A)[2]==3){
      plot.args$M1 <- mshape(A)
      plot.args$M2 <- picked.shapes[[p]]
      class(plot.args$M2) <- "predshape.k3"
      if(plot.args$method == "TPS"){
        view3d(phi = 0, fov = 30, interactive = FALSE)
        } else {
          view3d(phi = 0, fov = 30, interactive = TRUE)
        }
      do.call(plotRefToTarget,  args = plot.args)
    }
    ans <- readline("Save deformation grid as png file (y/n)? ")
    if(ans=="y") {
      file.name <- readline("Please provide file name for saving deformation grid (without quotes) ")
      rgl.snapshot(file = file.name)
    }
    if(ans=="n"){
      try(rgl.close(), silent=T)
    }
    continue <- readline("Do you want to pick another point (y/n)? ")
    p = p + 1
  } 
  out <- list(points = picked.pts, shapes = picked.shapes)
  invisible(out)
}

  
  

