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
#' @param x a geomorph plot object of class plot.gm.prcomp, plot.procD.lm, plot.pls, or plotAllometry 
#' @param ... other arguments passed to \code{\link{plotRefToTarget}}
#' @return A list with the following components:
#' \item{points}{A list with the xy coordinates of the selected points.}
#' \item{shapes}{A list with the corresponding estimated shapes.}
#' @keywords visualization
#' @export
#' @author Antigoni Kaliontzopoulou & Emma Sherratt
#' @seealso  \code{\link{shape.predictor}}, \code{\link{plotRefToTarget}}
#' @seealso  \code{\link[rgl]{rgl-package}} (used in 3D plotting)
#' @examples
#' 
#' ### Because picknplot requires user decisions, the following examples
#' ### are not run (but can be with removal of #s)
#' 
#' # 2d
#' # data(plethodon) 
#' # Y.gpa <- gpagen(plethodon$land)
#' # pleth.pca <- gm.prcomp(Y.gpa$coords)
#' # pleth.pca.plot <- plot(pleth.pca)
#' # picknplot.shape(pleth.pca.plot) 
#' # May change arguments for plotRefToTarget
#' # picknplot.shape(plot(pleth.pca), method = "points", mag = 3, links=plethodon$links)
#' 
#' # 2d with phylogeny
#' # data(plethspecies) 
#' # Y.gpa <- gpagen(plethspecies$land)
#' # gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
#' # pleth.phylo <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
#' # pleth.phylomorphospace <- plot(pleth.phylo, phylo = TRUE, cex = 2, pch = 22, 
#' # bg = gps, phylo.par = list(edge.color = "blue", edge.width = 2, edge.lty = 2,
#' # node.pch = 22, node.bg = "black"))
#' # links.species <- plethodon$links[-11,]
#' # links.species[11, 1] <- 11
#' # picknplot.shape(pleth.phylomorphospace, method = "points", links = links.species)
#' 
#' # 3d
#' # data(scallops)
#' # Y3d.gpa <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
#' # scallops.pca <- gm.prcomp(Y3d.gpa$coords)
#' # scallops.pca.plot <- plot(scallops.pca)
#' # picknplot.shape(scallops.pca.plot) 

picknplot.shape <- function(x, ...){
  if(!inherits(x, c("gm.prcomp", "plot.procD.lm", "plotAlometry", "plot.pls"))){
    stop("Class of plot object not compatible with picknplot.shape. \nPlease see the help file for allowed plot objects\n",
         call. = FALSE)
  }
  
  do.call(plot, x$plot.args)
  continue <- "y"
  p = 1
  picked.pts <- picked.shapes <- list()
  
  prt.arg.names <- formalArgs(plotRefToTarget)
  prt.arg.names <- prt.arg.names[- (length(prt.arg.names))]
  prt.args <- vector("list", length(prt.arg.names))
  names(prt.args) <- prt.arg.names
  dots <- list(...)
  if(length(dots) >= 1) {
    m <- match(prt.arg.names, names(dots), 0L)
    for(i in 1:length(m)) if(m[i] > 0) prt.args[i] <- dots[m[i]]
  }
  if(is.null(prt.args$method)) prt.args$method <- "TPS"
  if(is.null(prt.args$mag)) prt.args$mag <- 1
  if(is.null(prt.args$label)) prt.args$label <- FALSE
  if(is.null(prt.args$axes)) prt.args$axes <- FALSE
  if(is.null(prt.args$useRefPts )) prt.args$useRefPts <- FALSE
  
  if(!inherits(x, "plot.pls")) {
    type <- "LM"
    A1 <- x$GM$fitted
    if(is.null(A1)) stop("No shape data provided\n", call. = FALSE)
    A2 <- NULL
    if(!is.null(x$CAC)) A1 <- A1 + x$GM$residuals
  } else {
    type <- "PLS"
    A1 <- x$A1
    A2 <- x$A2
    if(length(dim(A1)) != 3) A1 <- NULL
    if(length(dim(A2)) != 3) A2 <- NULL
    if(is.null(A1) && is.null(A2)) stop("No shape data provided\n", call. = FALSE)
  }
  
  while(continue == "y"){
    cat("Pick a point in the shape space", "\n")
    picked.pts[[p]] <- unlist(locator(n = 1, type = "p", pch = 20, col = "red", cex = 1))
    cat("Picked point coordinates are:", "\n")
    cat(picked.pts[[p]], "\n")

  
    if(type == "LM") {
      X <- as.matrix(cbind(x$plot.args$x, x$plot.args$y))
      picked.shapes[[p]] <- shape.predictor(A1, X, 
                                            pred1 = picked.pts[[p]])$pred1 
    }

    if(type == "PLS") {
      X <- as.matrix(cbind(x$plot.args$x, x$plot.args$y))
      picked.shapes[[p]] <- list(P1 = shape.predictor(A1, X, 
                                                      pred1 = picked.pts[[p]])$pred1,
                                 P2 = shape.predictor(A2, X,  
                                                      pred1 = picked.pts[[p]])$pred1)
    }

    if(type == "LM") {
      
      if (dim(A1)[2] == 2) {
        prt.args$M1 <- cbind(mshape(A1), 0)
        prt.args$M2 <- cbind(picked.shapes[[p]], 0)
        class(prt.args$M2) <- "predshape.k2"
        view3d(phi = 0, fov = 30, interactive = FALSE) 
        do.call(plotRefToTarget, prt.args)
      }
      if (dim(A1)[2] == 3){
        prt.args$M1 <- mshape(A1)
        prt.args$M2 <- picked.shapes[[p]]
        class(prt.args$M2) <- "predshape.k3"
        if(prt.args$method == "TPS"){
          view3d(phi = 0, fov = 30, interactive = FALSE)
        } else {
          view3d(phi = 0, fov = 30, interactive = TRUE)
        }
        do.call(plotRefToTarget,  prt.args)
      }
      
    }
    
    if(type == "PLS") cat("Hold on for a moment...\n")

    ans <- readline("Save deformation grid as png file (y/n)? ")
    if(ans=="y") {
      file.name <- readline("Please provide file name for saving deformation grid (without quotes) ")
      rgl.snapshot(filename = file.name)
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

  
  

