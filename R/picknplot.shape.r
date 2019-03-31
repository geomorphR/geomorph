#' Pick points in geomorph scatterplots to visualize shape variation
#'
#' Function plots the shape corresponding to a clicked point in the area of a geomorph plot
#' 
#' THIS FUNCTION IS A BIT EXPERIMENTAL!  
#' 
#' This function recycles scatter plots generated by \code{\link{plot.gm.prcomp}}, \code{\link{plot.procD.lm}}, 
#' \code{\link{plot.pls}}, or \code{\link{plotAllometry}}, and makes them interactive to visualize shape variation by selecting one or more points in morphospace. 
#' The function uses \code{\link{shape.predictor}} 
#' to estimate the shape corresponding to the selected point(s) based on the prediction underlying the scatterplot, and it plots the estimated 
#' shape as compared to the consensus landmark configuration using \code{\link{plotRefToTarget}}. The user is then prompted as to whether the plotted
#' shape is to be saved as a png file, in which case the name of the file needs to be provided (without quotation marks).
#' Interactive plots are at present available for plots produced by \code{\link{plot.gm.prcomp}}.  The function is limited in terms of the options for 
#' \code{\link{plotRefToTarget}} (because of the complexity of graphics); using \code{\link{shape.predictor}} and \code{\link{plotRefToTarget}}, directly, 
#' will always offer more flexibility.
#' 
#' IF YOU EXPERIENCE AN ERROR, please use \code{\link{shape.predictor}} and \code{\link{plotRefToTarget}}, directly.  (But please alert the 
#' geomorph package maintainer.)
#' 
#' 
#' @param x a geomorph plot object of class plot.gm.prcomp, plot.procD.lm, plot.pls, or plotAllometry 
#' @param ... other arguments passed to \code{\link{plotRefToTarget}}
#' @return A list with the following components:
#' \item{points}{A list with the xy coordinates of the selected points.}
#' \item{shapes}{A list with the corresponding estimated shapes.}
#' @keywords visualization
#' @export
#' @author Antigoni Kaliontzopoulou, Emma Sherratt, & Michael Collyer
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
#' # 2d allometry 
#' # gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
#' # species = plethodon$species) 
#' # fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, print.progress = FALSE)
#' # Predline
#' # PA <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", pch = 19)
#' # picknplot.shape(PA)
#' 
#' # 3d and two-b-pls
#' # data("scallops")
#' # Y.gpa <- gpagen(scallops$coorddata, curves = scallops$curvslide, 
#' #              surfaces = scallops$surfslide)
#' # PLS <- two.b.pls(Y.gpa$coords, Y.gpa$Csize)
#' # PLS.plot = plot(PLS)
#' # picknplot.shape(PLS.plot)

picknplot.shape <- function(x, ...){
  if(!inherits(x, c("plot.gm.prcomp", "plot.procD.lm", "plotAllometry", "plot.pls"))){
    stop("Class of plot object not compatible with picknplot.shape. \nPlease see the help file for allowed plot objects\n",
         call. = FALSE)
  }
  
  do.call(plot, x$plot.args)
  if(!is.null(x$phylo)){
    phylo.par <- x$phylo$phylo.par
    phy <- x$phylo$phy
    phy.pcdata <- x$phylo$phy.pcdata
    for (i in 1:nrow(phy$edge)) {
      dt.xy <- xy.coords(phy.pcdata[phy$edge[i,], ])
      plot.xy(dt.xy, type="l", col = phylo.par$edge.color, 
              lwd = phylo.par$edge.width, lty = phylo.par$edge.lty)
    }
    plot.xy(xy.coords(phy.pcdata[1:length(phy$tip),]), type="p",...)
    plot.xy(xy.coords(phy.pcdata[(length(phy$tip)+1):nrow(phy.pcdata),]), type="p",
            pch = phylo.par$node.pch, cex = phylo.par$node.cex, bg = phylo.par$node.bg)
    
  }
  
  continue <- "y"
  p = 1
  picked.pts <- picked.shapes <- list()
  prt.args <- list(...)
  if(is.null(prt.args$method)) prt.args$method <- "TPS"
  if(prt.args$method == "surface") {
    cat("method = 'surface' is not possible.",
        "\n Use shape.predictor and plotRefToTarget for this option.\n")
  }

  if(!inherits(x, "plot.pls")) {
    type <- "PC"
    if(is.null(x$GM$A)) A1 <- x$GM$fitted + x$GM$residuals else
      A1 <- x$GM$A
    if(is.null(A1)) stop("No shape data provided\n", call. = FALSE)
    A2 <- NULL
    
    if(!is.null(x$PredLine)) {
      if(!is.null(x$CAC)) {
        if(identical(x$plot.args$y, x$CAC)) {
          A1 <- A1 + x$GM$residuals
          type <- "regression2"
        }
      }
      if(identical(x$plot.args$y, x$PredLine) || identical(x$plot.args$y, x$RegScore))
        type <- "regression2"
      if(length(dim(A1)) != 3) stop("No shape data provided\n", call. = FALSE)
    }
   
  } else {
    type <- "PLS"
    A1 <- x$A1
    A2 <- x$A2
    if(length(dim(A1)) != 3) A1 <- NULL
    if(length(dim(A2)) != 3) A2 <- NULL
    if(is.null(A1) && is.null(A2)) stop("No shape data provided\n", call. = FALSE)
    if(is.null(A1) && !is.null(A2)) type <- "regression2"
    if(!is.null(A1) && is.null(A2)) type <- "regression1"
    if(type == "regression2") {
      A1 <- A2
      A2 <- NULL
    }
  }
  
  if(!is.null(A1) && dim(A1)[2] == 3) {
    cat("Only method = 'points' or 'vector' can be used for 3D data in this function.",
        "\nSwitching to method = 'points'.  Use shape.predictor and plotRefToTarget for more options.\n\n")
    prt.args$method <- "points"
  } 
  if(!is.null(A1) && dim(A1)[2] != 3 && !is.null(A2)  && dim(A2)[2] == 3) {
    cat("Only method = 'points' or 'vector' can be used for 3D data in this function.",
        "\nSwitching to method = 'points'.  Use shape.predictor and plotRefToTarget for more options.\n\n")
    prt.args$method <- "points"
  }
  
  while(continue == "y"){
    cat("Pick a point in the shape space", "\n")
    picked.pts[[p]] <- unlist(locator(n = 1, type = "p", pch = 20, col = "red", cex = 1))
    cat("Picked point coordinates are:", "\n")
    cat(picked.pts[[p]], "\n")

  
    if(type == "PC") {
      X <- as.matrix(cbind(x$plot.args$x, x$plot.args$y))
      picked.shapes[[p]] <- shape.predictor(A1, X, 
                                       pred1 = picked.pts[[p]])$pred1 
    }
    if(type == "regression2") {
      h <- picked.pts[[p]][2]
      abline(h = h, col = "red")
      X <- x$plot.args$y
      picked.shapes[[p]] <- shape.predictor(A1, X, 
                                            pred1 = h)$pred1 
    }
    
    if(type == "regression1") {
      v <- picked.pts[[p]][1]
      abline(v = v, col = "red")
      X <- x$plot.args$x
      picked.shapes[[p]] <- shape.predictor(A1, X, 
                                            pred1 = v)$pred1 
    }

    if(type == "PLS") {
      X <- as.matrix(cbind(x$plot.args$x, x$plot.args$y))
      picked.shapes[[p]] <- list(P1 = shape.predictor(A1, X, 
                                  pred1 = picked.pts[[p]])$pred1,
                                 P2 = shape.predictor(A2, X,  
                                   pred1 = picked.pts[[p]])$pred1)
    }

    if(type == "PC" || type == "regression1" || type == "regression2") {
      
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
    
    if(type == "PLS") {
      
      if (dim(A1)[2] == 2) {
        
        prt.args$M1 <- cbind(mshape(A1), 0)
        prt.args$M2 <- cbind(picked.shapes[[p]][[1]], 0)
        class(prt.args$M2) <- "predshape.k2"
        view3d(phi = 0, fov = 30, interactive = FALSE) 
        mfrow3d(1, 2)
        do.call(plotRefToTarget, prt.args)
        
      }
      if (dim(A1)[2] == 3){
        prt.args$main = "PLS Block 1"
        prt.args$M1 <- mshape(A1)
        prt.args$M2 <- picked.shapes[[p]][[1]]
        class(prt.args$M2) <- "predshape.k3"
        if(prt.args$method == "TPS"){
          open3d()
          mfrow3d(1, 2)
        } else {
          open3d()
          mfrow3d(1, 2)
        }
        do.call(plotRefToTarget,  prt.args)
        
      }
      
      if (dim(A2)[2] == 2) {
        prt.args$main = "PLS Block 2"
        prt.args$M1 <- cbind(mshape(A2), 0)
        prt.args$M2 <- cbind(picked.shapes[[p]][[2]], 0)
        class(prt.args$M2) <- "predshape.k2"
        do.call(plotRefToTarget, prt.args)
        
      }
      if (dim(A2)[2] == 3){
        prt.args$main = "PLS Block 2"
        prt.args$M1 <- mshape(A2)
        prt.args$M2 <- picked.shapes[[p]][[2]]
        class(prt.args$M2) <- "predshape.k3"
        do.call(plotRefToTarget,  prt.args)
        
      }
      
    }

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

  
  

