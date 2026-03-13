#' Plot shape differences between a reference and target specimen
#'
#' Function plots shape differences between a reference and target specimen
#'
#' The function generates a plot of the shape differences of a target specimen relative to a reference 
#'  specimen. The option "mag" allows the user to indicates the degree of magnification to be used when 
#'  displaying the shape difference. The function will plot either two- or three-dimensional data. 
#'  
#'  For two-dimensional data and thin-plate spline deformation plots, the user may also supply boundary 
#'  curves of the object, which will be deformed from the reference to the target specimen using the 
#'  thin-plate spline. Such curves are often useful in describing 
#'  the biological shape differences expressed in the landmark coordinates.  Note that to utilize this option, 
#'  a boundary curve from a representative specimen must first be warped to the reference specimen using
#'   \code{\link{warpRefOutline}}.
#'   
#'  Additionally, if a matrix of links is provided, the landmarks will be connected by lines.  
#'  The link matrix is an M x 2 matrix, where M is the desired number of links. Each row of the link matrix 
#'  designates the two landmarks to be connected by that link.
#'   
#'  Four distinct methods for plots are available:
#'  \enumerate{
#'  \item {TPS} a thin-plate spline deformation grid is generated. For 3D data, 
#'  this method will generate thin-plate spline deformations in the x-y and x-z planes. 
#'  \item {vector}: a plot showing the vector displacements between corresponding landmarks in the reference 
#'  and target specimen is shown. 
#'  \item {points} a plot is displayed with the landmarks in the target overlaying 
#'  those of the reference.  
#'  \item {surface} a mesh3d surface is warped using thin-plate spline (for 3D data only). 
#'  Requires mesh3d object in option "mesh", made using \code{\link{warpRefMesh}}. 
#'  }
#'  This function combines numerous plotting functions found in Claude (2008).
#'
#' @param M1 Matrix of landmark coordinates for the first (reference) specimen
#' @param M2 Matrix of landmark coordinates for the second (target) specimen
#' @param mesh A mesh3d object for use with method = "surface"
#' @param outline An x,y curve or curves warped to the reference (2D only)
#' @param method Method used to visualize shape difference; see below for details
#' @param mag The desired magnification to be used when visualizing the shape difference (e.g., mag = 2)
#' @param links An optional matrix defining for links between landmarks
#' @param label A logical value indicating whether landmark numbers will be plotted
#' @param axes A logical value indicating whether the box and axes should be plotted (points and vector only)
#' @param gridPars An optional object made by \code{\link{gridPar}}
#' @param useRefPts An option (logical value) to use reference configuration points rather than target configuration points (when method = "TPS")
#' @param ... Additional parameters not covered by \code{\link{gridPar}} to be passed to plotting routines
#' or \code{\link[rgl]{shade3d}}
#' @return If using method = "surface", function will return the warped mesh3d object.
#' @keywords visualization
#' @export
#' @author Dean Adams, Emma Sherratt, Antigoni Kaliontzopoulou & Michael Collyer
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York. 
#' @seealso  \code{\link{gridPar}}
#' @seealso  \code{\link{define.links}}
#' @seealso  \code{\link{warpRefMesh}}
#' @seealso  \code{\link{warpRefOutline}}
#' @examples
#' \dontrun{
#' 
#' # Two dimensional data
#'  data(plethodon) 
#'  Y.gpa <- gpagen(plethodon$land)    #GPA-alignment
#'  ref <- mshape(Y.gpa$coords)
#'  plotRefToTarget(ref, Y.gpa$coords[,,39])
#'  plotRefToTarget(ref, Y.gpa$coords[,,39], mag = 2, outline = plethodon$outline)   
#' 
#' #magnify by 2X
#'  plotRefToTarget(ref, Y.gpa$coords[,,39], method = "vector", mag = 3)
#'  plotRefToTarget(ref, Y.gpa$coords[,,39], method = "points", 
#'  outline = plethodon$outline)
#'  plotRefToTarget(ref, Y.gpa$coords[,,39], method = "vector", 
#'  outline = plethodon$outline, mag = 2.5)
#'  plotRefToTarget(ref, Y.gpa$coords[,,39], 
#'  gridPars = gridPar(pt.bg = "green", pt.size = 1),
#'  method = "vector", mag = 3)
#'
#' # Three dimensional data
#'  data(scallops)
#'  Y.gpa <- gpagen(A = scallops$coorddata, curves = scallops$curvslide, 
#'  surfaces = scallops$surfslide)
#'  ref <- mshape(Y.gpa$coords)
#'  plotRefToTarget(ref, Y.gpa$coords[,,1], method = "points")
#'  scallinks <- matrix(c(1,rep(2:16, each=2),1), nrow = 16, byrow = TRUE)
#'  plotRefToTarget(ref, Y.gpa$coords[,,1],
#'  gridPars = gridPar(tar.pt.bg = "blue", tar.link.col="blue",
#'  tar.link.lwd = 2), method = "points", links = scallinks)
#' }
plotRefToTarget<-function(M1, M2, mesh= NULL, outline=NULL, 
                          method=c("TPS","vector","points","surface"),
                          mag=1.0, links=NULL, label=FALSE, axes=FALSE, 
                          gridPars=NULL, useRefPts=FALSE,...){
  method <- match.arg(method)
  if(any(is.na(M1))){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(any(is.na(M2))){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(is.null(gridPars)) gP <- gridPar() else gP <- gridPars
  k <- dim(M1)[2]
  p <- dim(M1)[1]
  mag <- (mag-1)
  M2 <- M2 + (M2-M1)*mag
  limits <- function(x,s){ 
    r <- range(x)
    rc <- scale(r, scale=FALSE)
    l <- mean(r)+s*rc
  }
  if(k==2){
    if(method=="TPS"){
      tps(M1, M2, gP$n.col.cell, sz=gP$tar.pt.size, pt.bg=gP$tar.pt.bg, grid.col=gP$grid.col, 
          grid.lwd=gP$grid.lwd, grid.lty=gP$grid.lty, refpts=useRefPts)
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1], M2[links[i,1],2], M2[links[i,2],1], M2[links[i,2],2],
                   col=linkcol[i], lty=linklty[i], lwd=linklwd[i])
        }
      }
      if(label){
        text(M2, label=paste(1:dim(M2)[1]), adj=gP$txt.adj,
                             pos=gP$txt.pos, cex=gP$txt.cex, col=gP$txt.col)
        }
      if(!is.null(outline)){
        curve.warp <- xy.coords(tps2d(outline, M1, M2))
        plot.xy(curve.warp, type="p", pch=19, cex=gP$tar.out.cex, col=gP$tar.out.col) 
      }
      if(!useRefPts){
        plot.xy(xy.coords(M2), type="p", pch=21, cex=gP$tar.pt.size, bg=gP$tar.pt.bg)
      } else {
        plot.xy(xy.coords(M1), type="p", pch=21, cex=gP$pt.size, bg=gP$pt.bg)}
    }
    if(method=="vector"){
      plot.new()
      if(axes){
        plot.window(limits(M1[,1], 1.25), limits(M1[,2], 1.25),
                    xlab="x", ylab="y", asp = 1)
      }
      if(!axes){
        plot.window(limits(M1[,1], 1.25), limits(M1[,2],1.25),
                    xlab="", ylab="", asp = 1, xaxt="n", yaxt="n")
      }
      if(label){
        text(M1, label=paste(1:dim(M1)[1]),adj=gP$txt.adj,
             pos=gP$txt.pos, cex=gP$txt.cex, col=gP$txt.col)
      }
      if(!is.null(outline)){
        curve.warp <- tps2d(outline, M1, M2)
        plot.xy(xy.coords(outline), type="p", pch=19, cex=gP$out.cex, col=gP$out.col) 
        plot.xy(xy.coords(curve.warp), type="p", pch=19, cex=gP$tar.out.cex, col=gP$tar.out.col) 
      }
      if(!is.null(links)){
        linkcol <- rep(gP$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1], M2[links[i,1],2], M2[links[i,2],1], M2[links[i,2],2],
                   col=linkcol[i], lty=linklty[i], lwd=linklwd[i])
        }
      }
      
      arrows(M1[,1], M1[,2], M2[,1], M2[,2], length=0.075,lwd=2)
      plot.xy(xy.coords(M1), type="p", pch=21, bg=gP$pt.bg, cex=gP$pt.size)
    }
    if(method=="points"){
      plot.new()
      if(axes){
        plot.window(limits(M1[,1], 1.25), limits(M1[,2], 1.25),
                    xlab="x", ylab="y", asp = 1)
      }
      if(!axes){
        plot.window(limits(M1[,1], 1.25), limits(M1[,2],1.25),
                    xlab="", ylab="", asp = 1, xaxt="n", yaxt="n")
      }
      if(label){
        text(M1, label=paste(1:dim(M1)[1]), adj=gP$txt.adj,
             pos=gP$txt.pos, cex=gP$txt.cex, col=gP$txt.col)
        }
      if(!is.null(outline)){
        curve.warp <- tps2d(outline, M1, M2)
        plot.xy(xy.coords(outline), type="p", pch=19, cex=gP$out.cex, col=gP$out.col) 
        plot.xy(xy.coords(curve.warp), type="p", pch=19, cex=gP$tar.out.cex, col=gP$tar.out.col) 
      }
      if(!is.null(links)){
        linkcol <- rep(gP$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty,nrow(links))[1:nrow(links)]
        tarlinkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        tarlinklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M1[links[i,1],1],M1[links[i,1],2],
                   M1[links[i,2],1],M1[links[i,2],2],col=linkcol[i],
                   lty=linklty[i],lwd=linklwd[i])
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],
                   M2[links[i,2],2],col=tarlinkcol[i],
                   lty=tarlinklty[i],lwd=tarlinklwd[i])
        }
      }
      plot.xy(xy.coords(M2), type="p", pch=21, bg=gP$tar.pt.bg, cex=gP$tar.pt.size)
      plot.xy(xy.coords(M1), type="p", pch=21, bg=gP$pt.bg, cex=gP$pt.size)
    }
    if(method=="surface"){
      stop("Surface plotting for 3D landmarks only.")
    }      
  }
  if(k==3){
    if(method=="TPS" && any(class(M2) == "matrix")){
      old.par <- par(no.readonly = TRUE)
      par(mfrow = c(1,2))
      par(mar=c(1,1,1,1))
      tps(M1[,1:2],M2[,1:2],gP$n.col.cell, sz=gP$tar.pt.size, pt.bg=gP$tar.pt.bg, grid.col=gP$grid.col, 
          grid.lwd=gP$grid.lwd, grid.lty=gP$grid.lty, refpts=useRefPts)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          linkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
          linklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
          linklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
          segments(M2[links[i,1],1],M2[links[i,1],2],
                   M2[links[i,2],1],M2[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      if (label) {
        text(M2, label = paste(1:dim(M2)[1]), adj = gP$txt.adj, 
             pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
      }
      title("X,Y tps grid")
      b<-c(1,3)
      tps(M1[,b],M2[,b],gP$n.col.cell, sz=gP$tar.pt.size, pt.bg=gP$tar.pt.bg, grid.col=gP$grid.col, 
          grid.lwd=gP$grid.lwd, grid.lty=gP$grid.lty, refpts=useRefPts)
      if(is.null(links)==FALSE){
        linkcol <- rep(gP$tar.link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],3],
                   M2[links[i,2],1],M2[links[i,2],3],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      if (label) {
        text(M2[,b], label = paste(1:dim(M2)[1]), adj = gP$txt.adj, 
             pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
      }
      title("Y,Z tps grid")
      par(mfrow = c(1,1))
      on.exit(par(old.par))
    }
    if(method!="TPS"){
      fig <- plot_ly()
      if(!axes){
        fig <- fig |>
          plotly::layout(scene = list(aspectmode = "data",
            xaxis = list(title = '', 
              showgrid = F, visible = F,showticklabels = F, 
              zeroline = F, showbackground = F),
            yaxis = list(title = '', showgrid = F, 
              visible = F, showticklabels = F, zeroline = F, 
              showbackground = F),
            zaxis = list(title = '', showgrid = F, 
              visible = F, showticklabels = F, zeroline = F, 
              showbackground = F)),showlegend = FALSE)
      }  
      if(axes){
        fig <- fig |>
          plotly::layout(scene = list(aspectmode = "data", showlegend = FALSE,
            xaxis = list(title = "X"),
            yaxis = list(title = "Y"),
            zaxis = list(title = "Z")))
      }
      if(method!="surface"){
        if(label == TRUE){ 
          fig <- fig |>
            add_trace(x = ~M1[,1], y = ~M1[,2], z = ~M1[,3],
                type = "scatter3d", mode = "text",
                text = seq_len(p), 
                textfont = list(color = gP$txt.col, size = gP$txt.cex),
                textposition = "top center",
                showlegend = FALSE,
                inherit = FALSE)
        }
        if(is.null(links)==FALSE){ 
          dash_map <- c("solid", "dash", "dot", "dashdot", 
              "longdash",  "longdashdot")
          line_df <- data.frame()
          for (i in 1:nrow(links)) {
            line_df <- rbind(
             line_df,
             M1[links[i, 1], ],
             M1[links[i, 2], ],
             data.frame(X = NA, Y = NA, Z = NA)
            )
          }
          fig <- fig |>
            add_trace(x = ~line_df$X, y = ~line_df$Y, z = ~line_df$Z,
              type = "scatter3d", mode = "lines",      
              line = list(color = gP$link.col, width = gP$link.lwd,
              dash = dash_map[gP$link.lty]),
              showlegend = FALSE,
              inherit = FALSE)
          if(method!="vector"){
            line_df2 <- data.frame()
            for (i in 1:nrow(links)) {
              line_df2 <- rbind(
                line_df2,
                M2[links[i, 1], ],
                M2[links[i, 2], ],
                data.frame(X = NA, Y = NA, Z = NA)
              )
            }
            fig <- fig |>
              add_trace(x = ~line_df2$X, y = ~line_df2$Y, z = ~line_df2$Z,
                type = "scatter3d", mode = "lines",      
                line = list(color = gP$tar.link.col, width = gP$tar.link.lwd,
                dash = dash_map[gP$tar.link.lty]),
                showlegend = FALSE,
                inherit = FALSE)
          }
        }  
      }
      if(method=="points"){
        fig <- fig |>
          add_trace(x = ~M1[,1], y = ~M1[,2], z = ~M1[,3], 
             type = "scatter3d", mode = "markers", 
             name = "Mean", marker = list(color = gP$pt.bg, 
             size = gP$pt.size*2),
             showlegend = FALSE,
             inherit = FALSE) 
        fig <- fig |>
          add_trace(x = ~M2[,1], y = ~M2[,2], z = ~M2[,3], 
             type = "scatter3d", mode = "markers", 
             name = "Mean", marker = list(color = gP$tar.pt.bg, 
             size = gP$tar.pt.size*2),
             showlegend = FALSE,
             inherit = FALSE) 
      }
      if(method=="vector"){
        fig <- fig |>
          add_trace(x = ~M1[,1], y = ~M1[,2], z = ~M1[,3], 
             type = "scatter3d", mode = "markers", 
             name = "Mean", marker = list(color = gP$pt.bg, 
             size = gP$pt.size*2),
             showlegend = FALSE,
             inherit = FALSE) 
        x <- as.vector(t(cbind(M1[,1], M2[,1], NA)))
        y <- as.vector(t(cbind(M1[,2], M2[,2], NA)))
        z <- as.vector(t(cbind(M1[,3], M2[,3], NA)))
        fig <- fig |>
          add_trace(type = "scatter3d", mode = "lines",
            x = x, y = y, z = z, 
            line = list(color = gP$tar.link.col, width = gP$tar.link.lwd),
            inherit = FALSE
          )
      }
      if(method=="surface"){
        if(is.null(mesh)){
          stop("Surface plotting requires a template mesh3d object (see 'warpRefMesh').")
        }
        warp.PLY <- mesh
        if (is.null(mesh$material$color)){mesh$material$color <- "gray"} 
        vb <- as.matrix(t(mesh$vb)[,-4])
        cat("\nWarping mesh\n")
        warp <- tps2d3d(vb, M1, M2)
        warp.PLY$vb <- rbind(t(warp), 1)
        fig <- fig |>
          add_trace(type = "mesh3d",
              x = warp.PLY$vb[1,], y = warp.PLY$vb[2,], z = warp.PLY$vb[3,],
              i = warp.PLY$it[1,] - 1, j = warp.PLY$it[2,] - 1, k = warp.PLY$it[3,] - 1,
              vertexcolor = t(grDevices::col2rgb(warp.PLY$material$color, alpha = TRUE)),
              inherit = FALSE)
      }
    fig
    }
  } 
}
