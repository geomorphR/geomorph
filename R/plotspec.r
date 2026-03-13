#' Plot 3D specimen, fixed landmarks and surface semilandmarks
#'
#' A function to plot three-dimensional (3D) specimen along with its landmarks.
#' 
#' Function to plot 3D specimens along with their digitized "fixed" landmarks and semilandmarks
#' "surface sliders" and "curve sliders". If specimen is a 3D surface (class shape3d/mesh3d) mesh is plotted.
#' Data collected outside geomorph should likely be read using centered = "FALSE". 
#' 3D coordinate data collected using geomorph versions prior to build. 1.1-6 were centered by 
#' default, and thus use centered = "TRUE". The function assumes the fixed landmarks are listed at the beginning of
#' the coordinate matrix  (digitspec).
#' 
#' @param spec An object of class shape3d/mesh3d, or matrix of 3D vertex coordinates.
#' @param digitspec Name of data matrix containing 3D fixed and/or surface sliding coordinates.
#' @param fixed Numeric The number of fixed template landmarks (listed first in digitspec)
#' @param fixed.pt.col The color for plotting fixed template landmarks (if any)
#' @param fixed.pt.size The size for plotting fixed template landmarks (if any)
#' @param mesh.pt.col The color for plotting mesh template landmarks
#' @param mesh.pt.size The size to plot the mesh template points    
#' @param centered Logical Whether the data matrix is in the surface mesh coordinate system (centered = "FALSE") or
#' if the data were collected after the mesh was centered (centered = "TRUE")- see details.
#' @param ... additional parameters which will be passed
#' @export
#' @keywords visualization
#' @seealso \code{\link{read.ply}}
#' @examples
#' \dontrun{
#' 
#'  data(scallopPLY)
#'  ply <- scallopPLY$ply
#'  digitdat <- scallopPLY$coords
#'  plotspec(spec = ply, digitspec = digitdat, fixed = 16, 
#'  centered = TRUE, fixed.pt.col = "red", 
#'  fixed.pt.size = 15, col = "blue", size = 5)
#'  }
#' @author Dean Adams, Antigoni Kaliontzopoulou, & Michael Collyer

plotspec <- function (spec, digitspec, fixed = NULL, 
                      fixed.pt.col = "red", fixed.pt.size = 4, 
                      mesh.pt.col = "black", mesh.pt.size = 2,
                      centered = FALSE, ...) {
  
  if(centered){
    specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
    spec$vb <- rbind(t(specimen), 1)
  }
  if (is.null(spec$material$color)){spec$material$color <- "gray"} 
  colnames(digitspec) <- c("X","Y","Z")
  
  fig <- plot_ly()
  fig <- fig |>
    add_trace(type = "mesh3d",
    x = spec$vb[1,], y = spec$vb[2,], z = spec$vb[3,],
    i = spec$it[1,] - 1, j = spec$it[2,] - 1, k = spec$it[3,] - 1,
    vertexcolor = t(grDevices::col2rgb(spec$material$color, alpha = TRUE)),
    inherit = FALSE
  ) |> 
    plotly::layout(scene = list(aspectmode = "data",
          xaxis = list(title = '', showgrid = F, visible = F,
                       showticklabels = F, zeroline = F, showbackground = F),
          yaxis = list(title = '', showgrid = F, visible = F, 
                       showticklabels = F, zeroline = F, showbackground = F),
          zaxis = list(title = '', showgrid = F, visible = F, 
                       showticklabels = F, zeroline = F, showbackground = F)),
          showlegend = FALSE)

  #Fixed points
  if(!is.null(fixed)){ 
  fig <- fig |>
    add_trace(x = ~digitspec$X[1:fixed], y = ~digitspec$Y[1:fixed], 
              z = ~digitspec$Z[1:fixed],
              type = "scatter3d", mode = "markers", 
              name = "fixed", marker = list(color = fixed.pt.col, 
                                           size = fixed.pt.size),
              inherit = FALSE) 
  }
  #Mesh points (sliders)
  if(is.null(fixed)){ fixed = 0}
  fig <- fig |>
    add_trace(x = ~digitspec$X[(fixed+1):nrow(digitspec)], 
              y = ~digitspec$Y[(fixed+1):nrow(digitspec)], 
              z = ~digitspec$Z[(fixed+1):nrow(digitspec)],
              type = "scatter3d", mode = "markers", 
              name = "sliding", marker = list(color = mesh.pt.col, 
                                            size = mesh.pt.size),
              inherit = FALSE) 
  fig
}