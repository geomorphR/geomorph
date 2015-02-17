#' Define links between landmarks
#' 
#' An interactive function to define which landmarks should be linked to aid visualization.
#' 
#' Function takes a matrix of digitized landmark coordinates (e.g. from \code{\link{mshape}}) and 
#' allows the user to define pairs of landmarks to be linked, for visualization purposes. The output is a matrix 
#' to be used by \code{\link{plotAllSpecimens}} & \code{\link{plotRefToTarget}} option 'links='.  
#'  
#'  \subsection{Selection}{ 
#'  In the plot window select two landmarks that will be linked. In the console, the user will be prompted
#'  to continue linking landmarks to build the wireframe or end the session (typing y or n respectively).
#'  For 2D data in plot window, use LEFT mouse button to select landmarks.
#'  For 3D data in rgl window, use RIGHT mouse button (or command+LEFT for mac) to select landmarks.
#' }
#' 
#' @param spec Name of specimen, as an object matrix containing 2D or 3D landmark coordinates
#' @param ptsize Numeric Size to plot the landmarks
#' @return Function returns a matrix of which landmarks will be linkes
#' to be used by \code{\link{plotAllSpecimens}} & \code{\link{plotRefToTarget}} option 'links='.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{plotAllSpecimens}}
#' @seealso  \code{\link{plotRefToTarget}}
#' @author Emma Sherratt
#' 
define.links <- function(spec, ptsize=1){
  spec.name <- deparse(substitute(spec))
  checkmat <- is.matrix(spec)
  if (checkmat == FALSE) {
    stop("Input must be a p-x-k matrix of landmark coordinates") }
  checkdim <- dim(spec)[2]
  links <- NULL
  # 2D
  if (checkdim == 2) {
    plot(spec[, 1], spec[, 2], cex = ptsize, pch=21, 
         xlim = range(spec[, 1]), ylim = range(spec[, 2]), asp = 1,
         xlab="x", ylab="y")
    text(spec[, 1], spec[, 2], label = paste(1:dim(spec)[1]), 
         adj = 0.5, pos = 1)
    repeat{
            sel <- ans <- NULL
            cat("Select landmarks to link","\n")
              sel <- identify(spec, n=2, plot=FALSE)
              segments(x0=spec[sel[1],1], y0=spec[sel[1],2], x1=spec[sel[2],1], y1=spec[sel[2],2])
              links <- rbind(links, sel) 
            cat(paste("link made between landmarks", paste(sel, collapse=" & "), "\n", "Continue? type y/n"))
            ans <- readline()
            if(ans == "y"){ sel <- NULL}
            if(ans == "n"){ break}
            }
      }
  # 3D
  if (checkdim == 3) {
    if(is.null(size)){ size = 5 }
    ids <- plot3d(spec[, 1], spec[, 2],spec[, 3], size = ptsize,
         xlim = range(spec[, 1]), ylim = range(spec[, 2]), zlim = range(spec[, 3]), 
         asp = 1, box=F, axes=F, xlab="", ylab="", zlab="")
    text3d(spec[, 1], spec[, 2], spec[, 3], texts = paste(1:dim(spec)[1]), 
         adj = 1.3, pos = 4)
    repeat{
            rgl.bringtotop(stay = FALSE)
            cat("Select landmarks to link","\n")
            sel1 <- sel2 <- ans <- NULL
            sel1 <- selectpoints3d(ids["data"], value= FALSE, button = "right")[2]
            points3d(spec[sel1,1],spec[sel1,2],spec[sel1,3],size=10,color="yellow",add=TRUE)
            sel2 <- selectpoints3d(ids["data"], value= FALSE, button = "right")[2]
            points3d(spec[sel2,1],spec[sel2,2],spec[sel2,3],size=10,color="yellow",add=TRUE)
            segments3d(rbind(spec[sel1,], spec[sel2,]), lwd = 1)
            links <- rbind(links, c(sel1,sel2)) 
            cat(paste("link made between landmarks", sel1, "&", sel2, "\n", "Continue? type y/n"))
            ans <- readline()
            if(ans == "y"){ sel <- NULL
                            rgl.bringtotop(stay = FALSE)}
            if(ans == "n"){ break}
            }
    }
  return(as.matrix(links))
}