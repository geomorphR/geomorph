#' Set up parameters for grids, points, and links in plotRefToTarget
#'
#' Function produces a list of parameter changes within \code{\link{plotRefToTarget}}.
#'
#' The function allows users to vary certain plotting parameters to produce different
#' graphical outcomes for \code{\link{plotRefToTarget}}.  Not all parameters need to be adjusted to use this function, as the defaults above will be used.
#' 
#'   
#' @param pt.bg Background color of reference configuration points (single value or vector of values)
#' @param pt.size Scale factor for reference configuration points (single value or vector of values)
#' @param link.col The color of links for reference configurations (single value or vector of values)
#' @param link.lwd The line weight of links for reference configurations (single value or vector of values)
#' @param link.lty The line type of links for reference configurations (single value or vector of values)
#' @param out.col The color of outline for reference configurations (single value or vector of values)
#' @param out.cex The size of plotting symbol of outline for reference configurations (single value or vector of values)
#' @param tar.pt.bg Background color of target configuration points (single value or vector of values)
#' @param tar.pt.size Scale factor for target configuration points (single value or vector of values)
#' @param tar.link.col The color of links for target configurations (single value or vector of values)
#' @param tar.link.lwd The line weight of links for target configurations (single value or vector of values)
#' @param tar.link.lty The line type of links for target configurations (single value or vector of values)
#' @param tar.out.col The color of outline for target configurations (single value or vector of values)
#' @param tar.out.cex The size of plotting symbol of outline for target configurations (single value or vector of values)
#' @param n.col.cell The number of square cells (along x axis) for grids (single numerical value)
#' @param grid.col The color of grid lines (single value)
#' @param grid.lwd Scale factor for the weight of grid lines (single numerical value)
#' @param grid.lty The line type for grid lines (single numerical value, as in base R \code{\link{plot}})
#' @param txt.adj The adjustment value of the landmark label (one or two values, as in base R \code{\link{text}}) 
#' @param txt.pos The position of the landmark label (single numerical value, as in base R \code{\link{text}}) 
#' @param txt.cex The size of the landmark label text (single numerical value, as in base R \code{\link{text}})
#' @param txt.col The color of the landmark label text (single numerical value, as in base R \code{\link{text}})
#' @keywords utilities
#' @keywords visualization
#' @export
#' @author Michael Collyer & Emma Sherratt
#' @seealso  \code{\link{plotRefToTarget}}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' ref<-mshape(Y.gpa$coords)
#' plotRefToTarget(ref,Y.gpa$coords[,,39]) # default settings
#' 
#' # Altering points and links
#' GP1 <- gridPar(pt.bg = "red", pt.size = 1, link.col="blue", link.lwd=2, n.col.cell=50)
#' plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP1, mag=2, 
#' links=plethodon$links, method="TPS")
#' 
#' # Altering point color
#' GP2 <- gridPar(pt.bg = "green", pt.size = 1) 
#' plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP2, mag=3, method="vector")
#' 
#' # Altering ref and target points
#' GP3 <- gridPar(pt.bg = "blue", pt.size = 1.5, tar.pt.bg = "orange", tar.pt.size = 1) 
#' plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP3, mag=3, method="points")
#' 
#' # Altering outline color
#' GP4 <- gridPar(tar.out.col = "red", tar.out.cex = 0.3) 
#' plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP4, mag=3, 
#' outline=plethodon$outline, method="TPS")
#' 
#' # Altering text labels
#' GP5 <- gridPar(txt.pos = 3, txt.col = "red") 
#' plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP5, mag=3, method="vector", label=TRUE)
gridPar <- function(pt.bg = "gray", 
                    pt.size = 1.5,
                    link.col = "gray",
                    link.lwd = 2,
                    link.lty = 1,
                    out.col = "gray",
                    out.cex = 0.1,
                    tar.pt.bg = "black",
                    tar.pt.size = 1,
                    tar.link.col = "black",
                    tar.link.lwd = 2,
                    tar.link.lty = 1,
                    tar.out.col = "black",
                    tar.out.cex = 0.1,
                    n.col.cell = 20,
                    grid.col = "black",
                    grid.lwd = 1,
                    grid.lty = 1,
                    txt.adj = 0.5,
                    txt.pos = 1, 
                    txt.cex = 0.8,
                    txt.col = "black"
){
  list(pt.bg=pt.bg,pt.size=pt.size,n.col.cell=n.col.cell,
       grid.col=grid.col,grid.lwd=grid.lwd,grid.lty=grid.lty,
       txt.adj=txt.adj,txt.pos=txt.pos,txt.cex=txt.cex,txt.col=txt.col,
       link.col=link.col,link.lwd=link.lwd,link.lty=link.lty,
       out.col=out.col,out.cex=out.cex,
       tar.pt.bg=tar.pt.bg,tar.pt.size=tar.pt.size,
       tar.link.col=tar.link.col,tar.link.lwd=tar.link.lwd,tar.link.lty=tar.link.lty,
       tar.out.col=tar.out.col,tar.out.cex=tar.out.cex)
}