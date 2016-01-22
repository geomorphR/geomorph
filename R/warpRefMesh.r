#' Creates a mesh3d object warped to the mean shape
#' 
#' A function to take a 3D mesh and use thin-plate spline method to warp the file
#' into the estimated mean shape for a set of aligned specimens
#'
#' Function takes a 3D mesh (class mesh3d or shape3d, e.g. from \code{\link{read.ply}}) and its digitized landmark coordinates 
#' and uses the thin-plate spline method (Bookstein 1989) to warp the mesh into the shape 
#' defined by a second set of landmark coordinates, usually those of the 
#' mean shape for a set of aligned specimens. It is highly recommended that the mean shape is used as the 
#' reference for warping (see Rohlf 1998). The workflow is as follows:
#' \enumerate{
#' \item {Calculate the mean shape using \code{\link{mshape}}}
#' \item{Choose an actual specimen to use for the warping. The specimen used as the template for this warping 
#' is recommended as one most similar in shape to the average of the sample, but can be any reasonable 
#' specimen - do this by eye, or use \code{\link{findMeanSpec}} }
#' \item{Warp this specimen into the mean shape using \code{\link{warpRefMesh}} }
#' \item{Use this average mesh where it asks for a mesh= in the analysis functions and visualization functions  }
#' }
#' 
#' For landmark coordinates digitized with geomorph digitizing functions, centered = TRUE. This refers to the
#' specimen being centered prior to landmark acquisition in the RGL window. For landmark data collected outside
#' of geomorph, centered=FALSE will usually be the case. The returned mesh3d object is for use in geomorph
#' functions where shape deformations are plotted (\code{\link{plotTangentSpace}}, 
#' \code{\link{two.b.pls}}, \code{\link{bilat.symmetry}}, and \code{\link{plotRefToTarget}}). 
#' 
#' @param mesh A mesh3d object (e.g. made by \code{\link{read.ply}})
#' @param mesh.coord A p x k matrix of 3D coordinates digitized on the ply file.
#' @param ref A p x k matrix of 3D coordinates made by \code{\link{mshape}}
#' @param color Color to set the ply file $material. If the ply already has color, use NULL. 
#' For ply files without color, color=NULL will be plotted as grey.
#' @param centered Logical If the data in mesh.coords were collected from a centered mesh (see details).
#' @export
#' @seealso \code{\link{findMeanSpec}}
#' @keywords utilities
#' @keywords visualization
#' @author Emma Sherratt
#' @return Function returns a mesh3d object, which is a list of class mesh3d (see rgl for details)
#' @references  Bookstein, F. L. 1989 Principal Warps: Thin-Plate Splines and the Decomposition
#' of Deformations. IEEE Transactions on Pattern Analysis and Machine Intelligence 11(6):567-585.
#' @references  Rohlf, F. J. 1998. On Applications of Geometric Morphometrics to Studies of Ontogeny and Phylogeny. Systematic Biology. 47:147-158.
warpRefMesh <- function(mesh, mesh.coord, ref, color=NULL, centered=FALSE){
  if (inherits(mesh, "mesh3d") == FALSE){
    stop ("File is not a mesh3d object or xyz matrix") }
  open3d(); shade3d(mesh) ; title3d(main="Imported Mesh")
  mesh.vb <- as.matrix(t(mesh$vb)[,-4])
    if (centered == TRUE){ mesh.vb <- scale(mesh.vb, scale = FALSE) }
  checkmat <- is.matrix(mesh.coord)
    if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(mesh.coord)[2]
    if (checkdim==2) {stop("Input must be a p-x-k matrix of three-dimensional landmark coordinates") }
  coord <- scale(mesh.coord, scale=F) 
  sc.mat <- matrix(rep(1,nrow(mesh.vb)), ncol=1) %*% apply(mesh.coord,2,mean)
  mesh.vb <- mesh.vb - sc.mat 
  warp <- tps2d3d(mesh.vb, coord, ref)
  mesh$vb[1:3,] <- t(warp)
     if(is.null(color)==FALSE){ mesh$material <- color }
     if(is.null(color)==TRUE && is.null(mesh$material)==TRUE) { mesh$material <- "gray" }
  if(!is.null(mesh$normals)){ mesh <- addNormals(mesh)}
  open3d(); shade3d(mesh); title3d(main="Warped Ref Mesh")
  return(mesh)
}
