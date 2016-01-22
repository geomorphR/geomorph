#' Read mesh data (vertices and faces) from ply files
#'
#' A function to read ply files, which can be used for digitizing landmark coordinates or for shape warps.
#'
#' Function reads three-dimensional surface data in the form of a single ply file
#'  (Polygon File Format; ASCII format only, from 3D scanners such as NextEngine and David scanners). 
#'  Vertices of the surface may then be used to digitize three-dimensional points, 
#'  and semilandmarks on curves and surfaces. The surface may also be used as a mesh for visualizing 3D deformations (\code{\link{warpRefMesh}}).
#'  The function opens the ply file and plots the mesh,
#'  with faces rendered if file contains faces, and colored if the file contains vertex color.
#'  Vertex normals allow better visualization and more accurate digitizing with \code{\link{digit.fixed}}.
#'
#' @param file An ASCII ply file
#' @param ShowSpecimen logical Indicating whether or not the ply file should be displayed
#' @param addNormals logical Indicating whether or not the normals of each vertex should be calculated (using \code{\link[rgl]{addNormals}})
#' @export
#' @keywords IO
#' @author Dean Adams & Emma Sherratt
#' @return Function returns the following components:
#'   \item{mesh3d}{list of class mesh3d- see rgl for details}
#' @examples 
#' # If the file has no mesh color, or color is undesirable, user can assign this as follows:
#' # Using the example scallop PLY
#' data(scallopPLY) 
#' myply <- scallopPLY$ply
#' myply$material <- "gray" # using color word
#' myply$material <- "#FCE6C9" # using RGB code
read.ply <- function (file, ShowSpecimen = TRUE, addNormals = TRUE) 
{
  plyfile <- scan(file = file, what = "char", sep = "\n", strip.white = TRUE, 
                  quiet = TRUE)
  is.ply <- grep("ply", plyfile)
  if ((length(is.ply) ==0) ) stop ("File is not a PLY file")
  format <- unlist(strsplit(grep(c("format "), plyfile, value = TRUE), " "))
  if (format[2] != "ascii")
    stop("PLY file is not ASCII format: ","format = ", format[2:length(format)])
  face <- NULL
  material <- NULL
  xline <- unlist(strsplit(grep(c("element vertex "), plyfile, value = TRUE), " "))
  nvertices <- as.numeric(xline[grep(c("vertex"), xline) + 1])
  yline <- unlist(strsplit(grep(c("element face"), plyfile, value = TRUE), " "))
  nface <- as.numeric(yline[grep(c("face"), yline) + 1])
  headerend <- grep(c("end_header"), plyfile)
  ncolpts <- (length(grep(c("property float"), plyfile)))
  cols <- grep(c("property float"), plyfile, value = TRUE)
  x <- grep(c(" x"), cols)
  y <- grep(c(" y"), cols)
  z <- grep(c(" z"), cols)
  nuchar <- (length(grep(c("property uchar"), plyfile)))
  points <- as.matrix(as.numeric(unlist(strsplit(plyfile[(headerend + 
                                                            1):(headerend + nvertices)], " "))))
  points <- matrix(points, nrow = nvertices, byrow = T)
  xpts <- points[,x]
  ypts <- points[,y]
  zpts <- points[,z]
  vertices <- rbind(xpts, ypts, zpts, 1)
  if (yline[3] == 0) {print("Object has zero faces")}
  if (yline[3] != 0) {
    face <- as.matrix(as.numeric(unlist(strsplit(plyfile[(headerend + 
                                                            nvertices + 1):(headerend + nvertices + nface)], " "))))
    face <- t(matrix(face, nrow = nface, byrow = T))
    face <- face[2:4,]
    face = face +1 
    }
  if (nuchar !=0) {
    color <- rgb(points[,4], points[,5], points[,6], maxColorValue = 255)
    material$color <- matrix(color[face], dim(face))}   
  
  mesh <- list(vb = vertices, it = face, primitivetype = "triangle", 
               material = material)
  class(mesh) <- c("mesh3d", "shape3d")
  if(addNormals==TRUE){ mesh <- addNormals(mesh)}
  if(ShowSpecimen==TRUE){ 
    clear3d()
    if (length(face) == 0) { dot3d(mesh) }
    if (length(material) != 0){ shade3d(mesh) }
    shade3d(mesh, color="gray") }
  return(mesh)
}