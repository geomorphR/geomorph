#' Build 3D surface template 
#'
#' An interactive function to build a template for the digitization across 
#' specimens of three dimensional (3D) surface sliding semilandmarks. Input 
#' for the function is either a matrix of vertex coordinates
#'  defining a 3D surface object or a mesh3d object as obtained 
#'  from \code{\link{read.ply}}.
#'
#' Function constructs a template of surface slider semilandmarks. If desired, the user can simultaneously 
#' digitize the fixed points (see digitizing below), however these may have been previously digitized separately 
#' using \code{\link{digit.fixed}}.
#' 
#' The function finds surface semilandmarks, chosen automatically from the mesh at roughly equal distances,
#' using the nearest-neighbor algorithm outlined in Gunz et al. (2005) and Mitteroecker and Gunz (2009). 
#' 
#' The set of fixed and surface slider landmarks are saved in the working directory as a txt file named "template". 
#' This file will then be used to extract a set of similarly numbered surface semilandmarks on subsequent specimens 
#' using the function \code{\link{digitsurface}}. Because template matching is based on the correspondence of fixed 
#' landmark points in the template and the target specimen, a minimum of four fixed landmarks must be used. However,  
#' to ensure a strong match between the scan and the template, it is recommended that a 
#' higher number of fixed points is used.
#' 
#' For more details see the vignette: \code{vignette("geomorph.digitize3D")}.
#' 
#'  NOTE: Function centers the mesh before digitizing by default (center = TRUE). If one chooses not to 
#'  center, the specimen may be difficult to manipulate in the rgl window.
#' 
#' \subsection{Digitizing}{
#' Digitizing of fixed landmarks is interactive. Once a point is selected, the user is asked if the system should keep or discard the 
#' selection (y/n). If "y", the user is asked to continue to select the next landmark. If "n" the removes the last chosen
#' landmark, and the user is asked to select it again. This can be repeated until the user is comfortable with the landmark
#' chosen. 
#'
#' To digitize with a standard 3-button (PC):
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to select points to be digitized,
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. 
#' Macs using platform specific single button mice, XQuartz must be configured: go to 
#' Preferences > Input > tick "Emulate three button mouse":
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing COMMAND key to select vertex to be used as a landmark),
#'  \item press button while pressing OPTION key to adjust mesh perspective.
#'  \item the mouse SCROLLER or trackpad two finger scroll is used to zoom in an out.
#'  }
#' 
#' NOTE: there is no pan (translate) functionality in rgl library for all platforms at this time.
#' The template can be edited using function \code{\link{editTemplate}}. 
#' }
#' 
#' \subsection{AUTO mode}{ 
#' The function as described above (for interactive mode) calls \code{\link{digit.fixed}}, prompting 
#' the user to select fixed landmarks in the rgl window. However if the user has digitized these fixed 
#' landmark elsewhere (e.g., in other software), then the input for parameter 'fixed' can be a 
#' p-x-k matrix of 3D coordinates. In this case, the function will automatically use these landmarks to build the 
#' template of sliding semilandmarks.
#' }
#'
#' @param spec An object of class shape3d/mesh3d, or matrix of 3D vertex coordinates
#' @param fixed Either a numeric value designating the number of fixed landmarks to be selected by \code{\link{digit.fixed}}, or a matrix of 3D coordinates collected previously
#' @param surface.sliders The number of desired surface sliders 
#' @param ptsize Size of mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes 
#' @param center Should the object 'spec' be centered prior to digitizing?
#' @seealso \code{\link{read.ply}} 
#' @seealso \code{\link{digit.fixed}} 
#' @seealso \code{\link{digitsurface}} 
#' @export
#' @keywords digitizing
#' @author Erik Otarola-Castillo & Emma Sherratt
#' @return The function writes to the working directory three files: an NTS file with the name of the 
#' specimen and .nts suffix containing the landmark coordinates, "template.txt" containing the same 
#' coordinates for use with the function \code{\link{digitsurface}}, and "surfslide.csv", a file 
#' containing the address of the landmarks defined as "surface sliders" for use with \code{\link{gpagen}}.
#' The function also returns to the console an n x 3 matrix containing the x,y,z coordinates of the 
#' digitized landmarks. 
#' @references Gunz P, Mitteroecker P, & Bookstein FJ (2005) Semilandmarks in Three Dimensions. Modern Morphometrics in Physical 
#' Anthropology, ed Slice DE (Springer-Verlag, New York), pp 73-98.
#' @references Mitteroecker P & Gunz P (2009) Advances in Geometric Morphometrics. Evolutionary Biology 36(2):235-247.

buildtemplate <- function(spec, 
                          fixed, 
                          center = TRUE,
                          surface.sliders, 
                          ptsize = 1) {
  if(length(fixed)==1 && fixed<4){stop ("Number of fixed points is not sufficient.")}
  spec.name<-deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    if (center == TRUE){
      specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    if (center == FALSE){ specimen <- as.matrix(t(spec$vb)[,-4]) }
    mesh <- spec 
    if (is.null(mesh$material)) mesh$material$color <- "gray" 
    if (is.null(mesh$material$color)) mesh$material$color <- "gray" 
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    if (center == TRUE){ specimen <- scale(spec, scale = FALSE) }
    if (center == FALSE){ specimen <- spec }
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  # Manual Mode
  if(length(fixed)==1){
    selected<-digit.fixed(spec,fixed,index=TRUE,ptsize, center)
    lmk.add<-selected$fix
    lmk.coords<-selected$selected }
  # Automatic mode
  if(inherits(fixed, "matrix") == TRUE){
    lmk.coords <- fixed
    lmk.add <- NULL
    for(i in 1:nrow(fixed)){
      lmk.add <- rbind(lmk.add, which.min(sqrt((fixed[i,1]-specimen[,1])^2+(fixed[i,2]-specimen[,2])^2+(fixed[i,3]-specimen[,3])^2))[1])}
    fixed <- nrow(fixed)} 
  surfs<-specimen[-lmk.add,]
  options(warn=-1)
    template<-rbind(lmk.coords,kmeans(x=surfs,centers=surface.sliders,iter.max=25)$centers)
  options(warn=0)
  points3d(template[-(1:fixed),1],template[-(1:fixed),2],template[-(1:fixed),3],size=10,col="blue")
  writeland.nts(template, spec.name, comment="Landmark coordinates for digitized template")
  write.table(template,file="template.txt",row.names=F,col.names=TRUE)
  write.csv(seq(from=(fixed+1),to=(surface.sliders+fixed)),file="surfslide.csv",row.names=FALSE)  
  rownames(template)<-seq(from=1,to=nrow(template))
  return(template)
}
