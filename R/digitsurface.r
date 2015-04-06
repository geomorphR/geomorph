#' Digitize 3D fixed landmarks and surface semilandmarks.
#' 
#' An interactive function to digitize three-dimensional (3D) landmarks on a surface lacking known landmarks.
#' Input for the function is either a matrix of vertex coordinates defining a 3D surface object
#' or a mesh3d object as obtained from \code{\link{read.ply}}. 
#' 
#' Function for digitizing fixed 3D landmarks and placing "surface sliders", semilandmarks that slide over a surface.
#'  Following selection of fixed points (see digitizing below), function finds surface semilandmarks following 
#'  algorithm outlined in Gunz et al. (2005) and Mitteroecker and Gunz (2009). digitsurface finds the same number of surface 
#'  semilandmarks as the template (created by \code{\link{buildtemplate}}) by downsampling scanned mesh, registering template with 
#'  current specimen via GPA. A nearest neighbor algorithm is used to match template surface landmarks to current specimen's. 
#'  To use function digitsurface, the template must be constructed first, and 'template.txt' be in the working directory. Because template 
#'  matching is based on the correspondence of fixed landmark points in the template and the specimen, a minimum of four fixed landmarks must be used. 
#' 
#' Some of the "fixed" landmarks digitized with digitsurface can be later designated as "curve sliders" using function \code{\link{define.sliders.3d}} 
#'  if required (see details in \code{\link{digit.fixed}}).
#'  NOTE: Function centers the mesh before digitizing by default (center=TRUE). If one chooses not to center,
#'  specimen may be difficult to manipulate in rgl window.
#' 
#' \subsection{Digitizing}{
#' Digitizing is interactive between landmark selection using a mouse (see below for instructions), 
#' and the R console. Once a point is selected, the user is asked if the system should keep or discard the 
#' selection #'(y/n). If "y", the user is asked to continue to select the next landmark. If "n" the removes the last chosen
#' landmark, and the user is asked to select it again. This can be repeated until the user is comfortable with the landmark
#' chosen. 
#' 
#' To digitize with a standard 3-button (PC):
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to select points to be digitized,
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. Macs using platform 
#' specific single button mice, XQuartz must be configured: go to Preferences > Input > tick "Emulate three button mouse":
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing COMMAND key to select vertex to be used as a landmark,
#'  \item press button while pressing OPTION key to adjust mesh perspective.
#'  \item the mouse SCROLLER or trackpad two finger scroll is used to zoom in an out.
#'  }
#' 
#' NOTE: there is no pan (translate) functionality in rgl library for all platforms at this time.
#' }
#'
#' @param spec Name of surface file, as either an object of class shape3d/mesh3d, or matrix of three-dimensional vertex coordinates.
#' @seealso \code{\link{buildtemplate}}
#' @seealso \code{\link{read.ply}} 
#' @seealso \code{\link{digit.fixed}} 
#' @param fixed numeric: The number of fixed template landmarks 
#' @param ptsize numeric: Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes
#' @param center Logical Whether the object 'spec' should be centered prior to digitizing (default {center=TRUE})
#' @return Function writes to the working directory an NTS file with the name of the specimen and .nts suffix containing the landmark coordinates.   
#' @references Gunz P, Mitteroecker P, & Bookstein FJ (2005) Semilandmarks in Three Dimensions. Modern Morphometrics in Physical Anthropology, ed Slice DE (Springer-Verlag, New York), pp 73-98.
#' @references Mitteroecker P & Gunz P (2009) Advances in Geometric Morphometrics. Evolutionary Biology 36(2):235-247. 
#' @export
#' @keywords digitizing
#' @author Erik Otarola-Castillo & Emma Sherratt
digitsurface<-function(spec, fixed, ptsize = 1, center = TRUE)    {
  if(fixed<4){stop ("Number of fixed points is not sufficient.")}
  spec.name<-deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    if (center == TRUE){
      specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    if (center == FALSE){ specimen <- as.matrix(t(spec$vb)[,-4]) }
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" }
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    if (center == TRUE){ specimen <- scale(spec, scale = FALSE) }
    if (center == FALSE){ specimen <- spec }
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  selected<-digit.fixed(spec, fixed,index=T,ptsize, center)
  template<-as.matrix(read.table("template.txt",header=TRUE))
  specimen<-trans(as.matrix(specimen))
  template<-trans(template)*(csize(specimen[selected$fix,])[[1]]/csize(template[(1:fixed),])[[1]])  
  template<-template%*%(pPsup(template[(1:fixed),],specimen[selected$fix,]))[[3]] 
  template.tps<-tps2d3d(template[-(1:fixed),],template[(1:fixed),],specimen[selected$fix,])             
  spec.surfs<-specimen[-selected$fix,]
  nei<-numeric(dim(template.tps)[1])
  sliders<-matrix(NA,nrow=dim(template.tps)[1],ncol=3)
  for (i in 1:dim(template.tps)[1])     {
    nei[i]<-which.min(sqrt((template.tps[i,1]-spec.surfs[,1])^2+(template.tps[i,2]-spec.surfs[,2])^2+(template.tps[i,3]-spec.surfs[,3])^2))[1] #3D NN
    sliders[i,]<-spec.surfs[nei[i],]
    spec.surfs<-spec.surfs[-nei[i],]  
  }
  clear3d();plot3d(specimen[,1],specimen[,2],specimen[,3],size=ptsize,aspect=FALSE)
  if (!is.null(mesh)) { 
    mesh$vb <- rbind(t(specimen), 1)
    shade3d(mesh, add=TRUE) 
  }
  points3d(specimen[selected$fix,],col="red",size=10)
  points3d(template.tps,col="blue",size=10)
  points3d(sliders[,1:3],size=10,col="green")
  selected.out <- rbind(specimen[selected$fix,],sliders)
  writeland.nts(selected.out, spec.name, comment=NULL)
}