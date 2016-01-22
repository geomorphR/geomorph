#' Digitize 3D landmarks on mesh3d object
#'
#' An interactive function to digitize three-dimensional (3D) landmarks.
#' Input for the function is either a matrix of vertex coordinates defining a 3D surface object
#' or a mesh3d object as obtained from \code{\link{read.ply}}. 
#'
#' Function for digitizing "n" three-dimensional landmarks. The landmarks are "fixed" (traditional landmarks). They can be later
#'  designated as "curve sliders" (semilandmarks, that will "slide" along curves lacking known landmarks 
#'  if required. A sliding semi-landmark ("sliders") will slide between two designated points, along a line tangent to the specified curvature,
#'  and must be defined as "sliders" using function \code{\link{define.sliders}} or with similar format matrix made outside R. 
#' 
#' For 3D "surface sliders" (surface semilandmarks that slide over a surface) the function \code{\link{digitsurface}} should be used instead.
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
#' @param spec An object of class shape3d/mesh3d, or matrix of 3D vertex coordinates
#' @param fixed Numeric The number landmarks (fixed, and curve sliders if desired)
#' @param index Logical Whether selected landmark addresses should be returned (internal use only)
#' @param ptsize Numeric Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes 
#' @param center Logical Whether the object 'spec' should be centered prior to digitizing (default {center=TRUE})
#' @return Function writes to the working directory an NTS file with the name of the specimen and .nts suffix containing the landmark coordinates.
#' If index=FALSE function returns to the console an n x 3 matrix containing the x,y,z coordinates of the digitized landmarks.  
#' If index=TRUE, function returns a list:
#'  \item{selected}{a matrix containing the x,y,z coordinates of the digitized landmarks} 
#'  \item{fix}{a matrix of addresses for landmarks that are "fixed" (for internal use)}
#' @seealso \code{\link{read.ply}} 
#' @export
#' @keywords digitizing
#' @author Erik Otarola-Castillo & Emma Sherratt
digit.fixed <- function(spec, fixed, index=FALSE, ptsize = 1, center = TRUE)    {
  spec.name<-deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    if (center == TRUE){
      specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    if (center == FALSE){specimen <- as.matrix(t(spec$vb)[,-4])}
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" }
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    if (center == TRUE){ specimen <- scale(spec, scale = FALSE) }
    if (center == FALSE){ specimen <- spec }
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  clear3d(); ids <- plot3d(specimen[,1],specimen[,2],specimen[,3],size=ptsize,aspect=FALSE)
  if (!is.null(mesh)) { shade3d(mesh, add=TRUE) }
  selected<-matrix(NA,nrow=fixed,ncol=3);fix<-NULL    
  for (i in 1:fixed)      {
    keep<-ans<-NULL
    keep <- selectpoints3d(ids["data"], value= FALSE, button = "right")[2]
    points3d(specimen[keep,1],specimen[keep,2],specimen[keep,3],size=10,color="red",add=TRUE)
    selected[i,]<-as.numeric(specimen[keep,])
    fix<-c(fix,keep)   
    cat(paste("Keep Landmark ",i,"(y/n)?"),"\n")
    ans<-readLines(n=1)
    if(ans=="y" & length(fix)!=fixed) {
      cat("Select Landmark ",i+1,"\n")
      rgl.bringtotop(stay = FALSE)
    } 
    if(ans=="n" ) {
      cat(paste("Select Landmark ",i," Again"),"\n")
    }
    while(ans=="n") {
      selected[i,]<-NA
      fix<-fix[-i] 
      rgl.bringtotop(stay = FALSE)
      clear3d();ids <- plot3d(specimen[,1],specimen[,2],specimen[,3],size=ptsize,aspect=FALSE)
      if (!is.null(mesh)) { shade3d(mesh, add=TRUE) }
      if(sum(1-is.na(selected))>0){
        points3d(selected[,1],selected[,2],selected[,3],size=10,color="red",add=TRUE)
      }      
      keep<-ans<-NULL
      keep <- selectpoints3d(ids["data"], value= FALSE, button = "right")[2]
      points3d(specimen[keep,1],specimen[keep,2],specimen[keep,3],size=10,color="red",add=TRUE)
      selected[i,]<-as.numeric(specimen[keep,])
      fix<-c(fix,keep)   
      cat(paste("Keep Landmark ",i,"(y/n)?"),"\n")
      ans<-readLines(n=1)
      if(ans=="y" & length(fix)!=fixed) {
        cat("Select Landmark ",i+1,"\n")
      } 
      if(ans=="n") {
        cat(paste("Select Landmark ",i," Again"),"\n")
      }
    } 
  } 
  if(index==FALSE){
    writeland.nts(selected, spec.name, comment=NULL)
    rownames(selected)<-seq(from=1,to=nrow(selected))
    colnames(selected)<-c("xpts","ypts","zpts")
    return(selected) 
  }
  if(index==TRUE){
    return(list(selected=selected,fix=fix))
  }
}