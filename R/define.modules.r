#' Define modules (landmark partitions)
#' 
#' An interactive function to define which landmarks should be assigned to each module (landmark partition).
#' 
#' Function takes a matrix of digitized landmark coordinates (e.g. from \code{\link{mshape}}) and allows the user to assign 
#' landmarks to each module. The output is a list of which 
#' landmarks belong in which partition, to be used by \code{\link{modularity.test}} or \code{\link{integration.test}}.  
#'  
#'  \subsection{Selection in 2D}{ 
#' Choosing which landmarks will be included in each module involves landmark selection using a mouse in 
#' the plot window. The user is prompted to select each landmark in turn to be assigned to module 1: using the LEFT mouse button 
#' (or regular button for Mac users), click on the hollow circle to choose the landmark. Selected landmarks 
#' will be filled in. When all landmarks for module 1 are chosen, press 'esc', and then start selecting
#' landmarks for module 2. Repeat until all modules are defined.
#' }
#' 
#'  \subsection{Selection in 3D}{ 
#' Choosing which landmarks will be included in each module involves landmark selection using a mouse in 
#' the rgl plot window. The user is prompted to select one or more landmarks. To do so, use the RIGHT mouse button 
#' (or command + LEFT button for Mac users), draw a rectangle around landmarks to select.
#' Selected landmarks will be colored yellow. Then type into the console a letter (e.g. 1, 2, 3...) to assign selected landmark(s) 
#' to this module. Repeat until all landmarks are assigned to modules.
#' }
#' 
#' @param spec Name of specimen, as an object matrix containing 2D or 3D landmark coordinates
#' @param nmodules Number of modules to be defined
#' @return Function returns a vector of which landmarks belong in which module (e.g. 1,1,1,2,2,3,3,3,2) to be used
#' with \code{\link{modularity.test}} or \code{\link{integration.test}}.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{modularity.test}} and \code{\link{integration.test}} 
#' @author Emma Sherratt
#' 
define.modules <- function(spec, nmodules){
  spec.name <- deparse(substitute(spec))
  checkmat <- is.matrix(spec)
  if (checkmat == FALSE) {
    stop("Input must be a p-x-k matrix of landmark coordinates") }
  checkdim <- dim(spec)[2]
  if(nmodules > dim(spec)[1]){ 
    stop("Number of modules exceeds number of landmarks") }
  selected <- matrix(NA, nrow=nrow(spec), ncol=1)
  module <- cbind(c(1:nmodules), rainbow(nmodules))
  # 2D
  if (checkdim == 2) {
    plot(spec[, 1], spec[, 2], cex = 1, pch = 21, bg = "white", 
         xlim = range(spec[, 1]), ylim = range(spec[, 2]), asp = 1,
         xlab="x", ylab="y")
    text(spec[, 1], spec[, 2], label = paste(1:dim(spec)[1]), 
         adj = 0.5, pos = 1)
    for (i in 1:nmodules){
      cat("Select landmarks in module ",i,"\n",sep="")
      cat("Press esc when finished ","\n",sep="")
      select <- identifyPch(spec, col=module[i,2])
      selected[select] <- module[i,1] }
  }
  # 3D
  if (checkdim == 3) {
    plot3d(spec[, 1], spec[, 2],spec[, 3], size = 5,
         xlim = range(spec[, 1]), ylim = range(spec[, 2]), zlim = range(spec[, 3]), 
         asp = 1, box=F, axes=F, xlab="", ylab="", zlab="")
    text3d(spec[, 1], spec[, 2], spec[, 3], texts = paste(1:dim(spec)[1]), 
         adj = 1.3, pos = 4)
    rgl.bringtotop(stay = FALSE)
    while(anyNA(selected)==TRUE){
      cat("Select landmarks","\n")
      f<-select<-ans<-NULL
      f<-select3d(button="right")
      select<-f(spec)
      if(anyNA(spec[which(select==TRUE)[1],])==TRUE){
            cat("No vertex selected, try again","\n") 
            f<-select<-NULL
            f<-select3d(button="right")
            select<-f(spec)}
      points3d(spec[which(select==TRUE),1],spec[which(select==TRUE),2],spec[which(select==TRUE),3],
               size=12,color="yellow")
      cat(paste("Assign landmarks", paste(which(select==TRUE), collapse=","), "to which module?","\n")) 
      ans<-readLines(n=1)
      if(!is.null(ans)){
        selected[select]<-ans
        plot3d(spec[, 1], spec[, 2],spec[, 3], size = 5,
               xlim = range(spec[, 1]), ylim = range(spec[, 2]), zlim = range(spec[, 3]), 
               asp = 1, box=F, axes=F, xlab="", ylab="", zlab="")
        text3d(spec[, 1], spec[, 2], spec[, 3], texts = paste(1:dim(spec)[1]), 
               adj = 1.3, pos = 4)
        for(i in which(!is.na(selected))){
        points3d(spec[i,1],spec[i,2],spec[i,3],
                 size=10,color=module[which(selected[i]==module[,1]),2], add=T)}
      }
    }
  } 
  return(as.vector(selected))
}
