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
#'  \subsection{Notes for geomorph 4.1}{ 
#'  Starting with geomorph version 4.1, interactive module selection in 3D is no longer
#'  supported, as RGL is not supported under MacOS Tahoe.
#'  }
#' 
#' @param spec A p x k matrix containing landmark coordinates of a single specimen (2D)
#' @param nmodules Number of modules to be defined
#' @return Function returns a vector of which landmarks belong in which module (e.g. 1, 1, 1, 2, 2, 3, 3, 3, 2) to be used
#' with \code{\link{modularity.test}} or \code{\link{integration.test}}.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{modularity.test}} and \code{\link{integration.test}} 
#' @author Dean Adams and Emma Sherratt
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
  if (checkdim == 3) {
    stop("3D no longer supported.") 
  } 
  return(as.vector(selected))
}
