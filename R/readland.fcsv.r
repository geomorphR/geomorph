#' Read landmark data matrix from fcsv file
#'
#' Read landmark data for a single specimen from an fcsv file obtained from SlicerMorph. 
#' 
#' This function merely extracts x, y, z coordinates from an fcsv file.  
#'
#' @param file the name of a *.fcsv file containing three-dimensional landmark data to be read in
#' @keywords IO
#' @export
#' @author Murat Maga and Michael Collyer
#' @return Function returns a p x 3 matrix of x, y, z coordinates for p landmarks.

readland.fcsv = function (file = NULL)  {
  
  res <- read.csv(file = file, skip = 2, header = T)[, 1:4]
  res.ind <- res[, 1]
  res.nms <- try(
    unlist(strsplit(res.ind, "vtkMRMLMarkups")), silent = TRUE)
  if(inherits(res.nms, "try-error")) res.nms <- res.ind else
    res.nms <- res.nms[seq(2, length(res.nms), 2)]
  res <- res[,-1]
  rownames(res) <- res.nms
  return(as.matrix(res))
}
