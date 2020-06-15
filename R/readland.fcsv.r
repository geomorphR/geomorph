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
  res.nms <- res[, 1]
  res <- res[,-1]
  res.nms <- unlist(strsplit(res.nms, "vtkMRMLMarkups"))
  res.nms <- res.nms[seq(2, length(res.nms), 2)]
  rownames(res) <- res.nms
  return(as.matrix(res))
}
