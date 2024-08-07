#' Read landmark data matrix from fcsv file
#'
#' Read landmark data for a single specimen from an fcsv file obtained from SlicerMorph. 
#' 
#' This function merely extracts x, y, z coordinates from an fcsv file.  
#'
#' @param file The name of a *.fcsv file containing three-dimensional landmark data to be read in
#' @param header_no The number of header lines in the fscv file, needed for knowing these lines are not data.
#' @keywords IO
#' @export
#' @author Murat Maga, Michael Collyer and Dean Adams
#' @return Function returns a p x 3 matrix of x, y, z coordinates for p landmarks.

readland.fcsv = function (file = NULL, header_no = 3)  {
  testfile <- scan(file=file, what="char", quote="", sep="\n", strip.white=TRUE, comment.char="\"", quiet=TRUE)
  header_no <- length(grep("#", testfile))
  res <- read.csv(file = file, skip = header_no, header = F)[, 1:4]
  res.ind <- res[, 1]
  res.nms <- try(
    unlist(strsplit(res.ind, "vtkMRMLMarkups")), silent = TRUE)
  if(inherits(res.nms, "try-error")) res.nms <- res.ind else
    res.nms <- res.nms[seq(2, length(res.nms), 2)]
  res <- res[,-1]
  rownames(res) <- res.nms
  return(as.matrix(res))
}
