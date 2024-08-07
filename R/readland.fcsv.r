#' Read landmark data matrix from fcsv file
#'
#' Read landmark data for a single specimen from an fcsv file obtained from SlicerMorph. 
#' 
#' This function merely extracts x, y, z coordinates from an fcsv file.  
#'
#' @param file The name of a *.fcsv file containing three-dimensional landmark data to be read in
#' @keywords IO
#' @export
#' @author Murat Maga, Michael Collyer and Dean Adams
#' @return Function returns a p x 3 matrix of x, y, z coordinates for p landmarks.
#' 
#' Note: to read in multiple files the following is useful:
#' 
#' filelist <- list.files(path = "PATH TO FOLDER with FILES", pattern = "*.fcsv", full.names = TRUE)
#'
#' mydata <- simplify2array(lapply(filelist, readland.fcsv))

readland.fcsv = function (file = NULL)  {
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
