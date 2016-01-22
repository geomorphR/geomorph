#' Read landmark data from tps file
#'
#' Read *.tps file to obtain landmark coordinates
#'
#' This function reads a *.tps file containing two- or three-dimensional landmark coordinates. 
#'   Tps files are text files in one of the standard formats for geometric morphometrics (see Rohlf 2010). 
#'   Two-dimensional landmarks coordinates are designated by the identifier "LM=", while three-dimensional 
#'   data are designated by "LM3=". Landmark coordinates are multiplied by their scale factor if this is 
#'   provided for all specimens. If one or more specimens are missing the scale factor, landmarks are treated 
#'   in their original units.  
#'   
#'   Missing data may be present in the file. In this case, they must be designated by 'NA'. The 
#'   positions of missing landmarks may then be estimated using estimate.missing.
#' 
#' The user may specify whether specimen names are to be extracted from the 'ID=' field or 'IMAGE=' field 
#' and included in the resulting 3D array. 
#' e.g., for 'ID=' use (file, specID = "ID") and for 'IMAGE=' use (file, specID = "imageID").
#' The default is specID="None".
#' 
#' If there are curves defined in the file (i.e., CURVES= fields), the option 'readcurves' should be used.
#' When readcurves = TRUE, the coordinate data for the curves will be returned as semilandmarks and will be appended to
#' the fixed landmark data. Then the user needs to use \code{\link{define.sliders}} or \code{\link{define.sliders}}
#' to create a matrix designating how the curve points will slide (used by 'curves=' in \code{\link{gpagen}}).
#' When readcurves = FALSE, only the landmark data are returned.
#' 
#' NOTE: At present, all other information that can be contained in tps files (comments, variables, radii, etc.)
#'   is ignored. 
#'
#' @param file A *.tps file containing two- or three-dimensional landmark data
#' @param specID a character specifying whether to extract the specimen ID names from the ID or IMAGE lines (default is "None").
#' @param readcurves A logical value stating whether CURVES= field and associated coordinate data will be read as semilandmarks (TRUE)
#' or ignored (FALSE).
#' @param warnmsg A logical value stating whether warnings should be printed
#' @export
#' @keywords IO
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the image names in the *.tps file. 
#'   
#' @references  Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology 
#'   and Evolution, State University of New York at Stony Brook, Stony Brook, NY.

readland.tps <- function (file, specID = c("None", "ID", "imageID"), 
                          readcurves = FALSE, warnmsg = T) 
{
  ignore.case = TRUE
  specID <- match.arg(specID)
  tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
  lmdata <- grep("LM=", tpsfile, ignore.case)
  if (length(lmdata !=0)) {
    nland <- as.numeric(sub("LM=", "", tpsfile[lmdata], ignore.case))
    k <- 2
  }
  if (length(lmdata) == 0) {
    lmdata <- grep("LM3=", tpsfile, ignore.case)
    nland <- as.numeric(sub("LM3=", "", tpsfile[lmdata], ignore.case))
    k <- 3
  }
  n <- nspecs <- length(lmdata)
  if (max(nland) - min(nland) != 0) {
    stop("Number of landmarks not the same for all specimens.")
  }
  p <- nland[1]
  imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE", 
                                                       tpsfile, ignore.case)], ignore.case))
  if (is.null(imscale)) {
    imscale = array(1, nspecs)
  }
  if (warnmsg == T) {
    if (length(imscale) != nspecs) {
      print("Not all specimens have scale. Assuming landmarks have been previously scaled.")
    }
  }
  if (length(imscale) != nspecs) {
    imscale = array(1, nspecs)
  }
  crvs <- grep("CURVES=", tpsfile, ignore.case)
  if(length(crvs)>0){
    if (readcurves == TRUE && length(crvs) == 0){ stop("No CURVES= field present in file") } 
    ncurve <- as.numeric(sub("CURVES=", "", tpsfile[crvs], ignore.case))
    ncurvepts <- as.numeric(sub("POINTS=", "", tpsfile[grep("POINTS=", tpsfile, ignore.case)], ignore.case))
      if (max(ncurve) - min(ncurve) != 0) {
        stop("Number of curves not the same for all specimens.") }
      if (warnmsg == T && readcurves==T) {print(paste("Landmarks 1:", p, " are fixed landmarks.", sep=""))
                         print(paste("Landmarks ", p+1, ":", p+sum(ncurvepts[1:ncurve[1]]), " are semilandmarks.", sep=""))}
      p <- nland[1] + sum(ncurvepts[1:ncurve[1]]) 
  }    
  tmp <- tpsfile[-(grep("=", tpsfile))]
  options(warn = -1)
  tmp <- matrix(as.numeric(unlist(strsplit(tmp,"\\s+"))),ncol = k, byrow = T)
 
  if (warnmsg == T) {
    if (sum(which(is.na(tmp) == TRUE)) > 0) {
      print("NOTE.  Missing data identified.")
    }
  }
  coords <- aperm(array(t(tmp), c(k, p, n)), c(2, 1, 3))
  imscale <- aperm(array(rep(imscale, p * k), c(n, k, p)), 
                   c(3, 2, 1))
  coords <- coords * imscale
  if (readcurves==F){coords<-coords[1:nland,,] 
      if(n==1) coords <- array(coords, c(nland,k,n))}
  if (specID == "None") {
      if (warnmsg == T) {print("No Specimen names extracted")
    }
  }
  if (specID == "imageID") {
    imageID <- (sub("IMAGE=", "", tpsfile[grep("IMAGE", tpsfile, ignore.case)], 
                    ignore.case))
    if (length(imageID) != 0) {
      imageID <- sub(".jpg", "", imageID, ignore.case)
      imageID <- sub(".tif", "", imageID, ignore.case)
      imageID <- sub(".bmp", "", imageID, ignore.case)
      imageID <- sub(".tiff", "", imageID, ignore.case)
      imageID <- sub(".jpeg", "", imageID, ignore.case)
      imageID <- sub(".jpe", "", imageID, ignore.case)
      dimnames(coords)[[3]] <- as.list(imageID)
      if (warnmsg == T) {
        print("Specimen names extracted from line IMAGE=")
      }
    }
    if (length(imageID) == 0) {
      if (warnmsg == T) {
        print("No name given under IMAGE=. Specimen names not extracted")
      }
    } 
  }
  if (specID == "ID") {
    ID <- sub("ID=", "", tpsfile[grep("ID", tpsfile, ignore.case)], ignore.case)
    if (length(ID) == 0) {
      if(warnmsg ==T){
        print("No name given under ID=. Specimen names not extracted")
        }
      }
    if (length(ID) != 0) {
      dimnames(coords)[[3]] <- as.list(ID)
      if (warnmsg == T) {
        print("Specimen names extracted from line ID=")
      }
    }
  }
return(coords = coords)                    
}
