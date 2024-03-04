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
#'  Missing data may be present in the file. If data were digitized in geomorph or StereoMorph, these are 
#'  automatically identified. If, instead, data digitizing took place in a software package that records missing
#'  values as negative coordinates (e.g. tpsDig), the user needs to specify whether negative values should be
#'  transformed to 'NAs' through the argument neg.NA = TRUE. The positions of missing landmarks may 
#'  then be estimated using \code{\link{estimate.missing}}.
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
#' @param negNA A logical value indicating whether negative landmark coordinates in the tps file should be imported as missing 
#' values and coded as 'NA' (TRUE), or imported as such.
#' @param readcurves A logical value stating whether CURVES= field and associated coordinate data will be read as semilandmarks (TRUE)
#' or ignored (FALSE).
#' @param warnmsg A logical value stating whether warnings should be printed
#' @export
#' @keywords IO
#' @author Dean Adams, Emma Sherratt, Michael Collyer & Antigoni Kaliontzopoulou
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the image names in the *.tps file. 
#'   
#' @references  Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology 
#'   and Evolution, State University of New York at Stony Brook, Stony Brook, NY.

readland.tps <- function (file, specID = c("None", "ID", "imageID"), negNA = FALSE,  
                          readcurves = FALSE, warnmsg = TRUE) {
  lmi <- .readland.tps(file, specID, negNA, readcurves, warnmsg)
  tbl <- data.frame(id = names(lmi), 
                    p = sapply(lmi, NROW),
                    k = sapply(lmi, NCOL))
  rownames(tbl) <- NULL
  
  lm.check <- apply(tbl[, -1], 2, unique)
  
  if(is.list(lm.check)){
    cat("\nEither the number of landmarks (p) or the landmark dimensions (k) differ\n")
    cat("among specimens.  An array is not returned but the following table is\n")
    cat("provided so that discrepencies can be investigated.\n\n")
    
    print(tbl)
    lmo <- tbl
    
  } else {
    lmo <- simplify2array(lmi)
  }
  
  invisible(lmo)
  
}
