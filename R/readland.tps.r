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
  tpsf <- scanTPS(file)
  n <- length(tpsf)
  specID <- match.arg(specID)
  if(specID == "ID") id <- sapply(1:n, function(j) tpsf[[j]]$id) else
    if(specID == "imageID") id <- sapply(1:n, function(j) tpsf[[j]]$image) else {
      if(warnmsg) cat("\nNo specID provided; specimens will be numbered 1, 2, 3 ...\n")
      id <- 1:n
    }
  if(warnmsg){
    if(specID == "ID" && length(id[[1]]) == 0) {
      cat("\nWarning: specID = 'ID' did not produce reliable ID information;")
      cat("\nspecimens will be numbered 1, 2, 3 ...\n")
      id <- 1:n
    }
    if(specID == "imageID" && length(id[[1]]) == 0) {
      cat("\nWarning: specID = 'imageID' did not produce reliable ID information;")
      cat("\nspecimens will be numbered 1, 2, 3 ...\n")
      id <- 1:n
    }
    
    pcv.check <- sapply(1:n, function(j) tpsf[[j]]$pcv)
    pcv.unique <- unique(pcv.check)
    if(length(pcv.unique) == 1 && pcv.unique == 0) {
      cat("\nNo curves detected; all points appear to be fixed landmarks.\n")
    } else if(length(pcv.unique) == 1 && pcv.unique > 0) {
      cat("\n", pcv.unique, "curve points detected per specimen and are appended to fixed landmarks.\n")
    } else if(length(pcv.unique) > 1){
      cat("\nCurve points detected but numbers vary among specimens.\n")
      cat("\nCurve point frequencies:\n")
      print(table(factor(pcv.check)))
    }
  }
  
  kcheck <- sapply(1:n, function(j) length(tpsf[[j]]$k))
  k.error <- which(kcheck > 1)
  if(length(k.error) == 0) k.error <- NULL else
    cat("\nWarning: improper landmark number or formatting appear for specimen(s):", k.error,"\n")
  
  pcheck <- sapply(1:n, function(j) tpsf[[j]]$p)
  if(all(pcheck==0)) {
    stop(paste("File", file, "does not contain landmark coordinates", sep = " "))
  }
  
  p.unique <- unique(pcheck)
  if(length(p.unique) > 1) {
    names(pcheck) <- id
    cat(paste("Number of landmarks per specimen in", file, sep = " "))
    print(as.matrix(pcheck))
    stop("\nDifferent numbers of landmarks among specimens")
  } else p.error <- NULL
  
  scale.list <- unlist(lapply(1:n, function(j) tpsf[[j]]$scale))
  if(length(scale.list) != n && warnmsg) {
    cat("\nWarning: not all specimens have scale adjustment (perhaps because they are already scaled);")
    cat("\nno rescaling will be performed in these cases\n")
  }
  
  if(!readcurves) {
    lmo <- lapply(1:n, function(j) {
      x <- tpsf[[j]]
      l <- x$lm
      k <- max(x$k)
      p <- x$plm
      lm <- matrix(NA, p, k)
      for(i in 1:p) {
        pts <- unlist(l[[i]])
        kk <- length(pts)
        if(kk > 0) lm[i,1:kk] <- pts
      }
      
      if(length(x$scale) == 0) x$scale = 1
      lm*x$scale
    })
  } else {
    lmo <- lapply(1:n, function(j) {
      x <- tpsf[[j]]
      l <- c(x$lm, x$curve.lm)
      k <- max(x$k)
      p <- x$plm + x$pcv
      lm <- matrix(NA, p, k)
      for(i in 1:p) {
        pts <- unlist(l[[i]])
        kk <- length(pts)
        if(kk > 0) lm[i,1:kk] <- pts
      }
      if(length(x$scale) == 0) x$scale = 1
      lm*x$scale
    })
  }
  
  lmo <- try(simplify2array(lmo), silent = TRUE)
  lmo <- two.d.array(lmo)
  
  if(any(na.omit(lmo) < 0)){
    if(negNA == TRUE){
      lmo[which(lmo < 0)] <- NA
    } else {
      cat("\nNegative landmark coordinates have been identified and imported as such.") 
      cat("\nIf you want to treat them as NAs please set negNA = TRUE")
    }
  }
  lmo <- arrayspecs(lmo, p.unique, ncol(lmo)/p.unique)
  
  if(!is.null(p.error) && warnmsg) {
    target <- as.numeric(names(sort(p.error, decreasing = TRUE))[1])
    p.error <- pcheck != target
    badspec <- unique(c(id[k.error], id[p.error]))
    cat("\nThere was a problem because of imbalanced data!")
    cat("\nBased on the specID argument used,")
    cat("\ncheck the following specimens for landmark issues:", badspec)
  } 
  if(!is.null(k.error)){
    cat("\nThere appears to be missing or superfluous data.\n")
    cat("Check these specimens for mistakes:", id[k.error], "\n")
  }
  if(n==1){if(is.array(lmo)) dimnames(lmo)[[3]] <- list(id)} else {if(is.array(lmo)) dimnames(lmo)[[3]] <- id} #added check for N=1 in file
  if(is.list(lmo)){
    cat("\n\nNote that the landmark array may not be properly formatted,")
    cat("\nin which case a list of landmarks by specimen is available for inspection.\n")
    names(lmo) <- id
  }
  
  if(warnmsg){
    if(length(pcv.unique) > 1 || length(pcv.unique) > 1){
      cat("\n\nA break down of fixed and curve points per specimen:\n")
      sp.list <- sapply(1:n, function(j) {
        x <- tpsf[[j]]
        c(x$plm, x$pcv)
      })
      rownames(sp.list) <- c("nFixedPts", "nCurvePts")
      colnames(sp.list) <- id
      print(sp.list)
    }
  }
  
  invisible(lmo)
}
