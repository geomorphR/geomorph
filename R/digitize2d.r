#' Digitize 2D landmarks on .jpg files
#'
#' An interactive function to digitize two-dimensional(2D) landmarks from .jpg files.
#'
#' This function may be used for digitizing 2D landmarks from jpeg images (.jpg). The user provides 
#' a list of image names, the number of landmarks to be digitized, and the name of an output
#' TPS file.  An option is included to allow the user to digitize a scale on each image to convert 
#' the landmark coordinates from pixels into meaningful units. Landmarks to be digitized can include 
#' both fixed landmarks and semi-landmarks, the latter of which are to be designated as "sliders" 
#' for subsequent analysis (see the function \code{\link{define.sliders}}).
#' 
#' \subsection{The Digitizing Session}{
#' Digitizing landmarks from 2D photos requires that a scale bar is placed in the image in order to scale the coordinate data. The 'scale' option requires: a single number (e.g. 10) which means that the scale to be measured in all images is a 10mm scale bar; OR a vector the same length as the filelist containing a number for the scale of each image. If scale=NULL, then the digitized coordinates will not be scaled. This option is NOT recommended.
#' 
#' Users may digitize all specimens in one session, or may return at a later time to complete digitizing. In the latter case, the user provides the same filelist and TPS file and the function will determine where the user left off. 
#' 
#' If specimens have missing landmarks, these can be incorporated during the digitizing process using the 'a' option as described below (a=absent).  
#' }
#' 
#' 
#' \subsection{Specimen Digitizing}{ 
#' Digitizing landmarks involves landmark selection using a mouse in the plot window, 
#' using the LEFT mouse button (or regular button for Mac users):
#' \enumerate{
#'  \item Digitize the scale bar (if requested) by selecting the two end points. Use a single click for start and end points. The
#'   user is asked whether the system should keep or discard the digitized scale bar. 
#'  \item Digitize each landmark with single click and the landmark is shown in red. 
#' }
#' If verbose = TRUE, digitizing is interactive between landmark selection using a mouse and the R console. 
#' Once a landmark is selected, the user is asked if the system should keep or discard the 
#' selection (y/n/a). If "y", the user is asked to continue to select the next landmark. If "n", the user is 
#' asked to select it again.
#' 
#'  To digitize a missing landmark, simply click on any location in the image. Then, when 
#'  prompted to keep selection, choose 'a' (for absent).  Missing landmarks can only be included during the digitizing process when verbose=TRUE. 
#'  
#' If verbose = FALSE the digitizing of landmarks is continuous and uninterrupted. Here the user will not be prompted to approve each landmark selection. 
#'  
#'   At the end of digitizing, the landmark coordinates are written to a TPS file. By default, the 
#'   x,y values are unscaled if a vector of scales is included, and the scale is returned on line 
#'   SCALE= after each specimen x,y data.  Optionally, one may have the coordinates pre-multiplied by
#'   scale by using the option MultScale=TRUE.
#'   }
#' 
#' @param filelist A list of names of jpeg images to be digitized. 
#' @param nlandmarks Number of landmarks to be digitized.
#' @param scale An vector containing the length of the scale to be placed on each image.
#' @param tpsfile The name of a TPS file to be created or read
#' @param MultScale A logical option indicating if the coordinates should be pre-multiplied by scale
#' @param verbose logical. User decides whether to digitize in verbose or silent format (see details), default is verbose
#' @return Function returns a tps file containing the digitized landmark coordinates.
#' @keywords digitizing
#' @export
#' @author Dean Adams, Erik Otarola-Castillo and Emma Sherratt
digitize2d <- function (filelist, nlandmarks, scale=NULL, tpsfile, MultScale=FALSE,verbose = TRUE) 
{
  flist <- dir()
  if (sum(which(flist == tpsfile)) == 0) {
    newdata <- array(0, c(nlandmarks, 2, length(filelist)))
    dimnames(newdata)[[3]] <- as.list(filelist)
    writeland.tps(newdata, tpsfile)
  }
  olddat <- readland.tps2(file=tpsfile, specID = "ID")
  
  newdata<-olddat$coords
  inscale<-olddat$scale
  names <- dimnames(newdata)[[3]]
  if (dim(newdata)[3] != length(filelist)) {
    stop("Filelist not the same length as input TPS file.")
  }
  if (length(na.omit(match(names, filelist))) != length(filelist)) {
    stop("Filelist does not contain the same specimens as TPS file.")
  }
  if(length(scale) != length(filelist)){
    if(length(scale)==1) { cat("Only 1 scale measure provided. Will use scale =", scale, " for all specimens.\n")
      scale = rep(scale, length(filelist))
    }
  }
  digitized <- apply(two.d.array(newdata), 1, sum)
  digstart <- min(which(digitized == 0))
  if(is.infinite(digstart)){
    stop("It appears that all specimens have already been digitized.")
  }
  scalebar=vector(length=length(filelist))
  if(digstart>1){scalebar[1:(digstart-1)]<-inscale[1:(digstart-1)]}  
  for (i in digstart:length(filelist)) {
    cat(paste("Digitizing specimen", i, "in filelist"), "\n")
    spec.name <- unlist(strsplit(basename(filelist[i]), "\\."))[1]
    specimen <- readJPEG(filelist[i], native = T)
    plot(seq(0, dim(specimen)[2], length.out = 10), seq(0, 
                                                        dim(specimen)[1], length.out = 10), type = "n", xlab = "x", 
         ylab = "y", asp = 1, tck = 0, xaxt = "n", yaxt = "n")
    rasterImage(specimen, 1, 1, dim(specimen)[2], dim(specimen)[1])
    if (is.null(scale)) {
      cat("Scale not provided! Proceed with caution.\n") 
      scalebar = NULL
    }
    if (!is.null(scale)) {
      cat("Set scale =", scale[i], "\n")
      scalebar[i] <- picscale(scale[i])
    }
    selected <- matrix(NA, nrow = nlandmarks, ncol = 2)
    fix <- NULL
    if (verbose == TRUE) {
      for (ii in 1:nlandmarks) {
        cat("Select landmark ", ii, "\n")
        keep <- ans <- NULL
        fix <- locator(n = 1, type = "p", col = "black", 
                       cex = 4, pch = 21, bg = "red")
        cat(paste("Keep Landmark ", ii, "(y/n/a)?"), 
            "\n")
        ans <- readLines(n = 1)
        if (ans == "y") {
          selected[ii, 1] <- fix$x
          selected[ii, 2] <- fix$y
        }
        if (ans == "a") {
          selected[ii, 1] <- NA
          selected[ii, 2] <- NA
          points(fix, type = "n", col = "black", cex = 1, pch = 21, bg = "red")
        }
        if (ans == "n") {
          cat(paste("Select Landmark ", ii, " Again"), "\n")
        }
        while (ans == "n") {
          fix <- locator(n = 1, type = "p", col = "black", 
                         cex = 4, pch = 21, bg = "red")
          cat(paste("Keep Landmark ", ii, "(y/n)?"),  "\n")
          ans <- readLines(n = 1)
          if (ans == "y") {
            selected[ii, 1] <- fix$x
            selected[ii, 2] <- fix$y
          }
          if (ans == "n") {
            cat(paste("Select Landmark ", ii, " Again"),   "\n")
          }
        }
      }
    }
    if (verbose == FALSE) {
      cat("Select landmarks 1:", nlandmarks, "\n", sep = "")
      for (ii in 1:nlandmarks) {
        fix <- locator(n = 1, type = "p", col = "black", 
                       cex = 4, pch = 21, bg = "red")
        selected[ii, 1] <- fix$x
        selected[ii, 2] <- fix$y
      }
    }
    
    newdata[, , i] <- selected
    if(MultScale==TRUE){newdata[,,i]<-selected*scalebar[i]}
    if(is.null(scalebar)){writeland.tps(newdata, tpsfile)}
    if(!is.null(scalebar) && MultScale==FALSE){
      writeland.tps(newdata, tpsfile, scale = scalebar)}
    if(!is.null(scalebar) && MultScale==TRUE){
      writeland.tps(newdata, tpsfile)}
    
    if (i < length(filelist)) {
      cat(paste("Continue to next specimen (y/n)?"), "\n")
      ans <- readLines(n = 1)
      if (ans == "n") {
        break
      }
    }
  }
  cat(paste("All specimens from session have been digitized to",  tpsfile), "\n")
}