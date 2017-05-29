#'  Read landmark data from Morphologika file
#'
#'  Read Morphologika file (*.txt) to obtain landmark coordinates and specimen information
#'
#'  This function reads a *.txt file in the Morphologika format containing two- or three-dimensional 
#'  landmark coordinates. Morphologika files are text files in one of the standard formats for 
#'  geometric morphometrics (see O'Higgins and Jones 1998,2006). 
#' 
#'  If the headers "[labels]", "[labelvalues]" and "[groups]" are present, then a data matrix containing all 
#'  individual specimen information is returned.
#'  If the header "[wireframe]" is present, then a matrix of the landmark addresses for the wireframe is
#'  returned (see \code{\link{plotRefToTarget}} option 'links')
#' 
#' @param file A Morphologika *.txt file containing two- or three-dimensional landmark data. 
#' @export
#' @keywords IO
#' @author Emma Sherratt & Erik Otarola-Castillo
#' @return Function returns a (p x k x n) array of the coordinate data. If other optional headers are present in
#' the file (e.g. "[labels]" or "[wireframe]") function returns a list containing the "coords" array, 
#' and data matrix of "labels" and or "wireframe".
#' @references O'Higgins P and Jones N (1998) Facial growth in Cercocebus torquatus: An application of three 
#' dimensional geometric morphometric techniques to the study of morphological  variation. Journal of Anatomy. 
#' 193: 251-272 
#' @references O'Higgins P and Jones N (2006) Tools for statistical shape analysis. Hull York Medical School.   
read.morphologika<-function(file){
  mfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\"",quiet=TRUE)
  if(length(grep("[",mfile, fixed=T))<1){
    stop("File does not appear to be a Morphologika format text file.")}
  n <- as.numeric(mfile[grep("individuals",mfile,ignore.case=T) + 1])
  p <- as.numeric(mfile[grep("landmarks",mfile,ignore.case=T) + 1])
  k <- as.numeric(mfile[grep("dimensions",mfile,ignore.case=T) + 1])
  labvalmat <- wiref <- names <- NULL
  rawdat <- mfile[grep("rawpoints",mfile,ignore.case=T) + 1:(n*p+n)]
  rawdat <- rawdat[-grep("'",rawdat)]
  landdata <- matrix(as.numeric(unlist(strsplit(rawdat,"\\s+"))),ncol = k, byrow = T)
  if (sum(which(is.na(landdata) == TRUE)) > 0) { print("NOTE.  Missing data identified.") }
  coords<-arrayspecs(landdata,p,k)
  names <- mfile[grep("names",mfile,ignore.case=T) + 1:n]
  if(!is.null(names)) dimnames(coords)[[3]]<-names
  if(length(grep("label",mfile,ignore.case=T))>0) {
    tmp <- unlist(strsplit(mfile[grep("labels",mfile,ignore.case=T) + 1],"\\s+"))
    labvals <- unlist(strsplit(mfile[grep("labelvalues",mfile,ignore.case=T) + 1:n],"\\s+"))
    labvalmat <- matrix(labvals, ncol=length(tmp), byrow=T, dimnames=list(names,tmp)) }
  if(length(grep("groups",mfile,ignore.case=T))>0) {
    tmp <- matrix(unlist(strsplit(mfile[grep("groups",mfile,ignore.case=T) + 1],"\\s+")), ncol=2, byrow=T)
    gpval <- NULL
    for(i in 1:nrow(tmp)){ gpval <- c(gpval, rep(tmp[i,1], tmp[i,2]))}
    labvalmat <- cbind(labvalmat, groups=gpval) }
  if(length(grep("wireframe",mfile,ignore.case=T))>0) {
    strtwf <- grep("wireframe",mfile, ignore.case=T)
    endwf <-strtwf + grep("[",mfile[strtwf+1:length(mfile)], fixed=T)[1]
      if(is.na(endwf)){ endwf <- length(mfile)}
    wiref <- matrix(as.numeric(unlist(strsplit(mfile[(strtwf+1):(endwf-1)],"\\s+"))),ncol=2, byrow=T) }  
  if(!is.null(wiref) && is.null(labvalmat)) return(list(coords = coords, wireframe = wiref)) 
  if(is.null(wiref) && !is.null(labvalmat)) return(list(coords = coords, labels = labvalmat)) 
  if(!is.null(wiref) && !is.null(labvalmat)) return(list(coords = coords, labels = labvalmat, wireframe = wiref)) 
  else return(coords=coords)
}  