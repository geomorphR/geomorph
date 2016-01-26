#' Read landmark data matrix from nts file
#'
#' Read *.nts file to obtain landmark coordinates for a set of specimens
#'
#' Function reads a *.nts file containing a matrix of two- or three-dimensional landmark coordinates. 
#'   NTS files are text files in one of the standard formats for geometric morphometrics (see Rohlf 2012). 
#'   The parameter line contains 5 or 6 elements, and must begin with a "1" to designate a rectangular 
#'   matrix. The second and third values designate how many specimens (n) and how many total variables 
#'   (p x k) are in the data matrix. The fourth value is a "0" if the data matrix is complete and a "1" 
#'   if there are missing values. If missing values are present, the '1' is followed by the arbitrary 
#'   numeric code used to represent missing values (e.g., -999). These values will be replaced with "NA" 
#'   in the output array. Subsequent analyses requires a full complement of data, see \code{\link{estimate.missing}}. 
#'   The final value of the parameter line denotes the dimensionality of the landmarks
#'   (2,3) and begins with "DIM=". If specimen and variable labels are included, these are designated placing 
#'   an "L" immediately following the specimen or variable values in the parameter file. The labels then 
#'   precede the data matrix.
#'   
#'   Missing data may also be represented by designating them using 'NA'. In
#'   this case, the standard NTSYS header is used with no numeric designation for missing data (i.e. the fourth value is '0').
#'   The positions of missing landmarks may then be estimated using estimate.missing.

#'
#' Function is for *.nts file containing landmark coordinates for multiple specimens. Note that *.dta files in the 
#' nts format written by Landmark Editor \url{http://graphics.idav.ucdavis.edu/research/projects/EvoMorph},
#' and *.nts files written by Stratovan Checkpoint \url{http://www.stratovan.com/} have incorrect 
#' header notation; every header is 1 n p-x-k 1 9999 Dim=3, rather than 1 n p-x-k 0 Dim=3, which denotes
#' that missing data is in the file even when it is not.
#'
#' @param file A *.nts file containing two- or three-dimensional landmark data
#' @keywords IO
#' @export
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the image names in the *.nts file. 
#' @references  Rohlf, F. J. 2012 NTSYSpc: Numerical taxonomy and multivariate analysis system. Version 
#'   2.2. Exeter Software, New York.
readland.nts<-function(file){    	
  ntsfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\"",quiet=TRUE)
  comment <- grep("\'", ntsfile)
  if (length(comment) != 0){
    ntsfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\'",quiet=TRUE)
  }
  header<-unlist(strsplit(ntsfile[1],"\\s+"))
  if(header[1]!=1){
    stop("NTS file not a rectangular matrix. First value in parameter line must be '1'.") }
  header<-casefold(header,upper=TRUE)
  dimval<-unlist(grep("DIM=",header))
  if(length(dimval)==0){
    stop("Header does not contain 'DIM=' designator.") }  
  labval<-unlist(grep("L",header))
  r.lab<-ifelse(is.element("2",labval)==TRUE,T,F)
  c.lab<-ifelse(is.element("3",labval)==TRUE,T,F)
  header<-sub("L","",header)
  header<-as.numeric(sub("DIM=","",header))
  missdata<-ifelse(header[4]!=0,T,F)
  if(missdata==TRUE){missval<-ifelse(dimval==6,header[5],header[6]) } 
  n<-header[2];k<-header[dimval];p<-header[3]/k;   
  tmp<-unlist(strsplit(ntsfile[-1],"\\s+"))
  speclab<-NULL; 
  if(r.lab==TRUE){
    #     speclab<-ntsfile[2:(1+n)]
    #     tmp <- tmp[c((length(tmp)-(p*k)+1):length(tmp))]
    speclab<-tmp[1:n]
    tmp<-tmp[-(1:length(speclab))]  
  }
  if(c.lab==TRUE){ tmp<-tmp[-(1:(p*k))] }
  if(missdata==TRUE){tmp[grep(missval,as.integer(tmp))] <- NA}
  options(warn=-1)
  landdata<-matrix(as.numeric(tmp),ncol=k,byrow=TRUE)
  if(sum(which(is.na(landdata)==TRUE))>0){print("NOTE.  Missing data identified.")}
  coords <- aperm(array(t(landdata), c(k,p,n)), c(2,1,3))
  dimnames(coords)[[3]]<-as.list(speclab)
  return(coords=coords)
}