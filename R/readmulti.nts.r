#' Read landmark data from multiple nts files
#'
#' Read a list of names for several *.nts files to obtain landmark coordinates for a set of specimens
#'
#' This function reads a list containing the names of multiple *.nts files, where each contains the 
#'   landmark coordinates for a single specimen. For these files, the number of variables (columns) of 
#'   the data matrix will equal the number of dimensions of the landmark data (k=2 or 3). When the function 
#'   is called a dialog box is opened, from which the user may select multiple *.nts files. These are then read 
#'   and concatenated into a single matrix for all specimens. The parameter line contains 5 or 6 elements, 
#'   and must begin with a "1" to designate a rectangular matrix. The second and third values designate how 
#'   many specimens (n) and how many total variables (p x k) are in the data matrix.
#'   The fourth value is a "0" if the data matrix is complete and a "1" 
#'   if there are missing values. If missing values are present, the '1' is followed by the arbitrary 
#'   numeric code used to represent missing values (e.g., -999). These values will be replaced with "NA" 
#'   in the output array. Subsequent analyses requires a full complement of data, see \code{\link{estimate.missing}}. 
#'
#'   Missing data may be present in the file by designating them using 'NA'. In
#'   this case, the standard NTSYS header is used with no numeric designation for missing data (i.e. the fourth value is '0').
#'   The positions of missing landmarks may then be estimated using estimate.missing.

#'
#' @param filelist A list of names for the *.nts files to be read by the function. The names in the list
#'   require quotes (") and .nts/.NTS suffix.
#' @keywords IO
#' @export
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is 
#'   the number of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension 
#'   of this array contains names for each specimen, which are obtained from the original file names. 
readmulti.nts<-function(filelist){   
  n<-length(filelist)
  names<-gsub(".nts","",filelist, ignore.case=T)
  landdata<-nind<-NULL
  for (i in 1:n){
    ntsfile<-scan(file=filelist[i],what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\"",quiet=TRUE)
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
    p<-header[2];k<-header[3]
    nind<-rbind(nind,p)
    if (min(nind)!=max(nind)) {
      stop("Number of landmarks not the same in all files.") } 
    tmp<-unlist(strsplit(ntsfile[-1],"\\s+"))
    rowlab<-NULL; 
    if(r.lab==TRUE){
      rowlab<-tmp[1:p]
      tmp<-tmp[-(1:length(rowlab))]   }
    if(c.lab==TRUE){ tmp<-tmp[-(1:k)] }
    if(missdata==TRUE){tmp[grep(missval,as.integer(tmp))] <- NA}
    options(warn=-1)
    data<-matrix(as.numeric(tmp),ncol=k,byrow=TRUE)
    landdata<-rbind(landdata,data)
  }
  coords<-arrayspecs(landdata,p,k)
  if(sum(which(is.na(landdata)==TRUE))>0){print("NOTE.  Missing data identified.")}
  dimnames(coords)[[3]]<- as.list(names)
  return(coords=coords)
}