#' Write landmark data to tps file
#'
#' Write *.tps file from obtain landmark coordinates in a 3-dimensional array 
#'
#' This function writes a *.tps file from a 3-dimensional array (p x k x n) 
#'  of landmark coordinates. 
#'
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param file Name of the *.tps file to be created
#' @param scale An optional vector containing the length of the scale for each specimen
#' @param specID A logical value stating whether specimen ID names should be saved to line ID=
#' @export
#' @keywords IO
#' @author Dean Adams
writeland.tps<-function(A, file, scale = NULL, specID = TRUE){
  n<-dim(A)[3]
  k<-dim(A)[2]
  p<-dim(A)[1]
  lmline<-ifelse(k==2,paste("LM=",p,sep=""), paste("LM3=",p,sep=""))  
  file.create(file, showWarnings=TRUE)
  if(!is.null(scale)){
    scaleline<-paste("SCALE", "=", scale, sep="")
  }
  for(i in 1:n){
    write(lmline,file,append = TRUE)
    write.table(A[,,i],file,col.names = FALSE, row.names = FALSE,append=TRUE)
    if(!is.null(scale)){
      if(length(scaleline) == 1){write(scaleline,file,append=TRUE)}
      if(length(scaleline) > 1){write(scaleline[i],file,append=TRUE)}
    }
    if(specID==TRUE){
      if(is.null(dimnames(A)[[3]])) dimnames(A)[[3]] <- c(1:dim(A)[3])
      idline<-paste("ID=",dimnames(A)[[3]][i],sep="")
      write(idline,file,append = TRUE)  
    }
    write("",file,append = TRUE)
  }
}
