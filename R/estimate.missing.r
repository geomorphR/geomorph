#' Estimate locations of missing landmarks
#'
#' A function for estimating the locations of missing landmarks 
#' 
#' The function estimates the locations of missing landmarks for incomplete specimens  in a set of landmark
#' configurations, where missing landmarks in the incomplete specimens are designated by NA in place 
#' of the x,y,z coordinates.  Two distinct approaches are implemented.
#' 
#' The first approach (method="TPS") uses the thin-plate spline to interpolate landmarks on a reference specimen to estimate 
#' the locations of missing landmarks on a target specimen. Here, a reference specimen is obtained from 
#' the set of specimens for which all landmarks are present, Next, each incomplete specimen is aligned to 
#' the reference using the set of landmarks common to both. Finally, the thin-plate spline is used 
#' to estimate the locations of the missing landmarks in the target specimen (Gunz et al. 2009).
#' 
#' The second approach (method="Reg") is multivariate regression. Here each landmark with missing values is
#' regressed on all other landmarks for the set of complete specimens, and the missing landmark values are
#' then predicted by this linear regression model. Because the number of variables can exceed the number of
#' specimens, the regression is implemented on scores along the first set of PLS axes for the complete and 
#' incomplete blocks of landmarks (see Gunz et al. 2009).
#'  
#'  One can also exploit bilateral symmetry to estimate the locations of missing landmarks. Several
#'   possibilities exist for implementing this approach (see Gunz et al. 2009).  Example R code for one 
#'   implementation is found in Claude (2008).
#'  
#' NOTE: Because all geometric morphometric analyses and plotting functions implemented in geomorph 
#' require a full complement of landmark coordinates, the alternative to estimating the missing 
#' landmark coordinates is to proceed with subsequent analyses EXCLUDING
#' specimens with missing values. To do this, see functions \code{\link[stats]{complete.cases}} (use: mydata[complete.cases(mydata),])
#' or \code{\link[stats]{na.omit}} (use: newdata <- na.omit(mydata)) to make a dataset of only the complete specimens.
#' These functions require the dataset to be a matrix in the form of a 2d array (see \code{\link{two.d.array}}).
#' 
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @param method Method for estimating missing landmark locations
#' @author Dean Adams
#' @keywords utilities
#' @return Function returns an array (p x k x n) of the same dimensions as input A, including coordinates for the target specimens 
#' (the original landmarks plus the estimated coordinates for the missing landmarks). These data need to be Procrustes Superimposed prior to analysis (seesee \code{\link{gpagen}}).
#' @export
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York.
#' @references  Bookstein, F. L., K. Schafer, H. Prossinger, H. Seidler, M. Fieder, G. Stringer, G. W. Weber, 
#' J.-L. Arsuaga, D. E. Slice, F. J. Rohlf, W. Recheis, A. J. Mariam, and L. F. Marcus. 1999. Comparing 
#' frontal cranial profiles in archaic and modern Homo by morphometric analysis. Anat. Rec. (New Anat.) 257:217-224.
#' @references Gunz, P., P. Mitteroecker, S. Neubauer, G. W. Weber, and F. L. Bookstein. 2009. Principles for 
#' the virtual reconstruction of hominin crania. J. Hum. Evol. 57:48-62.
#' @examples
#' data(plethodon)
#' plethland<-plethodon$land
#'   plethland[3,,2]<-plethland[8,,2]<-NA  #create missing landmarks
#'   plethland[3,,5]<-plethland[8,,5]<-plethland[9,,5]<-NA  
#'   plethland[3,,10]<-NA  
#'   
#' estimate.missing(plethland,method="TPS")
#' estimate.missing(plethland,method="Reg")
estimate.missing<-function(A,method=c("TPS","Reg")){ 
  if(any(is.na(A))==FALSE)  {stop("No missing data.")}
  method <- match.arg(method)
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(method=="TPS"){
    A2<-A
    spec.NA<-which(rowSums(is.na(two.d.array(A)))>0)
    Y.gpa<-gpagen(A[,,-spec.NA],PrinAxes=FALSE)
    p<-dim(Y.gpa$coords)[1]; k<-dim(Y.gpa$coords)[2]; n<-dim(Y.gpa$coords)[3]    
    ref<-mshape(arrayspecs(two.d.array(Y.gpa$coords)*Y.gpa$Csize,p,k))
    for (i in 1:length(spec.NA)){
      missing<-which(is.na(A2[,1,spec.NA[i]])== T)
      tmp<-tps2d3d(ref,ref[-missing,],A2[-missing,,spec.NA[i]])
      A2[,,spec.NA[i]]<-tmp
    }
  }
  if(method=="Reg"){
    spec.NA<-which(rowSums(is.na(two.d.array(A)))>0)    
    land.NA<-which(colSums(is.na(two.d.array(A)))>0)  
    p<-dim(A)[1]; k<-dim(A)[2];  
    A2<-A
    complete<-A[,,-spec.NA]
    incomplete<-A[,,spec.NA]
    Y.gpa<-gpagen(complete, PrinAxes=FALSE)
    ref<-mshape(arrayspecs(two.d.array(Y.gpa$coords)*Y.gpa$Csize,p,k))
    complete<-arrayspecs(two.d.array(Y.gpa$coords)*Y.gpa$Csize,p,k)
    if(length(dim(incomplete))>2){
      for (i in 1:dim(incomplete)[3]){
        missing<-which(is.na(incomplete[,1,i])== TRUE)
        lndmk<-which(is.na(incomplete[,1,i])!= TRUE)
        tmp <- apply.pPsup(center(ref[-missing,]), list(center(incomplete[-missing,,i])))[[1]]
        incomplete[lndmk,,i]<-tmp
      }      
    }
    if(length(dim(incomplete))==2){
      missing<-which(is.na(incomplete[,1])== T)
      lndmk<-which(is.na(incomplete[,1])!= T)
      tmp <- apply.pPsup(center(ref[-missing,]), list(center(incomplete[-missing,])))[[1]]
      incomplete[lndmk,]<-tmp
    }
    A2[,,-spec.NA]<-complete
    A2[,,spec.NA]<-incomplete
    A.2d<-two.d.array(A2)    
    for (i in 1:length(spec.NA)){
      missing<-which(is.na(A.2d[spec.NA[i],])== T)
      x<-A.2d[-spec.NA,-missing]
      y<-A.2d[-spec.NA,missing]
      XY.vcv<-cov(cbind(x,y))      
      S12<-XY.vcv[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]
      pls<-svd(S12)
      U<-pls$u; V<-pls$v
      XScores<-x%*%U; YScores<-y%*%V
      beta<-coef(lm(YScores[,1]~XScores[,1]))
      miss.xsc<-c(1,A.2d[spec.NA[i],-missing]%*%U[,1])
      miss.ysc<-c(miss.xsc%*%beta,(rep(0,(ncol(y)-1))))
      pred.val<-miss.ysc%*%t(V)
      for (j in 1:length(missing)){
        A.2d[spec.NA[i],missing[j]]<-pred.val[j]    
      }
    }
    A2<-arrayspecs(A.2d,dim(A)[1],dim(A)[2])
  }
  return(A2)
}
