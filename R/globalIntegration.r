#' Quantify global integration relative to self-similarity
#'
#' Function quantifies the overall level of morphological integration for a set of Procrustes shape variables
#'
#' The function quantifies the overall level of morphological integration for a set of 
#' Procrustes shape coordinates. It is assumed that the landmarks have previously been 
#' aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. 
#' Based on the set of aligned specimens, the function estimates the set of bending energies at various
#' spatial scales, and plots the log of the variance of the partial warps versus the log of their
#' corresponding bending energies (Bookstein 2015). The slope of a regression of these data provides information
#' regarding the degree of overall morphological integration (or lack thereof).
#' 
#' A slope of negative one corresponds to self-similarity, implying that patterns of shape variation are
#' similar across spatial scales. Steeper slopes (i.e., those more extreme than -1.0) correspond to data that are globally
#' integrated, while shallower slopes (between -1 and 0) correspond to data that are 'disintegrated (see Bookstein 2015). Isotropic data
#' will have an expected slope of zero. 
#'  
#' @param A 3D array (p1 x k x n) containing Procrustes shape variables 
#' @param ShowPlot A logical value indicating whether or not the plot should be returned
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @references  Bookstein, F. L. 2015. Integration, disintegration, and self-similarity: 
#' Characterizing the scales of shape variation in landmark data. Evol. Biol.42(4): 395-426.
#' 
#' @examples
#' \dontrun{
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)    #GPA-alignment    
#'
#' globalIntegration(Y.gpa$coords)
#' }

globalIntegration<-function(A,ShowPlot=TRUE){
  ref<-mshape(A)
  p<-dim(ref)[1]; k<-dim(ref)[2]  
  Pdist<-as.matrix(dist(ref))
  if(k==2){P<-Pdist^2*log(Pdist^2)}; if(k==3){P<- -Pdist}
  P[which(is.na(P))]<-0
  Q<-cbind(1, ref)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  Linv<-solve(L)
  L.be<-Linv[1:p,1:p]
  eig.L<-eigen(L.be, symmetric = TRUE)
  if(any(Im(eig.L$values)==0)) eig.L$values<-Re(eig.L$values)
  if(any(Im(eig.L$vectors)==0)) eig.L$vectors<-Re(eig.L$vectors)
  BE<-eig.L$values[1:(p-3)]; if(k==3){BE<-BE[1:(p-4)]}
  lambda <- zapsmall(eig.L$values)
  if(any(lambda == 0)){BE = lambda[lambda > 0]}
  Emat<- eig.L$vectors[,1:(length(BE))]
  BEval<-log(BE)
  Wmat<-NULL
  for (i in 1:dim(A)[[3]]){ 
    Wvec<-matrix(t(t(A[,,i])%*%Emat),nrow=1)
    Wmat<-rbind(Wmat,Wvec)
  }	
  Wvar<-diag(var(Wmat))
  Tmpvar<-matrix(Wvar,nrow=k,byrow=T); PWvar<-log(apply(Tmpvar,2,sum))
  start<-which.min(BEval)
  slope<-coef(lm(PWvar~BEval))[2]
  eq<-bquote(ObservedSlope [black] == .(slope))
  if(ShowPlot==TRUE){ 
  plot(BEval,PWvar,asp=1,main=eq)
  abline(lm(PWvar~BEval),lwd=2,col="black")
  lines(c(BEval[start],BEval[start]+10),c(PWvar[start],PWvar[start]-10),lty=3,lwd=2,col="red")}
  return(slope)
}
