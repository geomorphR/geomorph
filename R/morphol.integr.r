#' Quantify morphological integration between two modules
#'
#' Function quantifies the degree of morphological integration between two modules of Procrustes-aligned 
#'   coordinates
#'
#' The function quantifies the degree of morphological integration between two modules of shape data as 
#'   defined by landmark coordinates. It is assumed that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function may be used to assess
#'   the degree of morphological integration between two separate structures or between two modules defined within 
#'   the same landmark configuration. 
#'   
#'   Two analytical approaches are currently implemented to assess the degree of morphological 
#'   integration:
#'   \enumerate{
#'   \item {method="PLS"} (default) the function estimates the degree of morphological 
#'   integration using two-block partial least squares, or PLS. When used with landmark data, this analysis 
#'   is referred to as singular warps analysis (Bookstein et al. 2003). When {method="PLS"}, the scores along 
#'   the X & Y PLS axes are also returned, as is a plot of PLS scores from Block1 versus Block2 along the first set of
#'   PLS axes. Thin-plate spline deformation grids along these axes are also shown (if data were input as a 3D array).
#'   Note: deformation grids are displayed for each block of landmarks separately.  If the two blocks of landmarks
#'   are derived from a single structure (i.e. a single landmark configuration), one can plot the overall deformation
#'   along PLS1 using the procedure provided in the example below.   
#'   \item {method="RV"} the function estimates the degree of morphological integration using the RV coefficient 
#'   (Klingenberg 2009). Significance testing for both approaches is found by permuting the objects in one data 
#'   matrix relative to those in the other. A histogram of coefficients obtained via resampling is presented, 
#'   with the observed value designated by an arrow in the plot.
#'   }
#'   
#'  If evaluating an a priori hypothesis of modularity within a structure is of interest, one may use the 
#'    average RV coefficient as implemented in the function \code{\link{compare.modular.partitions}}. 
#'
#' @param A1 A matrix (n x [p1 x k]) or 3D array (p1 x k x n) containing GPA-aligned coordinates for the first module
#' @param A2 A matrix (n x [p2 x k]) or 3D array (p2 x k x n) containing GPA-aligned coordinates for the second module 
#' @param method Method to estimate morphological integration; see below for details
#' @param warpgrids A logical value indicating whether deformation grids for shapes along PLS1 should be displayed
#'  (only relevant if data for A1 or A2 [or both] were input as 3D array)
#' @param iter Number of iterations for significance testing
#' @param label An optional vector indicating labels for each specimen that are to be displayed
#' @param verbose A logical value indicating whether the output is basic or verbose ({method="PLS"} only) (see Value below)
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Function returns the a list with the following components: 
#'   \item{value}{The estimate of morphological integration: PLS.corr or RV }
#'   \item{pvalue}{The significance level of the observed signal}
#'   \item{Xscores}{PLS scores for the first block of landmarks ({method="PLS"} only when verbose=TRUE)}
#'   \item{Yscores}{PLS scores for the second block of landmarks ({method="PLS"} only when verbose=TRUE)}
#' @references  Bookstein, F. L., P. Gunz, P. Mitteroecker, H. Prossinger, K. Schaefer, and H. Seidler. 
#'   2003. Cranial integration in Homo: singular warps analysis of the midsagittal plane in ontogeny and 
#'   evolution. J. Hum. Evol. 44:167-187.
#' @references Klingenberg, C. P. 2009. Morphometric integration and modularity in configurations of 
#'   landmarks: tools for evaluating a priori hypotheses. Evol. Develop. 11:405-421.
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#'
#' #Morphological integration using PLS between two modules within a structure
#' morphol.integr(Y.gpa$coords[1:5,,],Y.gpa$coords[6:12,,],method="PLS",iter=99)
#'
#' #Morphological integration using RV between two modules within a structure
#' morphol.integr(Y.gpa$coords[1:5,,],Y.gpa$coords[6:12,,],method="RV",iter=99)
#' 
#' #Deformation plot for case when both blocks are derived from the same landmark configuration
#' res<-morphol.integr(Y.gpa$coords[1:5,,],Y.gpa$coords[6:12,,],method="PLS",iter=99,verbose=TRUE)
#'   ref<-mshape(Y.gpa$coords)   #overall reference
#'   plotRefToTarget(ref,Y.gpa$coords[,,which.min(res$x.scores)],method="TPS") #Min along PLS1
#'   plotRefToTarget(ref,Y.gpa$coords[,,which.max(res$x.scores)],method="TPS") #Max along PLS1
morphol.integr<-function(A1,A2,method=c("PLS","RV"),warpgrids=TRUE,iter=999, label=NULL,verbose=FALSE){
  method <- match.arg(method)
  if(any(is.na(A1))==T){
    stop("Data matrix 1 contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(any(is.na(A2))==T){
    stop("Data matrix 2 contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A1))==3){ 
    if(is.null(dimnames(A1)[[3]])){
      print("No specimen names in data matrix 1. Assuming specimens in same order.")  }
    x<-two.d.array(A1)}
  if (length(dim(A1))==2){ 
    if(is.null(dimnames(A1)[[1]])){
      print("No specimen names in data matrix 1. Assuming specimens in same order.")  }
    x<-A1 }
  if (length(dim(A2))==3){ 
    if(is.null(dimnames(A2)[[3]])){
      print("No specimen names in data matrix 2. Assuming specimens in same order.")  }
    y<-two.d.array(A2)}
  if (length(dim(A2))==2){ 
    if(is.null(dimnames(A2)[[1]])){
      print("No specimen names in data matrix 2. Assuming specimens in same order.")  }
    y<-A2 }
  if(nrow(x)!=nrow(y)){
    stop("Data matrices have different numbers of specimens.")  }
  if(is.null(rownames(x))==FALSE && is.null(rownames(y))==FALSE){
    mtch<-x[is.na( match(rownames(x),rownames(y)))]
    if (length(mtch)>0){stop("Specimen names in data sets are not the same.")  }
  }
  if(is.null(rownames(x))==FALSE && is.null(rownames(y))==FALSE){
    y<-y[rownames(x),]
  }
  XY.vcv<-cov(cbind(x,y))
  S12<-XY.vcv[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]; S21<-t(S12)
  S11<-XY.vcv[1:dim(x)[2],1:dim(x)[2]]
  S22<-XY.vcv[(dim(x)[2]+1):(dim(x)[2]+dim(y)[2]),(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]
  pls<-svd(S12)
  U<-pls$u; V<-pls$v
  XScores<-x%*%U[,1]; YScores<-y%*%V[,1]
  PLS.obs<-cor(XScores,YScores)
  RV.obs<- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22))) 
  integ.obs<-ifelse(method=="PLS",PLS.obs,RV.obs)
  P.val<-1
  integ.val<-rep(0,iter)
  for(i in 1:iter){
    y.r<-y[sample(nrow(y)),]  
    XY.vcv.r<-cov(cbind(x,y.r))
    S12.r<-XY.vcv.r[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2])]; S21.r<-t(S12.r)
    S11.r<-XY.vcv.r[1:dim(x)[2],1:dim(x)[2]]
    S22.r<-XY.vcv.r[(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2]),(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2])]
    pls.r<-svd(S12.r)
    U.r<-pls.r$u; V.r<-pls.r$v
    XScores.r<-x%*%U.r[,1]; YScores.r<-y.r%*%V.r[,1]
    PLS.r<-cor(XScores.r,YScores.r)
    RV.r<- sum(diag(S12.r%*%S21.r))/sqrt(sum(diag(S11.r%*%S11.r))*sum(diag(S22.r%*%S22.r))) 
    integ.r<-ifelse(method=="PLS",PLS.r,RV.r)
    integ.val[i]<-integ.r
    P.val<-ifelse(integ.r>=integ.obs, P.val+1,P.val) 
  }  
  integ.val[iter+1]=integ.obs
  P.val<-P.val/(iter+1)
  if(method=="PLS"){
    xsc<-0; ysc<-0
    if (length(dim(A1))==2 && length(dim(A2))==2){
      plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS Plot",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
      if(length(label!=0)){text(XScores[,1],YScores[,1],label,adj=c(-.7,-.7))}
    }
    if (length(dim(A1))==3){A1.ref<-mshape(A1); xsc<-1
                            pls1.min<-A1[,,which.min(XScores[,1])];pls1.max<-A1[,,which.max(XScores[,1])]}
    if (length(dim(A2))==3){A2.ref<-mshape(A2); ysc<-1
                            pls2.min<-A2[,,which.min(XScores[,1])];pls2.max<-A2[,,which.max(XScores[,1])]}
    if (dim(A1)[2] == 2 ||dim(A2)[2] == 2){
      if(xsc==1 && ysc==1){layout(matrix(c(4,1,5,1,1,1,1,2,1,1,1,1,1,1,1,3),4,4))}
      if(xsc==1 && ysc!=1){layout(matrix(c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,3),4,4))}
      if(xsc!=1 && ysc==1){layout(matrix(c(2,1,3,1,1,1,1,1,1,1,1,1,1,1,1,1),4,4))}
      plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
      if(length(label!=0)){text(XScores[,1],YScores[,1],label,adj=c(-.7,-.7))}
      if(warpgrids==TRUE){
        if (xsc==1){
          tps(A1.ref, pls1.min, 20)
          tps(A1.ref, pls1.max, 20,sz=.7)     
        }
        if (ysc==1){
          tps(A2.ref, pls2.min, 20,sz=.7)
          tps(A2.ref, pls2.max, 20,sz=.7)
        }
      }
      layout(1)
    }

    if (length(dim(A1))==3  && dim(A1)[2] == 3) {
      plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS Plot",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
      if(length(label!=0)){text(XScores[,1],YScores[,1],label,adj=c(-.7,-.7))}
      open3d()
      plot3d(pls1.min, type = "s", col = "gray", main = paste("PLS Block1 negative"),size = 1.25, aspect = FALSE)
      open3d()
      plot3d(pls1.max, type = "s", col = "gray", main = paste("PLS Block1 positive"),size = 1.25, aspect = FALSE)
    }  
    if (length(dim(A2))==3  && dim(A2)[2] == 3){
      open3d()
      plot3d(pls2.min, type = "s", col = "gray", main = paste("PLS Block2 negative"),size = 1.25, aspect = FALSE)
      open3d()
      plot3d(pls2.max, type = "s", col = "gray", main = paste("PLS Block2 positive"),size = 1.25, aspect = FALSE)
    }
    if(verbose==FALSE){ return(list(PLS.corr=integ.obs,pvalue=P.val)) }
    if(verbose==TRUE){ return(list(x.scores=XScores,y.scores=YScores,PLS.corr=integ.obs,pvalue=P.val)) }
  }
  if(method=="RV"){
    hist(integ.val,30,freq=TRUE,col="gray",xlab="RV Coefficient")
    arrows(integ.obs,50,integ.obs,5,length=0.1,lwd=2)
    return(list(RV=integ.obs,pvalue=P.val))
  }
}