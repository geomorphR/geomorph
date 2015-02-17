#' Plot allometric patterns in landmark data
#'
#' Function plots allometry curves for a set of specimens
#'
#' The function performs a regression of shape on size, and generates a plot that describes the 
#' multivariate relationship between size and shape 
#'   derived from landmark data (i.e., allometry). It is assumed that the landmarks have previously been 
#'   aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The abscissa 
#'   of the plot is log(centroid size) while the ordinate represents shape [NOTE: the function takes the 
#'   input size and performed log-transformation automatically by default (logsz = TRUE), as log(Size) should be used]. 
#'   Three complementary approaches can be implemented to visualize allometry: 
#'  \enumerate{
#'   \item {If "method=CAC" (the default) the function calculates the 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   within groups (Mitteroecker et al. 2004). The function also calculates the residual shape component (RSC) for 
#'   the data.}
#'   \item {If "method=RegScore" the function calculates shape scores 
#'   from the regression of shape on size, and plots these versus size (Drake and Klingenberg 2008). 
#'   For a single group, these shape scores are mathematically identical to the CAC (Adams et al. 2013).}
#'   \item {If "method=PredLine" the function calculates predicted values from a regression of shape on size, and 
#'   plots the first principal component of the predicted values versus size as a stylized graphic of the 
#'   allometric trend (Adams and Nistri 2010). }
#'   }
#'   For all methods, both centroid size and allometry scores are returned. Optionally, deformation grids can be 
#'   requested, which display the shape of the smallest and largest specimens relative to the average specimen (using 
#'   'warpgrids=T' or 'warpgrids=F'). 
#'   Finally, if groups are provided, the above approaches are implemented while 
#'   accounting for within-group patterns of covariation (see references for explanation). In this case,
#'   the regression is of the form: shape~size+groups (Note: to examine the interaction term use \code{\link{procD.lm}}).
#'   Specimens from each group are plotted using distinct colors based on the order in which the groups are
#'   found in the dataset, and using R's standard color palette: black, red, green, blue, cyan, magenta,
#'   yellow, and gray. NOTE: to change the colors of the groups, simply substitute a vector of the desired colors for 
#'   each specimen.
#'
#' @param A  matrix (n x [p1 x k]) or 3D array (p1 x k x n) containing GPA-aligned coordinates for the specimens
#' @param sz A vector of size measures for all specimens (centroid size or any other reasonable size measure)
#' @param groups An optional vector containing group labels for each specimen if available 
#' @param method Method for estimating allometric shape components; see below for details
#' @param warpgrids A logical value indicating whether deformation grids for small and large shapes 
#'  should be displayed (note: if groups are provided no TPS grids are shown)
#' @param iter Number of iterations for significance testing
#' @param label An optional vector indicating labels for each specimen that are to be displayed
#' @param mesh A mesh3d object to be warped to represent shape deformation of the directional and fluctuating components
#' of asymmetry if {warpgrids= TRUE} (see \code{\link{warpRefMesh}}).
#' @param logsz A logical value indicating whether log(size) is used 
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @keywords analysis
#' @keywords visualization
#' @export
#' @return Function returns an ANOVA table of statistical results for size: df, SS, MS, Prand.
#' If verbose=TRUE, function returns a list with the following components:
#'  \item{ProcDist.lm}{An ANOVA table as above}
#'  \item{allom.score}{ A matrix of the allometry shape scores}
#'  \item{logSize}{ A matrix of log size}
#'  \item{pred.shape}{A matrix containing the predicted shapes from the regression}
#'  \item{resid.shape}{ The residual shape component (RSC) of the data ("method=CAC" only)}
#' @author Dean Adams
#' @references Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field comes of age: geometric morphometrics 
#'   in the 21st century. Hystrix. 24:7-14. 
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#'   in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#'   transformation of skull shape in St Bernard dogs. Proceedings of the Royal Society B, Biological Sciences 275:71'76.
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#'   Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' @examples
#' data(ratland) 
#' Y.gpa<-gpagen(ratland,PrinAxes=FALSE)    #GPA-alignment
#' 
#' #Using CAC for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="CAC", iter=5)
#'
#' #Using Regression Scores for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="RegScore", iter=5)
#'
#' #Using predicted allometry curve for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="PredLine", iter=5)
plotAllometry<-function(A,sz,groups=NULL,method=c("CAC","RegScore","PredLine"),warpgrids=TRUE,
                        iter=249,label=NULL, mesh=NULL, logsz = TRUE, verbose=FALSE){
  method <- match.arg(method)
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A)) == 3) {
    if (is.null(dimnames(A)[[3]])) {
      print("No specimen names in data matrix. Assuming specimens in same order.")
    }
    y <- two.d.array(A)
  }
  if (length(dim(A)) == 2) {
    if (is.null(dimnames(A)[[1]])) {
      print("No specimen names in data matrix. Assuming specimens in same order.")
    }
    y <- as.matrix(A)
  }
  if(logsz == TRUE){size<-as.matrix(log(sz)); xlab<-"log(Size)" 
    print("Natural log of size is used.")} 
  if(logsz == FALSE) { size<-as.matrix(sz); xlab<-"Size" 
    print("Size has not been log transformed.")}
  n<-nrow(size)
  if(is.null(rownames(size))){
    print("No specimen names in size vector. Assuming specimens in same order.")  }
  if(nrow(y)!=nrow(size)){
    stop("Number of specimens differs from number of values in size vector.")  }
  if(is.null(rownames(y))==FALSE && is.null(rownames(size))==FALSE){
    mtch<-y[is.na( match(rownames(y),rownames(size)))]
    if (length(mtch)>0){stop("Specimen names in data set don't match those in size vector.")  }
  }
  if(is.null(rownames(y))==FALSE && is.null(rownames(size))==FALSE){
    size<-size[rownames(y),]
  }
  if(!is.null(groups)){
    groups<-as.factor(groups)    
    if(is.null(names(groups))){
      print("No specimen names in grouping variable. Assuming specimens in same order.")  }
  }
  if(is.null(rownames(y))==FALSE && is.null(names(groups))==FALSE){
    mtch<-y[is.na( match(rownames(y),names(groups)))]
    if (length(mtch)>0){stop("Specimen names in data set don't match those in grouping variable.")  }
  }
  if(is.null(rownames(y))==FALSE && is.null(names(groups))==FALSE){
    groups<-groups[rownames(y)]
  }
  if(is.null(groups)){lm.res<-procD.lm(y~size,iter=iter)}  
  if(!is.null(groups)){lm.res<-procD.lm(y~size+groups,iter=iter)
                       lm.res2<-procD.lm(y~size*groups,iter=iter)}
  if(is.null(groups)){
    y.mn<-predict(lm(y~1))
    B<-coef(lm(y~size))
    yhat<-predict(lm(y~size))
  }
  if(!is.null(groups)){
    y.mn<-predict(lm(y~groups))
    B<-coef(lm(y~size+groups))
    yhat<-predict(lm(y~size*groups))
    if(lm.res2[3,7]>0.05){
      yhat<-predict(lm(y~size+groups))      
    }
  } 
  asp = NULL
    if(lm.res[1,7]>0.05){ asp <- 1}
  y.cent<-y-y.mn
  a<-(t(y.cent)%*%size)%*%(1/(t(size)%*%size)); a<-a%*%(1/sqrt(t(a)%*%a))
  CAC<-y.cent%*%a  
    resid<-y.cent%*%(diag(dim(y.cent)[2])-a%*%t(a))
  RSC<-prcomp(resid)$x
  Reg.proj<-y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) 
  pred.val<-prcomp(yhat)$x[,1] 
  Ahat<-arrayspecs(yhat,dim(A)[1],dim(A)[2])
  ref<-mshape(A)
  if(method!="CAC"){
    layout(matrix(c(2,1,1,1,1,1,1,1,3),3,3))   
    if(method=="RegScore"){
      plot(size,Reg.proj,xlab=xlab, ylab="Shape (Regression Score)",pch=21,bg="black",cex=1.25, asp=asp)
      if(!is.null(groups)){points(size,Reg.proj,pch=21,bg=groups,cex=1.25)}
      if (length(label!=0)) {
        if(isTRUE(label)){text(size,Reg.proj,seq(1, n),adj=c(-0.7,-0.7)) }
        else{text(size,Reg.proj,label,adj=c(-0.1,-0.1))}
      }
      if(is.null(groups)){
        if(warpgrids==T && dim(A)[2]==2){
          arrows(min(size), (0.7*max(Reg.proj)), min(size), 0, length = 0.1,lwd = 2)
          arrows(max(size), (0.7 * min(Reg.proj)), max(size), 0, length = 0.1,lwd = 2)
        }
      }
    } 
    if(method=="PredLine"){
      plot(size,pred.val,xlab=xlab, ylab="Shape (Predicted)",pch=21,bg="black",cex=1.25, asp=asp)
      if(!is.null(groups)){points(size,pred.val,pch=21,bg=groups,cex=1.25)}
      if (length(label!=0)) {
        if(isTRUE(label)){text(size,pred.val,seq(1, n),adj=c(-0.7,-0.7)) }
        else{text(size,pred.val,label,adj=c(-0.1,-0.1))}
      }
      if(is.null(groups)){
        if(warpgrids==T && dim(A)[2]==2){
          arrows(min(size), (0.7*max(pred.val)), min(size), 0, length = 0.1,lwd = 2)
          arrows(max(size), (0.7 * min(pred.val)), max(size), 0, length = 0.1,lwd = 2)
        }
      }
    }
    if(is.null(groups)){
      if(warpgrids==T && dim(A)[2]==2){
        tps(ref,Ahat[,,which.min(size)],20)
        tps(ref,Ahat[,,which.max(size)],20)
      }
    }
    layout(1)    
  }
  if(method=="CAC"){
    layout(matrix(c(3,1,1,1,1,1,1,1,4,2,2,2,2,2,2,2,2,2),3,6))   
    plot(size,CAC,xlab=xlab, ylab="CAC",pch=21,bg="black",cex=1.25, asp=asp)
    if(is.null(groups)){
      if(warpgrids==T && dim(A)[2]==2){
        arrows(min(size), (0.7*max(CAC)), min(size), 0, length = 0.1,lwd = 2)
        arrows(max(size), (0.7 * min(CAC)), max(size), 0, length = 0.1,lwd = 2)
      }
    }
    if(!is.null(groups)){points(size,CAC,pch=21,bg=groups,cex=1.25)}
    if (length(label!=0)) {
      if(isTRUE(label)){text(size,CAC,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(size,CAC,label,adj=c(-0.1,-0.1))}
    }
    plot(CAC,RSC[,1], xlab="CAC",ylab="RSC 1", pch=21,bg="black",cex=1.25, asp=asp)
    if(!is.null(groups)){points(CAC,RSC[,1],pch=21,bg=groups,cex=1.25)}
    if (length(label!=0)) {
      if(isTRUE(label)){text(CAC,RSC,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(CAC,RSC,label,adj=c(-0.1,-0.1))}
    }
    if(is.null(groups)){
      if(warpgrids==T && dim(A)[2]==2){
        tps(ref,Ahat[,,which.min(size)],20)
        tps(ref,Ahat[,,which.max(size)],20)
      }
    }
    layout(1)
  }
  if(warpgrids==T && dim(A)[2]==3){
    if (is.null(mesh)==TRUE){
      open3d()
      plot3d(Ahat[,,which.min(size)],type="s",col="gray",main="Shape at minimum size",size=1.25,aspect=FALSE)
      open3d()
      plot3d(Ahat[,,which.max(size)],type="s",col="gray",main="Shape at maximum size",size=1.25,aspect=FALSE)
    }
    if(is.null(mesh)==FALSE){
      plotRefToTarget(ref, Ahat[,,which.min(size)], mesh, method = "surface")
      title3d(main="Shape at minimum size")
      plotRefToTarget(ref, Ahat[,,which.max(size)], mesh, method = "surface")
      title3d(main="Shape at maximum size")
    }
  }
  if(verbose==TRUE){ 
    if(method=="CAC"){return(list(allom.score=CAC,resid.shape=RSC,logSize=size,ProcDist.lm=lm.res,pred.shape=Ahat))}
    if(method=="RegScore"){return(list(allom.score=Reg.proj,logSize=size,ProcDist.lm=lm.res,pred.shape=Ahat))}
    if(method=="PredLine"){return(list(allom.score=pred.val,logSize=size,ProcDist.lm=lm.res2,pred.shape=Ahat))}
  }
  if(verbose==FALSE){ return(list(ProcDist.lm=lm.res))}
}