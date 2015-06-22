#' Assessing phylogenetic signal in morphometric data
#'
#' Function calculates the degree of phylogenetic signal from a set of Procrustes-aligned specimens
#'
#' The function estimates the degree of phylogenetic signal present in shape data for a given phylogeny. 
#' It is assumed that the landmarks have previously been aligned 
#'   using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].  
#'   Two approaches may be used to quantify phylogenetic signal. First, a multivariate version of the 
#'   K-statistic may be utilized (Kmult: Adams 2014). This value evaluates the degree of phylogenetic signal
#'   in a dataset relative to what is expected under a Brownian motion model of evolution. For geometric
#'   morphometric data, the approach is a mathematical generalization of the Kappa statistic (Blomberg et al. 
#'   2003) appropriate for highly multivariate data (see Adams 2014).
#'     
#'   The second approach estimates phylogenetic signal as the sum of squared changes (SSC) in 
#'   shape along all branches of the phylogeny (Klingenberg and Gidaszewski 2010). Significance testing 
#'   is found by permuting the shape data among the tips of the phylogeny. Note that this 
#'   method can be quite slow as ancestral states must be estimated for every iteration.
#'   
#'   For both approaches a plot of the specimens in tangent space with the phylogeny superimposed 
#'   is included (NOTE: if ancestral states are desired, run \code{\link{plotGMPhyloMorphoSpace}}). 
#'   
#'  The tree must have number of tips equal to number of taxa in the data matrix (e.g. \code{\link[ape]{drop.tip}}).
#' And, tip labels of the tree MUST be exactly the same as the taxa names in the landmark data matrix
#' (check using \code{\link[base]{match}}).
#' 
#' This function can be used with univariate data (i.e. centroid size) if imported as matrix with rownames
#' giving the taxa names.
#'
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param iter Number of iterations for significance testing
#' @param ShowPlot A logical value indicating whether or not the plot should be returned
#' @param method Method for estimating phylogenetic signal (Kmult or SSC)
#' @keywords analysis
#' @author Dean Adams
#' @export
#' @return Function returns a list with the following components: 
#'   \item{phy.signal}{The estimate of phylogenetic signal}
#'   \item{pvalue}{The significance level of the observed signal}
#' @references Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative 
#' data: behavioral traits are more labile. Evolution, 57:717-745.
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic signals 
#'   and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @references Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and 
#' other high-dimensional multivariate data. Systematic Biology.  63:685-697.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' #Test for phylogenetic signal in shape
#' physignal(plethspecies$phy,Y.gpa$coords,method="Kmult",iter=99)
#' 
#' #Test for phylogenetic signal in size
#' Csize <- matrix(Y.gpa$Csize, dimnames=list(names(Y.gpa$Csize))) # make matrix Csize with names
#' physignal(plethspecies$phy,Csize,method="Kmult",iter=99)
physignal<-function(phy,A,iter=999,ShowPlot=TRUE,method=c("Kmult","SSC")){
  method <- match.arg(method)
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A))==3){ 
    if(is.null(dimnames(A)[[3]])){
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
    x<-two.d.array(A)}
  if (length(dim(A))==2){ 
    if(is.null(rownames(A))){
      stop("Data matrix does not include taxa names as dimnames for rows.")  }
    x<-A }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  N<-length(phy$tip.label)
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(x))))) == T) {
    stop("Names do not match between tree and data matrix.") }
  x<-x[phy$tip.label,]
  if(is.null(dim(x)) == TRUE){ x <- matrix(x, dimnames=list(names(x))) }
  if (method=="Kmult"){
    x<-as.matrix(x)
    N<-length(phy$tip.label)
    ones<-array(1,N)
    C<-vcv.phylo(phy)
    C<-C[row.names(x),row.names(x)]
    det.C<-det(C)
    if(det.C>0){invC<-solve(C)}
    if(det.C==0){svd.C<-svd(C)
      Positive <- svd.C$d > max(1e-08 * svd.C$d[1L], 0)
      invC<- svd.C$v[, Positive, drop = FALSE] %*% ((1/svd.C$d[Positive]) *t(svd.C$u[, Positive, drop = FALSE]))}
    eigC <- eigen(C)
    lambda <- zapsmall(eigC$values)
    if(any(lambda == 0)){
      warning("Singular phylogenetic covariance matrix. Proceed with caution")
      lambda = lambda[lambda > 0]
    }
    eigC.vect = eigC$vectors[,1:(length(lambda))]
    D.mat <- solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect))
    Kmult<-function(x,invC,D.mat){
      a.obs<-colSums(invC)%*%x/sum(invC)  
       distmat<-as.matrix(dist(rbind(as.matrix(x),a.obs))) 
      MSEobs.d<-sum(distmat[(1:N),(N+1)]^2)   #sum distances root vs. tips
      dist.adj<-as.matrix(dist(rbind((D.mat%*%(x-(ones%*%a.obs))),0))) 
      MSE.d<-sum(dist.adj[(1:N),(N+1)]^2) #sum distances for transformed data)
      K.denom<-(sum(diag(C))- N*solve(t(ones)%*%solve(C)%*%ones)) / (N-1)
      K.stat<-(MSEobs.d/MSE.d)/K.denom
      return(K.stat)
    }
    K.obs<-Kmult(x,invC,D.mat)
    P.val <- 1
    K.val <- rep(0, iter)
    for (i in 1:iter){
      x.r<-as.matrix(x[sample(nrow(x)),])
      rownames(x.r)<-rownames(x)
      K.rand<-Kmult(x.r,invC,D.mat)
      P.val<-ifelse(K.rand>=K.obs, P.val+1,P.val)     
      K.val[i] <- K.rand
    }   
    P.val <- P.val/(iter + 1)
    K.val[iter + 1] = K.obs
    if(ShowPlot==TRUE && dim(x)[2]>1) {  plotGMPhyloMorphoSpace(phy,A,ancStates=FALSE) }
    return(list(phy.signal=K.obs,pvalue=P.val))     
  }
  if (method=="SSC"){
    anc.states<-NULL
    options(warn=-1)      
    for (i in 1:ncol(x)){
      tmp <- as.vector(fastAnc(phy, x[, i]))
      anc.states<-cbind(anc.states,tmp)   }
    colnames(anc.states)<-NULL
    dist.mat<-as.matrix(dist(rbind(as.matrix(x),as.matrix(anc.states)))^2)   
    SSC.o<-0
    for (i in 1:nrow(phy$edge)){
      SSC.o<-SSC.o+dist.mat[phy$edge[i,1],phy$edge[i,2]]    }
    P.val<-1
    SSC.val<-rep(0,iter)
    for(ii in 1:iter){
      x.r<-x[sample(nrow(x)),] 
      if(is.null(dim(x.r)) == TRUE){ x.r <- matrix(x.r) }
      row.names(x.r)<-row.names(x)
      anc.states.r<-NULL
      options(warn=-1)      
      for (i in 1:ncol(x.r)){
        tmp <- as.vector(fastAnc(phy, x.r[, i]))
        anc.states.r<-cbind(anc.states.r,tmp)   }
      colnames(anc.states.r)<-NULL
      dist.mat.r<-as.matrix(dist(rbind(as.matrix(x.r),as.matrix(anc.states.r)))^2)   
      SSC.r<-0
      for (i in 1:nrow(phy$edge)){
        SSC.r<-SSC.r+dist.mat.r[phy$edge[i,1],phy$edge[i,2]]    }
      P.val<-ifelse(SSC.r<=SSC.o, P.val+1,P.val) 
      SSC.val[ii]<-SSC.r
    }  
    P.val<-P.val/(iter+1)
    SSC.val[iter+1]=SSC.o
    if(ShowPlot==TRUE && dim(x)[2]>1) {plotGMPhyloMorphoSpace(phy,A,ancStates=FALSE)}
    return(list(phy.signal=SSC.o,pvalue=P.val)) 
  }
}