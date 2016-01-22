#' Assessing phylogenetic signal in morphometric data
#'
#' Function calculates the degree of phylogenetic signal from a set of Procrustes-aligned specimens
#'
#' The function estimates the degree of phylogenetic signal present in shape data for a given phylogeny. 
#' It is assumed that the landmarks have previously been aligned 
#'   using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].  
#'   The degree of phylogenetic signal in data is estimated using the multivariate version of the K-statistic 
#'   (Kmult: Adams 2014). This value evaluates the degree of phylogenetic signal
#'   in a dataset relative to what is expected under a Brownian motion model of evolution. For geometric
#'   morphometric data, the approach is a mathematical generalization of the Kappa statistic (Blomberg et al. 
#'   2003) appropriate for highly multivariate data (see Adams 2014).Significance testing 
#'   is found by permuting the shape data among the tips of the phylogeny. Note that this 
#'   method can be quite slow as ancestral states must be estimated for every iteration.
#' 
#' This function can also be used with univariate data (i.e. centroid size) if imported as matrix with rownames
#' giving the taxa names. In this case, the estimate of phylogenetic signal is identical to that found using the 
#' standard kappa statistic (Blomberg et al. 2003).
#' 
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{physignal}}.
#'  The generic function, \code{\link{plot}}, produces a histogram of random K statistics, associated with the resampling procedure.
#'  
#'  \subsection{Notes for geomorph 3.0}{ 
#' Compared to older versions of geomorph, the order of input variables has changed, so that it is consistent with other functions
#' in the program. Additionally, users should note that the function physignal no longer contains
#' multiple methods. Only Kmult is used. Thus, for older scripts method="" should be removed from the function call.
#' }
#' 
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @keywords analysis
#' @author Dean Adams
#' @export
#' @return Function returns a list with the following components: 
#'   \item{phy.signal}{The estimate of phylogenetic signal}
#'   \item{pvalue}{The significance level of the observed signal}
#'   \item{random.K}{Each random K-statistic from random permutations}
#'   \item{permutations}{The number of random permutations used in the resampling procedure}
#'   \item{call}{The matched call}
#' @references Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative 
#' data: behavioral traits are more labile. Evolution, 57:717-745.
#' @references Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and 
#' other high-dimensional multivariate data. Systematic Biology.  63:685-697.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' #Test for phylogenetic signal in shape
#' PS.shape <- physignal(A=Y.gpa$coords,phy=plethspecies$phy,iter=999)
#' summary(PS.shape)
#' plot(PS.shape)
#' 
#' #Test for phylogenetic signal in size
#' PS.size <- physignal(A=Y.gpa$Csize,phy=plethspecies$phy,iter=999)
#' summary(PS.size)
#' plot(PS.size)
physignal<-function(A,phy,iter=999, seed=NULL){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A))==3){ 
    if(is.null(dimnames(A)[[3]])){
      stop("Data array does not include taxa names as dimnames for 3rd dimension.")  }
    x<-two.d.array(A)}
  if (length(dim(A))==2){ 
    if(is.null(rownames(A))){
      stop("Data matrix does not include taxa names as dimnames for rows.")  }
    x<-A }
  if (is.vector(A)== TRUE){ 
    if(is.null(names(A))){
      stop("Data vector does not include taxa names as names.")  }
    x<-as.matrix(A) }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  N<-length(phy$tip.label)
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(x))))) == T) {
    stop("Names do not match between tree and data matrix.") }
#  x<-x[phy$tip.label,]
  if(is.null(dim(x)) == TRUE){ x <- matrix(x, dimnames=list(names(x))) }
  x<-as.matrix(x)
  phy.parts<-phylo.mat(x,phy)
  invC<-phy.parts$invC; D.mat<-phy.parts$D.mat;C = phy.parts$C
  Kmult<-function(x,invC,D.mat){
    ones<-matrix(1,nrow(x),1) 
    a.obs<-colSums(invC)%*%x/sum(invC)  
    distmat<-as.matrix(dist(rbind(as.matrix(x),a.obs))) 
    MSEobs.d<-sum(distmat[(1:N),(N+1)]^2)   
    dist.adj<-as.matrix(dist(rbind((D.mat%*%(x-(ones%*%a.obs))),0))) 
    MSE.d<-sum(dist.adj[(1:N),(N+1)]^2) 
    K.denom<-(sum(diag(C))- N*solve(t(ones)%*%solve(C)%*%ones)) / (N-1)
    K.stat<-(MSEobs.d/MSE.d)/K.denom
    return(K.stat)
  }
  K.obs<-Kmult(x,invC,D.mat)
  ind <- perm.index(nrow(x), iter, seed=seed)
  x.rand <-lapply(1:(iter+1), function(j) x[ind[[j]],])
  K.rand <- sapply(1:(iter+1), function(j) Kmult(as.matrix(x.rand[[j]]), invC,D.mat))
  p.val <- pval(K.rand)
  out <- list(phy.signal=K.obs, pvalue=p.val, random.K = K.rand,
              permutations = iter+1, call=match.call())
  class(out) <- "physignal"
  out
}